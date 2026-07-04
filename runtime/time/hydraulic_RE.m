%--------------------------------------------------------------------------
% hydraulic_RE — Loop temporal para a Equacao de Richards
%
% Resolve o problema transiente nao-linear de fluxo em meio nao saturado:
%
%   d(theta(h))/dt - div(K(h) grad(h)) = f
%
% Metodo:
%   - Discretizacao temporal implicita (Euler implicito)
%   - Linearizacao por iteracoes de Picard (dentro de ferncodes_solver)
%   - Atualizacao de kmap a cada passo de tempo (Richards nao-linear)
%
% Entradas:
%   env          — estrutura global (geometria, config, premethod, benchmark, metodo)
%   parms        — parametros fisicos (h_init, h_old, dt, theta_s, theta_r...)
%   source_wells — pocos injetores/produtores
%--------------------------------------------------------------------------
function hydraulic_RE(env, parms, source_wells)

centelem  = env.geometry.centelem;
dt        = parms.dt;
finaltime = env.config.totaltime;

%% ── Inicializacao ────────────────────────────────────────────────

time         = 0;
stopcriteria = 0;    % percentual do tempo total concluido (0 a 100)
count        = 1;    % contador de passos de tempo
h            = parms.h_init;

% theta no instante inicial (t=0) — usado para calculo de erros e MBE
theta_init = thetafunction(h, parms, env);

% storage de campos ao longo do tempo — coluna 1: coordenada z, coluna 2+: valores
% crescimento por colunas (mais eficiente que por linhas para acesso sequencial)
h_storage     = [centelem(:,2) h];
theta_storage = [centelem(:,2) theta_init];
kmap_storage  = env.config.auxkmap(:,2);
time_storage  = 0;
dtaux         = dt;

% ── extrai referencias fora do loop ──────────────────────────────
% evita lookup em env a cada iteracao (overhead de struct field access)
benchmark = env.benchmark;
metodo    = env.metodo;

% ── flag pre-calculado fora do loop ──────────────────────────────
% precisaAtualizarFlags e constante durante o loop para a maioria dos casos
% (so caso 436 retorna true — BC dependente do tempo)
flag_atualizaFlags = benchmark.precisaAtualizarFlags(0);

% ── inicializacao especifica do benchmark ─────────────────────────
% ex: caso 439 calcula pontos de monitoramento e series temporais t=0
[parms, extras] = benchmark.inicializar(env, parms, time);

tic
%% ── Loop temporal principal ──────────────────────────────────────
while stopcriteria < 100

    %% ── 1. Resolve o sistema linear em h_new ────────────────────
    % ferncodes_solver:
    %   a) monta a matriz global A (benchmark.montarSistema)
    %   b) adiciona termos de fonte e poco
    %   c) itera Picard ate convergencia (ferncodes_iterpicard)
    %   d) calcula flowrate com h convergido (metodo.calcularFlowrate)
    [h_new, flowrate, flowresult, flowratedif, faceaux, parms, env] = ...
        ferncodes_solver(env, parms, dtaux, source_wells, time);

    %% ── 2. Avanca o tempo ────────────────────────────────────────
    time         = time + dt;
    stopcriteria = time*100/finaltime;   % percentual concluido
    disp([num2str(stopcriteria) '% concluded']);
    count = count + 1;
    h     = h_new;
    dtaux = dt;

    %% ── 3. Calcula theta e armazena campos ───────────────────────
    % theta calculado UMA vez aqui e passado para atualizarEstado
    % evita recalculo dentro do benchmark (custo O(nelem))
    theta_n = thetafunction(h, parms, env);

    % armazena h, theta e kmap no passo atual
    % formato: [z, valor_t1, valor_t2, ...] — uma linha por elemento
    h_storage(:, 2*count-1:2*count)     = [centelem(:,2) h];
    theta_storage(:, 2*count-1:2*count) = [centelem(:,2) theta_n];
    kmap_storage(:, 2*count-1:2*count)  = [centelem(:,2) parms.auxperm(:,2)];
    time_storage(count,1)               = time;

    %% ── 4. Logica especifica do benchmark ────────────────────────
    % ex: caso 439 → atualiza h_old, chama postprocessor, armazena h_time
    % ex: caso 436 → calcula erro L2, atualiza fonte, atualiza h_init exato
    [parms, extras] = benchmark.atualizarEstado(env, parms, extras, ...
        h, flowrate, time, count);

    %% ── 5. Atualiza kmap e premethod para o proximo passo ────────
    % Richards nao-linear: K(h) muda a cada passo de tempo
    % PLUG_kfunction   → recalcula kmap com novo h (via benchmark.configurarPermeabilidade)
    % atualizarPremethod → recalcula Kde, Ded, Kn, Kt, pesos LPEW2 com novo kmap
    [env,parms] = PLUG_kfunction(env, parms, time);
    [env] = metodo.atualizarPremethod(env, parms);

    %% ── 6. Atualiza flags se necessario ──────────────────────────
    % Apenas caso 436 (BC dependente do tempo) retorna flag=true
    % Para caso 439: flag_atualizaFlags = false → bloco nao executa
    if flag_atualizaFlags
        [env] = ferncodes_calflag(env, parms, time);
    end

    %% ── 7. Criterio de parada especial do benchmark ──────────────
    % caso 431: para quando a frente de umidade atinge a altura alvo
    % caso 439: sempre false (para apenas pelo stopcriteria >= 100)
    if benchmark.deveParar(parms, env.premethod, extras, stopcriteria)
        break
    end

    %% ── 8. Log periodico ─────────────────────────────────────────
    % exibe progresso a cada 10 passos para nao poluir o console
    if mod(count, 10) == 0
        fprintf('%.1f%% concluded (t=%.4f)\n', stopcriteria, time);
    end

end % while
toc

%% ── Finalizacao — graficos e erros ──────────────────────────────
% delega ao benchmark: plots de h(t), theta(z), frentes de umidade, erros L2
% ex: caso 439 → figures 2,5,6,7,8
% ex: caso 437 → plots de massa total, erro relativo, MBE
env.benchmark.finalizar(env, extras, theta_n,theta_init);

%% ── Escrita de resultados em arquivo ─────────────────────────────
% trunca o storage ao tamanho real (count passos realizados)
h_storage     = h_storage(:,     1:2*count);
theta_storage = theta_storage(:, 1:2*count);
kmap_storage  = kmap_storage(:,  1:2*count);
time_storage  = time_storage(    1:count);

% delega ao benchmark: cada caso escreve os arquivos que precisa
% ex: caso 439 → h_steptime3.txt, WaterContent_steptime3.txt, ...
env.benchmark.escreverResultados(env, h_storage, theta_storage, ...
    kmap_storage, time_storage);

disp('------------------------------------------------');
disp('>> Global Hydraulic head extrema values [hmax hmin]:');
fprintf('hmax = %.6f | hmin = %.6f\n', max(h), min(h));
end