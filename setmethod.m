%--------------------------------------------------------------------------
% setmethod — Gerenciador do tipo de simulacao
%
% Seleciona e executa o solver adequado de acordo com env.config.phasekey:
%
%   phasekey = 1 → Escoamento monofasico (One-phase flow)
%   phasekey = 4 → Carga hidraulica (Groundwater / Hydraulic head)
%   phasekey = 5 → Transporte de contaminante com carga hidraulica
%   phasekey = 6 → Equacao de Richards (zona saturada e nao saturada)
%
% Observacoes:
%   - env.premethod e env.preGravity ja estao populados (via preprocessmethod)
%   - parms contem todos os parametros fisicos (h_init, dt, theta_s...)
%   - preprocessmethod NAO e chamado aqui — ja foi executado no main
%
% Entradas:
%   source_wells — pocos injetores/produtores (definidos por sim.definirFontes)
%   keywrite     — flag de escrita de VTK ('i' = escreve, outros = nao)
%   invh         — flag de inversao de h para pos-processamento
%   env          — estrutura global (geometria, config, benchmark, metodo...)
%   parms        — parametros fisicos do caso
%--------------------------------------------------------------------------
function setmethod(source_wells, keywrite, invh, env, parms)

switch env.config.phasekey

    %% ── Caso 1: Escoamento monofasico (One-phase flow) ──────────
    % Resolve o sistema linear estacionario A*p = b uma unica vez
    % e exibe os extremos de pressao no console
    case 1
        [pressure, flowrate] = env.metodo.resolver(env, env.premethod, ...
            parms, parms.h_init, parms.dt, 0);

        disp('>> One-phase extrema:');
        max(pressure)
        min(pressure)

        % salva VTK com campo de pressao e flowrate
        postprocessor(pressure, flowrate, 0, 0, env, 1, parms);

    %% ── Caso 4: Carga hidraulica (Groundwater) ───────────────────
    case 4

        if ismember(env.config.numcase, [336,334,335,337,338,339,340,341,347,341.1])
            %% Estado estacionario
            % Resolve o sistema linear A*h = b uma unica vez (sem loop temporal)
            % Usado para aquiferos confinados com BC fixas no tempo
            [pressure, flowrate] = env.metodo.resolver(env, env.premethod, ...
                parms, parms.h_init, parms.dt, 0);

            % pos-processamento: salva VTK com campo de carga hidraulica
            postprocessor(pressure, flowrate, 0, 1, 1, ...
                parms.overedgecoord, 1, keywrite, invh, env.config.normk, 0);

            % caso 333: exporta tambem para plotandwrite (comparacao com analitico)
            if env.config.numcase == 333
                plotandwrite(0, 0, pressure, 0, 0, 0, 0, 0, parms.overedgecoord);
            end

            disp('------------------------------------------------');
            disp('>> Global hydraulic head extrema values:');
            max_conval = max(pressure)
            min_conval = min(pressure)

        else
            %% Estado transiente
            % Loop temporal completo para aquiferos com BC ou K variaveis no tempo
            % ex: casos 342, 343 (carga hidraulica transiente)
            hydraulic(source_wells, parms);
        end

    %% ── Caso 5: Transporte de contaminante com carga hidraulica ──
    % Resolve o sistema acoplado IMHEC (Implicit Method for
    % Hydraulic head and Contaminant transport):
    %   1. Define elementos de injec/producao e saturacao nos contornos
    %   2. Identifica vertices/faces com saturacao conhecida (flagknownedge)
    %   3. Pre-processa a equacao de saturacao (preSaturation)
    %   4. Chama IMHEC que resolve h e C de forma acoplada
    case 5
        % define pocos e saturacao nos contornos
        [injecelem, producelem, satinbound, Con, wellsc] = ...
            wellsparameter(parms.wellsc, parms.Con, parms.klb);

        % identifica vertices e faces com saturacao/concentracao conhecida
        [~, ~, ~, flagknownedge] = getsatandflag(satinbound, injecelem, ...
            Con, env.config.nflag, env.config.nflagface, 0);

        % pre-processamento da equacao de saturacao (triangulacao, pesos...)
        [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,isonbound] = ...
            preSaturation(flagknownedge, injecelem, producelem);

        % resolve sistema acoplado carga hidraulica + concentracao
        IMHEC(source_wells, keywrite, invh, env, parms);

    %% ── Caso 6: Equacao de Richards (zona saturada e nao saturada) ─
    % Loop temporal nao-linear completo:
    %   - Linearizacao por Picard a cada passo de tempo
    %   - Atualizacao de K(h), theta(h) e capacidade hidrica C(h)
    %   - Criterio de parada por tempo final ou criterio especial do benchmark
    case 6
        hydraulic_RE(env, parms, source_wells);

    otherwise
        error('setmethod: phasekey %d nao reconhecido.', env.config.phasekey);

end
end