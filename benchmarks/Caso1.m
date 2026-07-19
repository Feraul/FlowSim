classdef Caso1 < SimulacaoBase
    %--------------------------------------------------------------------------
    % Caso439 — Processo de Recarga de Aquífero 2D
    %
    % Modelo físico: Richards não-linear com permeabilidade relativa
    %                de Brooks-Corey
    %
    % Parâmetros:
    %   theta_s = 0.3   (umidade de saturação)
    %   theta_r = 0.0   (umidade residual)
    %   c       = 40000 (parâmetro de forma de Brooks-Corey)
    %   D       = 2.9   (expoente de Brooks-Corey)
    %
    % Condições de contorno:
    %   flag 1 → Dirichlet fixo (bcflag)
    %   flag 2 → Dirichlet variável: h = 65 - z
    %   flag 3 → Dirichlet zero (superfície livre)
    %--------------------------------------------------------------------------

    properties
        Nome   = 'Linearity preserving verification: oblique drain '   % nome legível para log
        TipoID = 1                      % corresponde ao numcase do Start.dat
    end

    methods

        % ── 1. Permeabilidade
        function [env, parms] = configurarPermeabilidade(obj, env, parms, time)
            bedge=env.geometry.bedge;
            centelem=env.geometry.centelem;

            % troca de colunas 1 e 2 de bedge
            bedge(:,[1 2]) = bedge(:,[2 1]);

            x = centelem(:,1);
            y = centelem(:,2);
            n = size(centelem,1);
            alfa  = 0.2;
            theta = atand(alfa);
            phi1 = y - alfa*(x - 0.5) - 0.475;
            phi2 = phi1 - 0.05;

            % mascaras dos 3 ramos originais (condicoes estritas, com "gap" em phi1==0 / phi2==0)
            maskLow1  = phi1 < 0;
            maskHigh  = phi1 > 0 & phi2 < 0;
            maskLow2  = phi2 > 0;
            maskAny   = maskLow1 | maskHigh | maskLow2;
            maskLow   = maskLow1 | maskLow2;

            % matriz de rotacao (constante)
            R = zeros(2);
            R(1,1) = cosd(theta);
            R(1,2) = sind(theta);
            R(2,1) = -R(1,2);
            R(2,2) = R(1,1);
            A = inv(R);

            k_low  = A * [1 0; 0 0.1]  * R;
            k_high = A * [100 0; 0 10] * R;
            k_low_flat  = [k_low(1,1)  k_low(1,2)  k_low(2,1)  k_low(2,2)];
            k_high_flat = [k_high(1,1) k_high(1,2) k_high(2,1) k_high(2,2)];

            kmap = zeros(n,5);   % pre-alocacao explicita (recomendado)
            kmap(maskAny,1)    = find(maskAny);
            kmap(maskAny, 2:5) = maskLow(maskAny) .* k_low_flat + maskHigh(maskAny) .* k_high_flat;

            env.config.perm = kmap;
        end


        % ── 2. Condicoes de contorno de Dirichlet ─────────────────

        function bcattrib = configurarContorno(obj, vertices, flagptr, time, env, pR)

            bcattrib=[];
        end

        % ── 3. Flags de contorno (nflag e nflagface) ──────────────
        % Monta os vetores de flags para vertices e faces de contorno.
        % nflag(:,1)    → tipo de BC (Dirichlet < 200, Neumann >= 200)
        % nflag(:,2)    → valor prescrito de h no vertice
        % nflagface(:,1) → tipo de BC na face
        % nflagface(:,2) → valor prescrito de h na face
        function [nflag, nflagface] = configurarFlags(obj, env, pR, time)

            bedge=env.geometry.bedge;
            coord=env.geometry.coord;
            nflag=50000*ones(size(coord,1),2);
            vert  = bedge(:,1);
            x     = coord(vert,1);
            y     = coord(vert,2);
            delta = 0.2;

            nflag(vert,1) = 101;
            nflag(vert,2) = 2 - x - delta*y;
            nflagface=[];
        end

        % ── 4. Pre-processamento fisico do caso ───────────────────
        % Define os parametros fisicos e as condicoes iniciais:
        %   theta_s, theta_r: parametros de retencao hidrica
        %   h_init: carga hidraulica inicial (h = 65 - z)
        %   h_old:  chute inicial para o metodo iterativo de Picard
        %           h_old = +20 na zona saturada (z < 65)
        %           h_old = -30 na zona nao saturada (z > 65)
        %   dt: passo de tempo = 0.15 dias
        function [parms, env] = preprocessar(obj, env, parms)

        end

        % ── 5. Fontes e pocos ─────────────────────────────────────
        % Delega a definicao de pocos injetores/produtores para
        % a funcao padrao defineWells
        function wells = definirFontes(obj, env, pR)
            wells = [];
        end

        % ── 6a. Modelo de retencao hidrica — Brooks-Corey ─────────
        % Calcula o conteudo volumetrico de agua theta(h):
        %   theta = theta_s                              se h >= 0
        %   theta = theta_r + (theta_s-theta_r)*Se(h)   se h <  0
        % onde Se = c / (c + |h|^2.9)  (saturacao efetiva)
        function theta = calcularTheta(obj, h, parms)
            theta=[];
        end

        % ── 6b. Capacidade hidrica especifica — dtheta/dh ─────────
        % Derivada analitica de theta em relacao a h (Brooks-Corey):
        %   C(h) = dtheta/dh = -(Delta*c*D*|h|^(D-1)*sgn) / (c+|h|^D)^2
        % Usada na montagem da matriz de massa em soil_properties
        % para o metodo de Picard (Richards nao-linear)
        function dthetadh = calcularCapacidade(obj, h, parms)
            dthetadh = [];
        end

        % ── 7. Flowrate boundary ──────────────────────────────────
        % Caso 439 tem contribuicao gravitacional nas faces de contorno
        function flag = temFlowrateBoundary(obj)
            flag = false;
        end

        % Caso 439 nao inverte sinal do flowrate em nenhuma face
        function flowrateZ = ajustarFlowrate(obj, flowrateZ, bedge)
            % nao faz nada — fallback da base seria identico
        end

        % ── 9. Termo temporal de Richards ─────────────────────────
        % Adiciona a matriz de capacidade hidrica e o vetor de acumulacao
        % ao sistema linear — corresponde ao termo dtheta/dt na eq. de Richards
        % Delega para soil_properties que monta o bloco diagonal esparso
        function [M,I] = adicionarTermoTemporal(obj, M, I, parms, flowresultZ, env)
            % Caso1 nao adiciona termo temporal — problema sem retencao hidrica
            % (nao chama soil_properties); M e I retornam inalterados
        end

        % ── 10. Interpolacao de Neumann — LPEW2 ───────────────────
        % Calcula o termo "s" para nos de contorno com condicao de Neumann
        % usado na interpolacao dos pesos LPEW2 (Pre_LPEW_2_vect).
        % Vetorizado: calculado UMA vez fora do loop sobre nos
        %   s(No) = -(1/sum_lambda) * (r1*flux1 + r2*flux2)
        function s = calcularTermoNeumannVet(obj, r, sum_lambda, N, env)
            s=[];
        end

    end

    methods

        % ── 11. Inicializacao antes do loop temporal ───────────────
        % Localiza os elementos de monitoramento (6 pontos de observacao)
        % salva os indices em .mat para reutilizar em simulacoes futuras
        % e inicializa as series temporais de h e theta
        function [parms, extras] = inicializar(obj, env, parms, time)
            extras=false;
        end

        % ── 12. Atualizacao dentro do loop temporal ────────────────
        % A cada passo de tempo:
        %   1. Atualiza h_old com under-relaxation fisicamente motivada
        %      (h_old = +20 na zona saturada, -30 na nao saturada)
        %   2. Chama o pos-processador para salvar VTK
        %   3. Armazena h e theta nos pontos de monitoramento
        function [parms, extras] = atualizarEstado(obj, env, parms, extras, ...
                h, flowrate, time, count)

            extras=false;
        end

        % ── 13. Criterio de parada ────────────────────────────────
        % Caso 439 usa apenas stopcriteria >= 100 (tempo final atingido)
        % sem criterio especial de parada antecipada
        function parar = deveParar(obj, parms, premethod, extras, stopcriteria)
            parar = false;
        end

        % ── 14. Atualiza flags no loop ────────────────────────────
        % Caso 439 tem BC fixas — nao recalcula flags a cada passo
        function flag = precisaAtualizarFlags(obj, time)
            flag = false;
        end

        % ── 15. Finalizacao — graficos ────────────────────────────
        % Plota os resultados apos o loop temporal:
        %   fig 2: perfil de water content na coluna x=20
        %   fig 5-6: series temporais de h e theta na coluna x=11
        function finalizar(obj, env, options)

            arguments
                obj
                env
                options.extras = []
                options.theta_n = []
                options.theta_init_num = []
                options.p = []
                options.flowrate = []
            end

            p        = options.p;

            disp('>> One-phase extrema:');
            max(p)
            min(p)
  
            
            flowrate = options.flowrate;
  
            centelem = env.geometry.centelem;   % CORRIGIDO: geomerty->geometry, centlem->centelem
            bedge    = env.geometry.bedge;
            inedge   = env.geometry.inedge;
            coord    = env.geometry.coord;
            elemarea = env.geometry.elemarea;

            alfa = 0.2;   % ADICIONADO: faltava esta definicao

            % ── solucao analitica nos centroides (vetorizado) ──────────────
            x_c = centelem(:,1);
            y_c = centelem(:,2);
            analytical_sol = 2 - x_c - alfa.*y_c;

            % ── velocidade analitica nas faces (vetorizado) ─────────────────
            nb     = size(bedge,1);
            nFace  = nb + size(inedge,1);

            v1 = [bedge(:,1); inedge(:,1)];
            v2 = [bedge(:,2); inedge(:,2)];

            IJ    = coord(v2,:) - coord(v1,:);        % nFace x 3
            norma = sqrt(sum(IJ.^2, 2));               % nFace x 1

            Rrot     = [0 -1 0; 1 0 0; 0 0 0];         % matriz constante (renomeada p/ nao colidir)
            nij_all  = (Rrot * IJ') ./ norma';         % 3 x nFace

            auxpoint = (coord(v2,:) + coord(v1,:)) * 0.5;
            x_f = auxpoint(:,1);
            y_f = auxpoint(:,2);

            % epsilon so nas faces internas (inedge), igual ao original
            epsVec = zeros(nFace,1);
            epsVec(nb+1:end) = 1e-8;

            phi1 = y_f - alfa*(x_f - 0.5) - 0.475 + epsVec;
            phi2 = phi1 - 0.05;

            maskLow1 = phi1 < 0;
            maskHigh = phi1 > 0 & phi2 < 0;
            maskLow2 = phi2 > 0;
            maskAny  = maskLow1 | maskHigh | maskLow2;
            maskLow  = maskLow1 | maskLow2;

            % ── os dois unicos tensores possiveis (calculados uma vez) ─────
            theta = atand(alfa);
            R2 = zeros(2);
            R2(1,1) = cosd(theta); R2(1,2) = sind(theta);
            R2(2,1) = -R2(1,2);    R2(2,2) = R2(1,1);
            A2 = inv(R2);

            k_low  = A2 * [1 0; 0 0.1]   * R2;
            k_high = A2 * [100 0; 0 10]  * R2;

            KK_low  = [k_low(1,1)  k_low(1,2)  0; k_low(2,1)  k_low(2,2)  0; 0 0 0];
            KK_high = [k_high(1,1) k_high(1,2) 0; k_high(2,1) k_high(2,2) 0; 0 0 0];

            a_low  = -KK_low  * [-1; -alfa; 0];   % vetor 3x1 constante
            a_high = -KK_high * [-1; -alfa; 0];   % vetor 3x1 constante

            % ── produto escalar a.nij para as duas regioes, de uma vez ─────
            dot_low  = a_low'  * nij_all;   % 1 x nFace
            dot_high = a_high' * nij_all;   % 1 x nFace

            F = zeros(nFace,1);
            F(maskLow)  = dot_low(maskLow)';
            F(maskHigh) = dot_high(maskHigh)';

            if ~all(maskAny)
                warning('finalizar:phiGap', ...
                    '%d face(s) cairam exatamente em phi1==0/phi2==0 e ficaram com F=0.', ...
                    sum(~maskAny));
            end

            analytical_vel = F;

            %% O calculo destes erros foram adaptados de Gao an Wu 2010.
            % ── recupera as velocidades numericas (vetorizado) ────────────────
            nb    = size(bedge,1);
            nFace = nb + size(inedge,1);

            v1 = [bedge(:,1); inedge(:,1)];
            v2 = [bedge(:,2); inedge(:,2)];

            norma  = sqrt(sum((coord(v2,:) - coord(v1,:)).^2, 2));
            velnum = flowrate(1:nFace) ./ norma;

            % ── calcula o erro respeito a pressao (vetorizado) ────────────────
            erropressure = sqrt(sum((analytical_sol - p).^2 .* elemarea) / sum(elemarea));

            % ── calcula o erro respeito a velocidade (vetorizado) ─────────────
            Q = zeros(nFace,1);
            Q(1:nb)      = elemarea(bedge(:,3));
            Q(nb+1:end)  = elemarea(inedge(:,3)) + elemarea(inedge(:,4));

            e  = -analytical_vel - velnum;
            er = e.^2;
            errovelocity = sqrt(sum(Q.*er) / sum(Q));

            fprintf('\n>> Erros de verificacao (oblique drain):\n');
            fprintf('   Erro de pressao (L2 ponderado por area): %.6e\n', erropressure);
            fprintf('   Erro de velocidade (L2 ponderado por Q): %.6e\n', errovelocity);
        end
        % ── 16. Escrita de resultados em arquivo ──────────────────
        % Salva os campos h, theta, kmap e os centroides em .txt
        % para pos-processamento externo (ex: Python, MATLAB scripts)
        function escreverResultados(obj, env, varargin)
            % Problemas monofasicos estacionarios nao armazenam resultados em disco
            % (a analise fica so nos graficos/metricas gerados em finalizar)

        end
    end

    methods(Static)

        
    end
end
