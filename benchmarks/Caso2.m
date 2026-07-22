classdef Caso2 < SimulacaoBase
    

    properties
        Nome   = 'Convergence test: heterogeneous rotating anisotropy '   % Yuan and Sheng
        TipoID = 2                      % corresponde ao numcase do Start.dat
    end

    methods

        % ── 1. Permeabilidade
        function [env, parms] = configurarPermeabilidade(obj, env, parms, time)
            centelem=env.geometry.centelem;
            elem=env.geometry.elem;
            x = centelem(:,1);
            y = centelem(:,2);
            n = size(elem,1);

            a = 1 + 2*x.^2 + y.^2;
            b = 1 + x.^2 + 2*y.^2;

            theta = 5*pi/12;
            c = cos(theta);
            s = sin(theta);

            k11 = c^2*a + s^2*b;
            k12 = c*s*(a - b);       % = k21 (tensor simetrico)
            k22 = s^2*a + c^2*b;

            kmap = [(1:n)', k11, k12, k12, k22];

            env.config.perm = kmap;
        end
        % K e constante (nao depende de h) -- mas o esquema numerico pode ser
        % nao-linear mesmo assim, entao ainda passamos por Picard (FPI) sem
        % recalcular a permeabilidade a cada iteracao
        function flag = precisaAtualizarPermeabilidade(obj)
            flag = false;
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

            nflag(vert,1)=101;
            nflag(vert,2)=sin(pi*x).*sin(pi*y);
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
        function source_wells = definirFontes(obj, env, pR)
            centelem=env.geometry.centelem;
            elemarea=env.geometry.elemarea;
            x = centelem(:,1);
            y = centelem(:,2);

            t1 = 2.46711*(x.^2-y.^2).*cos(pi*x).*cos(pi*y) - pi^2*(1+1.0669*x.^2+1.93289*y.^2).*sin(pi*x).*sin(pi*y) + ...
                1.57061*x.*sin(pi*x).*cos(pi*y) + 6.70353*x.*cos(pi*x).*sin(pi*y);

            t2 = 2.46711*(x.^2-y.^2).*cos(pi*x).*cos(pi*y) - pi^2*(1+1.93289*x.^2+1.0669*y.^2).*sin(pi*x).*sin(pi*y) + ...
                6.70353*y.*sin(pi*x).*cos(pi*y) - 1.57061*y.*cos(pi*x).*sin(pi*y);

            source_wells.source = -(t1+t2) .* elemarea;
            source_wells.wells=[];
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
            centelem= env.geometry.centelem;
            x=centelem(:,1);
            y=centelem(:,2);
            
            analytical_sol= sin(pi*x) .* sin(pi*y);
            

            % ── velocidade analitica nas faces (vetorizado) ─────────────────
            nb    = size(bedge,1);
            nFace = nb + size(inedge,1);

            v1 = [bedge(:,1); inedge(:,1)];
            v2 = [bedge(:,2); inedge(:,2)];

            IJ    = coord(v2,:) - coord(v1,:);
            norma = sqrt(sum(IJ.^2, 2));

            Rrot    = [0 -1 0; 1 0 0; 0 0 0];
            nij_all = (Rrot * IJ') ./ norma';        % 3 x nFace

            auxpoint = (coord(v2,:) + coord(v1,:)) * 0.5;
            x = auxpoint(:,1);
            y = auxpoint(:,2);

            t0  = pi*(1+1.0670*x.^2+1.9330*y.^2).*cos(pi*x).*sin(pi*y) + pi*(x.^2-y.^2)*0.25.*sin(pi*x).*cos(pi*y);
            t01 = pi*(1+1.9330*x.^2+1.0670*y.^2).*sin(pi*x).*cos(pi*y) + pi*(x.^2-y.^2)*0.25.*cos(pi*x).*sin(pi*y);

            a_all = [-t0, -t01, zeros(nFace,1)];     % nFace x 3

            F = sum(a_all .* nij_all', 2);           % dot(a,nij) linha a linha

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

            fprintf('\n>> Erros de verificacao :\n');
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
