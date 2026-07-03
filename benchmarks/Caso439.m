classdef Caso439 < SimulacaoBase
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
        Nome   = 'Processo de Recarga'   % nome legível para log
        TipoID = 439                      % corresponde ao numcase do Start.dat
    end

    methods

        % ── 1. Permeabilidade relativa — Brooks-Corey ─────────────
        % Calcula o tensor de permeabilidade kmap em cada elemento
        % usando a permeabilidade relativa kr(h):
        %   kr = 1                              se h >= 0 (zona saturada)
        %   kr = c / (c + |h|^5)               se h <  0 (zona nao saturada)
        % Resultado armazenado em env.config.kmap e parms.auxperm
        function [env, parms] = configurarPermeabilidade(obj, env, parms, time)
            h   = parms.h_old;
            kr  = ones(env.utils.nelem, 1);
            neg = h < 0;
            kr(neg) = (2.99e6) ./ (2.99e6 + abs(h(neg)).^5);
            coef                   = env.config.perm(1,1) .* kr;
            env.config.kmap        = obj.iso(env, coef);      % tensor isotropico
            parms.auxperm          = env.config.kmap;         % copia para parms
            env.config.auxkmap     = obj.isoConst(env);       % K saturado (referencia)
            env.geometry.elem(:,5) = env.utils.idx;           % id de material = idx elemento
        end

        % ── 2. Condicoes de contorno de Dirichlet ─────────────────
        % Atribui o valor de h nas faces/vertices de contorno
        % de acordo com o flag de cada face:
        %   flag 1 → valor fixo do bcflag (condicao de entrada)
        %   flag 2 → h = 65 - z           (nivel piezometrico linear)
        %   flag 3 → h = 0                (dreno ou superficie livre)
        function bcattrib = configurarContorno(obj, vertices, flagptr, time, env, pR)

            % calcula o ponto medio de cada aresta (coordmid)
            if size(vertices,2) > 1
                % vertices tem 2 colunas → aresta (v1, v2)
                coordmid         = (env.geometry.coord(vertices(:,1),1:2) + ...
                    env.geometry.coord(vertices(:,2),1:2)) ./ 2;
                auxvertices(:,1) = 1:size(flagptr,1);
                vertices         = auxvertices;
            else
                % vertices tem 1 coluna → nos do contorno
                coordmid = env.geometry.coord(1:size(vertices,1), 1:2);
            end

            bcflag   = env.config.bcflag;
            bcattrib = zeros(size(vertices,1), 1);

            for i = 1:size(flagptr,1)
                vert = vertices(i,1);
                switch flagptr(i,1)
                    case 1,  bcattrib(vert,1) = bcflag(1,2);           % entrada fixa
                    case 2,  bcattrib(vert,1) = 65 - coordmid(vert,2); % linear com z
                    case 3,  bcattrib(vert,1) = 0;                      % dreno
                end
            end
        end

        % ── 3. Flags de contorno (nflag e nflagface) ──────────────
        % Monta os vetores de flags para vertices e faces de contorno.
        % nflag(:,1)    → tipo de BC (Dirichlet < 200, Neumann >= 200)
        % nflag(:,2)    → valor prescrito de h no vertice
        % nflagface(:,1) → tipo de BC na face
        % nflagface(:,2) → valor prescrito de h na face
        function [nflag, nflagface] = configurarFlags(obj, env, pR, time)
            nelem_nodes = size(env.geometry.coord, 1);
            nelem_faces = size(env.geometry.bedge, 1);

            % inicializa com valor sentinela (5000 = interior)
            nflag     = 5000 * ones(nelem_nodes, 2);
            nflagface = zeros(nelem_faces, 2);

            % mapeia bcflag → vertices e faces de contorno
            [vertex_idx, face_idx, bcflag_vertex, bcflag_face, ...
                bc_row_vertex, bc_row_face] = obj.prepararIndices(env);

            % preenche flags e valores nos vertices
            nflag(vertex_idx,1) = bcflag_vertex(:,1);
            mmmm                = obj.configurarContorno(vertex_idx, bc_row_vertex', time, env, pR);
            nflag(vertex_idx,2) = mmmm(vertex_idx);

            % preenche flags e valores nas faces
            nflagface(:,1) = bcflag_face(:,1);
            nflagface(:,2) = obj.configurarContorno(...
                [face_idx(:,1) face_idx(:,2)], bc_row_face', time, env, pR);
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
            centelem      = env.geometry.centelem;
            parms.theta_s = 0.3;
            parms.theta_r = 0.0;
            parms.dt      = 1.5e-1;

            % h inicial: coluna d'agua a partir da cota 65
            m            = 65 - centelem(:,2);
            parms.h_init = m .* ones(size(centelem,1), 1);

            % chute inicial para Picard — separado por zona
            h_old         = zeros(size(centelem,1), 1);
            mask          = centelem(:,2) > 65;   % zona nao saturada
            h_old(mask)   = -30;
            h_old(~mask)  =  20;
            parms.h_old   = h_old;
        end

        % ── 5. Fontes e pocos ─────────────────────────────────────
        % Delega a definicao de pocos injetores/produtores para
        % a funcao padrao defineWells
        function wells = definirFontes(obj, env, pR)
            wells = defineWells(env, pR);
        end

        % ── 6a. Modelo de retencao hidrica — Brooks-Corey ─────────
        % Calcula o conteudo volumetrico de agua theta(h):
        %   theta = theta_s                              se h >= 0
        %   theta = theta_r + (theta_s-theta_r)*Se(h)   se h <  0
        % onde Se = c / (c + |h|^2.9)  (saturacao efetiva)
        function theta = calcularTheta(obj, h, parms)
            theta_s = parms.theta_s;
            theta_r = parms.theta_r;
            c       = 40000;

            Se    = ones(size(h));
            theta = theta_s * ones(size(h));   % zona saturada: theta = theta_s

            idx_neg        = h < 0;
            Se(idx_neg)    = c ./ (c + abs(h(idx_neg)).^2.9);
            theta(idx_neg) = theta_r + (theta_s - theta_r) .* Se(idx_neg);
        end

        % ── 6b. Capacidade hidrica especifica — dtheta/dh ─────────
        % Derivada analitica de theta em relacao a h (Brooks-Corey):
        %   C(h) = dtheta/dh = -(Delta*c*D*|h|^(D-1)*sgn) / (c+|h|^D)^2
        % Usada na montagem da matriz de massa em soil_properties
        % para o metodo de Picard (Richards nao-linear)
        function dthetadh = calcularCapacidade(obj, h, parms)
            dthetadh = zeros(size(h));
            idx      = (h < 0);
            theta_s  = parms.theta_s;
            theta_r  = parms.theta_r;
            Delta    = theta_s - theta_r;   % amplitude de variacao de theta
            aps      = abs(h(idx));
            sgn      = sign(h(idx));
            c        = 40000;
            D        = 2.9;
            dthetadh(idx) = -(Delta .* c .* D .* aps.^(D-1) .* sgn) ./ (c + aps.^D).^2;
        end

        % ── 7. Flowrate boundary ──────────────────────────────────
        % Caso 439 tem contribuicao gravitacional nas faces de contorno
        function flag = temFlowrateBoundary(obj)
            flag = true;
        end

        % Caso 439 nao inverte sinal do flowrate em nenhuma face
        function flowrateZ = ajustarFlowrate(obj, flowrateZ, bedge)
            % nao faz nada — fallback da base seria identico
        end

        % ── 8. Permeabilidade na fronteira — Brooks-Corey ─────────
        % Ajusta K11/K22 nas faces de contorno com Dirichlet
        % usando a permeabilidade relativa avaliada em h_contorno:
        %   coef = 35 * kr(h)   com kr de Brooks-Corey
        % O fator 35 é a condutividade hidraulica saturada do solo (cm/dia)
        function [K11, K12, K21, K22] = ajustarKContorno(obj, env, parms, ...
                auxkmap, matid, h_contorno, maskT)
            coef     = zeros(length(maskT),1);
            mask_neg = h_contorno < 0;
            mask_pos = h_contorno >= 0;

            % zona nao saturada: permeabilidade relativa de Brooks-Corey
            coef(mask_neg) = 35 .* (2.99e6 ./ (2.99e6 + abs(h_contorno(mask_neg)).^5));
            % zona saturada: permeabilidade plena
            coef(mask_pos) = 35;

            % aplica apenas nas faces de Dirichlet (maskT = bedge(:,5) < 200)
            K11 = auxkmap(matid,2).*(~maskT) + coef;
            K12 = auxkmap(matid,3);
            K21 = auxkmap(matid,4);
            K22 = auxkmap(matid,5).*(~maskT) + coef;
        end

        % ── 9. Termo temporal de Richards ─────────────────────────
        % Adiciona a matriz de capacidade hidrica e o vetor de acumulacao
        % ao sistema linear — corresponde ao termo dtheta/dt na eq. de Richards
        % Delega para soil_properties que monta o bloco diagonal esparso
        function [M,I] = adicionarTermoTemporal(obj, M, I, parms, flowresultZ, env)
            [M,I] = soil_properties(M, I, parms, flowresultZ, env);
        end

        % ── 10. Interpolacao de Neumann — LPEW2 ───────────────────
        % Calcula o termo "s" para nos de contorno com condicao de Neumann
        % usado na interpolacao dos pesos LPEW2 (Pre_LPEW_2_vect).
        % Vetorizado: calculado UMA vez fora do loop sobre nos
        %   s(No) = -(1/sum_lambda) * (r1*flux1 + r2*flux2)
        function s = calcularTermoNeumannVet(obj, r, sum_lambda, N, env)
            nNodes = size(env.geometry.coord, 1);
            s      = zeros(nNodes, 1);
            bedge  = env.geometry.bedge;
            bcflag = env.config.bcflag;
            nb     = size(bedge, 1);
            ns2    = env.geometry.nsurn2;

            for No = 1:nNodes
                nns   = ns2(No+1) - ns2(No);   % numero de nos vizinhos
                face1 = N(No,1);               % primeira face de contorno
                face2 = N(No,nns);             % ultima face de contorno

                % no nao toca o contorno — pula
                if face1 > nb || face2 > nb, continue; end

                % verifica se o no e de contorno (bedge(:,1) == No)
                MM  = bedge(:,1) == No;
                MMM = find(MM, 1);
                if isempty(MMM), continue; end

                % verifica se tem condicao de Neumann (flag > 200)
                if ~(200 < bedge(MMM,4)), continue; end

                % localiza o fluxo prescrito em bcflag
                a  = bcflag(:,1) == bedge(face1,5);
                b  = bcflag(:,1) == bedge(face2,5);
                s1 = find(a, 1);
                s2 = find(b, 1);
                if isempty(s1) || isempty(s2), continue; end

                % formula eq. 26 do artigo LPEW2 (Agelas et al., 2010)
                s(No) = -(1/sum_lambda(No)) * ...
                    (r(No,1)*bcflag(s1,2) + r(No,2)*bcflag(s2,2));
            end
        end

    end

    methods

        % ── 11. Inicializacao antes do loop temporal ───────────────
        % Localiza os elementos de monitoramento (6 pontos de observacao)
        % salva os indices em .mat para reutilizar em simulacoes futuras
        % e inicializa as series temporais de h e theta
        function [parms, extras] = inicializar(obj, env, parms, time)
            centelem = env.geometry.centelem;

            % pontos de monitoramento — coluna x=11 (pontos 1,2,3)
            extras.centro1 = find((105<centelem(:,2) & centelem(:,2)<110) & ...
                (centelem(:,1)>10  & centelem(:,1)<15));
            extras.centro2 = find((130<centelem(:,2) & centelem(:,2)<135) & ...
                (centelem(:,1)>10  & centelem(:,1)<15));
            extras.centro3 = find((185<centelem(:,2) & centelem(:,2)<190) & ...
                (centelem(:,1)>10  & centelem(:,1)<15));

            % pontos de monitoramento — coluna x=161 (pontos 4,5,6)
            extras.centro4 = find((80<centelem(:,2)  & centelem(:,2)<85)  & ...
                (centelem(:,1)>160 & centelem(:,1)<165));
            extras.centro5 = find((115<centelem(:,2) & centelem(:,2)<120) & ...
                (centelem(:,1)>160 & centelem(:,1)<165));
            extras.centro6 = find((155<centelem(:,2) & centelem(:,2)<160) & ...
                (centelem(:,1)>160 & centelem(:,1)<165));

            % persiste indices para proxima simulacao (mesma malha) — cache em data/
            cacheDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data');
            if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end
            save(fullfile(cacheDir, 'indices_elementos_quadrilateral.mat'), '-struct', 'extras');

            % series temporais — linha t=0
            theta_init_num = thetafunction(parms.h_init, parms, env);
            extras.h_time1 = [time, parms.h_init(extras.centro1), theta_init_num(extras.centro1)];
            extras.h_time2 = [time, parms.h_init(extras.centro2), theta_init_num(extras.centro2)];
            extras.h_time3 = [time, parms.h_init(extras.centro3), theta_init_num(extras.centro3)];
            extras.h_time4 = [time, parms.h_init(extras.centro4), theta_init_num(extras.centro4)];
            extras.h_time5 = [time, parms.h_init(extras.centro5), theta_init_num(extras.centro5)];
            extras.h_time6 = [time, parms.h_init(extras.centro6), theta_init_num(extras.centro6)];
        end

        % ── 12. Atualizacao dentro do loop temporal ────────────────
        % A cada passo de tempo:
        %   1. Atualiza h_old com under-relaxation fisicamente motivada
        %      (h_old = +20 na zona saturada, -30 na nao saturada)
        %   2. Chama o pos-processador para salvar VTK
        %   3. Armazena h e theta nos pontos de monitoramento
        function [parms, extras] = atualizarEstado(obj, env, parms, extras, ...
                h, flowrate, time, count)

            % chute inicial para proxima iteracao de Picard
            p_old        = zeros(size(h));
            p_old(h >= 0) =  20;   % zona saturada
            p_old(h <  0) = -30;   % zona nao saturada
            parms.h_old  = p_old;
            parms.h_init = h;      % atualiza condicao inicial para proximo dt

            % flowresultZ: fluxo gravitacional acumulado por elemento
            if strcmp(env.config.pmethod,'tpfa')
                flowresultZ = env.premethod.TPFA.flowresultZ;
            else
                flowresultZ = env.premethod.MPFAD.flowresultZ;
            end

            % calcula theta e salva VTK
            theta_n = thetafunction(h, parms, env);
            postprocessor(h, theta_n, 0*flowresultZ, time, env, count, parms);

            % series temporais nos 6 pontos de monitoramento
            extras.h_time1(count,1:3) = [time, h(extras.centro1), theta_n(extras.centro1)];
            extras.h_time2(count,1:3) = [time, h(extras.centro2), theta_n(extras.centro2)];
            extras.h_time3(count,1:3) = [time, h(extras.centro3), theta_n(extras.centro3)];
            extras.h_time4(count,1:3) = [time, h(extras.centro4), theta_n(extras.centro4)];
            extras.h_time5(count,1:3) = [time, h(extras.centro5), theta_n(extras.centro5)];
            extras.h_time6(count,1:3) = [time, h(extras.centro6), theta_n(extras.centro6)];
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
        function finalizar(obj, env, extras, theta_n)
            centelem = env.geometry.centelem;

            figure(2);
            idx    = centelem(:,2)<200 & centelem(:,1)>20 & centelem(:,1)<25;
            centro = find(idx);
            plot(theta_n(centro), centelem(centro,2));
            xlabel('Water Content'); ylabel('Z'); grid;

            figure(5);
            plot(extras.h_time1(:,1), extras.h_time1(:,2)); hold on;
            plot(extras.h_time2(:,1), extras.h_time2(:,2)); hold on;
            plot(extras.h_time3(:,1), extras.h_time3(:,2));
            xlabel('Time'); ylabel('Water Pressure'); title('x=11'); grid;

            figure(6);
            plot(extras.h_time1(:,1), extras.h_time1(:,3)); hold on;
            plot(extras.h_time2(:,1), extras.h_time2(:,3)); hold on;
            plot(extras.h_time3(:,1), extras.h_time3(:,3));
            xlabel('Time'); ylabel('Water content'); title('x=11'); grid;
        end

        % ── 16. Escrita de resultados em arquivo ──────────────────
        % Salva os campos h, theta, kmap e os centroides em .txt
        % para pos-processamento externo (ex: Python, MATLAB scripts)
        function escreverResultados(obj, env, h_storage, theta_storage, ...
                kmap_storage, time_storage)
            centelem  = env.geometry.centelem;
            filepath  = env.mainpathfolders.path;
            tabfolder = env.mainpathfolders.tabfolder;
            fname     = [filepath '\' tabfolder '\'];

            writematrix(h_storage,     [fname 'h_steptime3.txt']);
            writematrix(theta_storage, [fname 'WaterContent_steptime3.txt']);
            writematrix(centelem,      [fname 'centrocell3.txt']);
            writematrix(time_storage,  [fname 'time_step3.txt']);
            writematrix(kmap_storage,  [fname 'condhydraulic_steptime3.txt']);
        end

    end
end