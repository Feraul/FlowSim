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
            nelem_nodes = size(env.geometry.coord, 1);
            nelem_faces = size(env.geometry.bedge, 1);

            bedge=env.geometry.bedge;
            bcflag=env.config.bcflag;

            % inicializa com valor sentinela (5000 = interior)
            nflag     = 5000 * ones(nelem_nodes, 2);
            nflagface = zeros(nelem_faces, 2);

            [tf, loc] = ismember(bedge(:,4), bcflag(:,1));

            % loc(k) = indice da linha de bcflag que corresponde a bedge(k,4)
            % (assume que todo bedge(:,4) tem correspondencia em bcflag(:,1))

            vert = bedge(:,1);
            nflag(vert,1) = bcflag(loc,1);
            nflag(vert,2) = bcflag(loc,2);
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
            const=env.config.perm(1,1);

            % zona nao saturada: permeabilidade relativa de Brooks-Corey
            coef(mask_neg) = const .* (2.99e6 ./ (2.99e6 + abs(h_contorno(mask_neg)).^5));
            % zona saturada: permeabilidade plena
            coef(mask_pos) = const;

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
        function finalizar(obj, env, extras, theta_n,theta_init_num)
            
        end

        % ── 16. Escrita de resultados em arquivo ──────────────────
        % Salva os campos h, theta, kmap e os centroides em .txt
        % para pos-processamento externo (ex: Python, MATLAB scripts)
        function escreverResultados(obj, env, h_storage, theta_storage, ...
                kmap_storage, time_storage)
            centelem  = env.geometry.centelem;
            filepath  = env.mainpathfolders.path;
            tabfolder = env.mainpathfolders.tabfolder;
            fname = fullfile(filepath, tabfolder);

            writematrix(h_storage,     [fname 'h_steptime3.txt']);
            writematrix(theta_storage, [fname 'WaterContent_steptime3.txt']);
            writematrix(centelem,      [fname 'centrocell3.txt']);
            writematrix(time_storage,  [fname 'time_step3.txt']);
            writematrix(kmap_storage,  [fname 'condhydraulic_steptime3.txt']);

        end
    end

    methods(Static)

        function idx = elemento_no_ponto(elem, coord, px, py)
            % Retorna o índice do elemento que contém o ponto (px,py)
            % elem  : matriz de conectividade (colunas 1:4 = nós, coluna 5 = material)
            %         triângulos têm elem(:,4) == 0
            % coord : coordenadas dos nós
            % px,py : coordenadas do ponto de interesse (escalares)

            n1 = elem(:,1);
            n2 = elem(:,2);
            n3 = elem(:,3);
            n4 = elem(:,4);

            isTri = (n4 == 0);
            n4(isTri) = n1(isTri);   % fecha o polígono no triângulo com aresta degenerada

            xv = [coord(n1,1), coord(n2,1), coord(n3,1), coord(n4,1)];
            yv = [coord(n1,2), coord(n2,2), coord(n3,2), coord(n4,2)];

            % vértices "seguintes" (wrap-around: 1->2->3->4->1)
            xv2 = xv(:, [2 3 4 1]);
            yv2 = yv(:, [2 3 4 1]);

            % algoritmo de ray casting (par-ímpar), vetorizado nas 4 arestas
            cond1  = (yv > py) ~= (yv2 > py);
            denom  = yv2 - yv;
            denom(denom == 0) = eps;              % evita divisão por zero (aresta horizontal/degenerada)
            xCross = (xv2 - xv) .* (py - yv) ./ denom + xv;
            cond2  = px < xCross;

            crossings  = cond1 & cond2;
            dentro     = mod(sum(crossings, 2), 2) == 1;

            idx = find(dentro, 1);   % mantém o comportamento original: primeiro elemento encontrado
            if isempty(idx)
                idx = [];
            end
        end


        function idx = elementos_centroide_na_caixa(elem, coord, ylim, xlim)
            % Retorna os indices dos elementos cujo centroide esta dentro da caixa
            % definida por ylim = [ymin ymax] e xlim = [xmin xmax]
            % elem  : matriz de conectividade (colunas 1:4 = nos, 0 = sem no / triangulo)
            % coord : coordenadas dos nos

            n1 = elem(:,1);
            n2 = elem(:,2);
            n3 = elem(:,3);
            n4 = elem(:,4);

            isTri = (n4 == 0);

            % indice seguro para indexacao (evita indice 0); contribuicao sera zerada
            n4safe = n4;
            n4safe(isTri) = n1(isTri);

            x4 = coord(n4safe,1);
            y4 = coord(n4safe,2);
            x4(isTri) = 0;
            y4(isTri) = 0;

            nNos       = 4 * ones(size(elem,1),1);
            nNos(isTri) = 3;

            cx = (coord(n1,1) + coord(n2,1) + coord(n3,1) + x4) ./ nNos;
            cy = (coord(n1,2) + coord(n2,2) + coord(n3,2) + y4) ./ nNos;

            mask = (cy > ylim(1) & cy < ylim(2)) & (cx > xlim(1) & cx < xlim(2));
            idx  = find(mask);
        end
    end
end
