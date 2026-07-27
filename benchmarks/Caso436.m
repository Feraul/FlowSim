classdef Caso436 < SimulacaoBase
    %--------------------------------------------------------------------------

    properties
        Nome   = 'Convergence test: Non-Linear Hydraulic Conductivity'   % nome legível para log
        TipoID = 436                     % corresponde ao numcase do Start.dat
    end

    methods

        % ── 1. Permeabilidade relativa — Brooks-Corey ─────────────

        function [env, parms] = configurarPermeabilidade(obj, env, parms, time)
            elem   = env.geometry.elem;
            nelem  = size(elem,1);
            theta  = ones(nelem,1);
            kr     = ones(nelem,1);
            n      = parms.nvg;
            alpha  = parms.alpha;
            h      = parms.h_old;
            neg    = h <= 0;

            theta(neg) = (1 + (-alpha*h(neg)).^n).^(-(n-1)/n);
            kr(neg)    = sqrt(theta(neg)) .* (1 - (1 - theta(neg).^(n/(n-1))).^((n-1)/n)).^2;

            coef                    = env.config.perm(1,1) .* kr;
            env.config.kmap         = obj.iso(env, coef);      % tensor isotropico
            parms.auxperm           = env.config.kmap;         % copia para parms
            env.config.auxkmap      = obj.isoConst(env);        % K saturado (referencia)
            env.geometry.elem(:,5)  = env.utils.idx;
        end
        function flag = precisaAtualizarPermeabilidade(obj)
            flag = true;
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
                coordmid = env.geometry.coord(:, 1:2);
            end

            bcattrib= -3.*time.*coordmid (:,1).*(1-coordmid (:,1)).*coordmid (:,2).*(1-coordmid (:,2))-1;

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
            parms.alpha=0.1844;
            parms.nvg=3;

            parms.h_init=-1*ones(size(env.geometry.elem,1),1);
            parms.h_old=-5*ones(size(env.geometry.elem,1),1);
            parms.dt= (1/64);
        end

        % ── 5. Fontes e pocos ─────────────────────────────────────
        % Delega a definicao de pocos injetores/produtores para
        % a funcao padrao defineWells
        function source_wells = definirFontes(obj, env, parms,time)
            source_wells = defineWells(env, parms,time);
        end

        % ── 6a. Modelo de retencao hidrica — Brooks-Corey ─────────
        % Calcula o conteudo volumetrico de agua theta(h):
        %   theta = theta_s                              se h >= 0
        %   theta = theta_r + (theta_s-theta_r)*Se(h)   se h <  0
        % onde Se = c / (c + |h|^2.9)  (saturacao efetiva)
        function theta = calcularTheta(obj, h, parms)
            alpha=parms.alpha;
            nvg=parms.nvg;

            theta = zeros(size(h));  % inicializa o vetor

            idx_neg = h < 0;          % índices onde h < 0

            % h < 0
            theta(idx_neg) = (1 + (-alpha*h(idx_neg)).^nvg).^(-(nvg-1)/nvg);
        end

        % ── 6b. Capacidade hidrica especifica — dtheta/dh ─────────
        % Derivada analitica de theta em relacao a h (Brooks-Corey):
        %   C(h) = dtheta/dh = -(Delta*c*D*|h|^(D-1)*sgn) / (c+|h|^D)^2
        % Usada na montagem da matriz de massa em soil_properties
        % para o metodo de Picard (Richards nao-linear)
        function dthetadh = calcularCapacidade(obj, h, parms)

            alpha=parms.alpha;
            nvg=parms.nvg;
            % calculo da capacidade da agua

            dthetadh = zeros(size(h));      % inicializa
            idx = (h <= 0);

            dthetadh(idx)  = alpha*(nvg-1) .* ((-alpha*h(idx)).^(nvg-1)) .* ...
                (((-alpha*h(idx)).^nvg + 1).^((1/nvg) - 2));
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
            s=false;
        end
    end

    methods

        % ── 11. Inicializacao antes do loop temporal ───────────────
        % Localiza os elementos de monitoramento (6 pontos de observacao)
        % salva os indices em .mat para reutilizar em simulacoes futuras
        % e inicializa as series temporais de h e theta
        function [parms, extras] = inicializar(obj, env, parms, time)
            extras.hnsum = 0;
            extras.vnsum = 0;

            centelem = env.geometry.centelem;
            x = centelem(:,1); y = centelem(:,2); t = time;
            u_exact0 = -3*t .* x.*(1-x) .* y.*(1-y) - 1;

            [theta0, kmap0] = obj.exataAuxiliares(u_exact0, parms, env);

            extras.exact_solution_storage  = [centelem(:,2), u_exact0];
            extras.theta_storage_analitica = [centelem(:,2), theta0];
            extras.kmap_storage_analitica  = [centelem(:,2), kmap0];
        end

        % ── 12. Atualizacao dentro do loop temporal ────────────────
        % A cada passo de tempo:
        %   1. Atualiza h_old com under-relaxation fisicamente motivada
        %      (h_old = +20 na zona saturada, -30 na nao saturada)
        %   2. Chama o pos-processador para salvar VTK
        %   3. Armazena h e theta nos pontos de monitoramento
        function [parms,extras] = atualizarEstado(obj, env, parms, ...
                h, theta_n, time, count,flowrate,extras)

            centelem=env.geometry.centelem;
            x = centelem(:,1);
            y = centelem(:,2);
            t = time;
            u_exact = -3*t .* x.*(1-x) .* y.*(1-y) - 1;

            p_oldaux1 = (h <= 0);
            p_oldaux2 = (h > 0);
            parms.p_old     = 10*p_oldaux2 - 5*p_oldaux1;

            parms.h_init=u_exact;

            %% calcula theta e salva VTK
            postprocessor(env, count, time, pressure=h, theta_n=theta_n, exact_sol=u_exact);
            %==============================================================

            centelem = env.geometry.centelem;
            bedge    = env.geometry.bedge;
            inedge   = env.geometry.inedge;
            coord    = env.geometry.coord;
            normals  = env.geometry.normals;   % confirme o nome exato desse campo no seu env
            elemarea = env.geometry.elemarea;

            alpha = parms.alpha;
            nvg   = parms.nvg;
            dt    = parms.dt;
            t     = time;

            sizebedge       = size(bedge,1);
            sizebedgeinedge = sizebedge + size(inedge,1);
            %------------------------------------------------------------------
            v1 = bedge(:,1); v2 = bedge(:,2);
            m_bedge = 0.5*(coord(v1,1:2)+coord(v2,1:2));
            areanormal = zeros(sizebedgeinedge,1);
            areanormal(1:sizebedge,1) = norm(coord(v1,1:2)-coord(v2,1:2));
            xbedge = m_bedge(:,1); ybedge = m_bedge(:,2);
            u_exact_bedge = -3*t .* xbedge.*(1-xbedge) .* ybedge.*(1-ybedge) - 1;

            dpdx = -3*t .* (1 - 2*xbedge) .* ybedge.*(1-ybedge);
            dpdy = -3*t .* xbedge.*(1-xbedge) .* (1 - 2*ybedge);

            ktheta = ones(size(u_exact_bedge));
            mask = u_exact_bedge <= 0;
            ktheta(mask) = (1 + (-alpha*u_exact_bedge(mask)).^nvg).^(-(nvg-1)/nvg);

            vx1 = -ktheta .* dpdx;  vy1 = -ktheta .* dpdy;
            vn = zeros(sizebedgeinedge,1);
            vn(1:sizebedge,1) = vx1 .* normals(1:sizebedge,1) + vy1 .* normals(1:sizebedge,2);
            %------------------------------------------------------------------
            v11 = inedge(:,1); v21 = inedge(:,2);
            m_inedge = 0.5*(coord(v11,1:2)+coord(v21,1:2));
            areanormal(sizebedge+1:sizebedgeinedge,1) = norm(coord(v11,1:2)-coord(v21,1:2));
            xinedge = m_inedge(:,1); yinedge = m_inedge(:,2);

            dpdx = -3*t .* (1 - 2*xinedge) .* yinedge.*(1-yinedge);
            dpdy = -3*t .* xinedge.*(1-xinedge) .* (1 - 2*yinedge);
            u_exact_inedge = -3*t .* xinedge.*(1-xinedge) .* yinedge.*(1-yinedge) - 1;

            ktheta = ones(size(u_exact_inedge));
            mask = u_exact_inedge <= 0;
            ktheta(mask) = (1 + (-alpha*u_exact_inedge(mask)).^nvg).^(-(nvg-1)/nvg);

            vx = -ktheta .* dpdx;  vy = -ktheta .* dpdy;
            vn(sizebedge+1:sizebedgeinedge,1) = vx .* normals(sizebedge+1:sizebedgeinedge,1) + vy .* normals(sizebedge+1:sizebedgeinedge,2);
            %--------------------------------------------------------------
            vn_numerico = flowrate ./ areanormal;
            %----------------------------------------------------------------
            wA = elemarea;
            diffp = u_exact - h;
            extras.hnsum = extras.hnsum + dt*(sum(wA .* diffp.^2) / sum(wA));

            Q = zeros(size(inedge,1) + size(bedge,1), 1);
            Q(1:size(bedge,1)) = elemarea(bedge(:,3));
            Q(size(bedge,1)+1:end) = 0.5.*(elemarea(inedge(:,3)) + elemarea(inedge(:,4)));

            e  = vn - vn_numerico;
            extras.vnsum = extras.vnsum + dt*(Q'*e.^2)/sum(Q');
            extras.u_exact=u_exact;
            extras.h=h;

            [theta_a, kmap_a] = obj.exataAuxiliares(u_exact, parms, env);

            extras.exact_solution_storage(:, 2*count-1:2*count)  = [centelem(:,2), u_exact];
            extras.theta_storage_analitica(:, 2*count-1:2*count) = [centelem(:,2), theta_a];
            extras.kmap_storage_analitica(:, 2*count-1:2*count)  = [centelem(:,2), kmap_a];

        end

        function [theta, kmap] = exataAuxiliares(obj, u_exact, parms, env)
            % theta(u_exact) e kmap(u_exact) via van Genuchten-Mualem — mesma formula do configurarPermeabilidade
            alpha = parms.alpha; nvg = parms.nvg;
            theta = ones(size(u_exact));
            kr    = ones(size(u_exact));
            neg   = u_exact <= 0;
            theta(neg) = (1 + (-alpha*u_exact(neg)).^nvg).^(-(nvg-1)/nvg);
            kr(neg)    = sqrt(theta(neg)) .* (1 - (1 - theta(neg).^(nvg/(nvg-1))).^((nvg-1)/nvg)).^2;
            kmap = env.config.perm(1,1) .* kr;
        end

        % ── 13. Criterio de parada ────────────────────────────────
        % Caso 439 usa apenas stopcriteria >= 100 (tempo final atingido)
        % sem criterio especial de parada antecipada
        function parar = deveParar(obj, parms, premethod, stopcriteria)
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
                options.theta_n = []
                options.theta_init_num = []
                options.p = []
                options.flowrate = []
                options.extras = []
            end
            if ~isempty(options.extras)
                %fprintf('Erro L2 pressao (integrado no tempo): %.6e\n', sqrt(options.extras.hnsum));
                %fprintf('Erro L2 velocidade (integrado no tempo): %.6e\n', sqrt(options.extras.vnsum));
                elemarea=env.geometry.elemarea;
                wA = elemarea;                         % pesos por elemento
                diffp = options.extras.u_exact-options.extras.h;                   % diferença analítica - numérica
                hnsum1=(sum(wA .* diffp.^2) / sum(wA) );

                errortotal1=(options.extras.hnsum+options.extras.vnsum); % Analysis of an Euler Implicit-Mixed Finite Element
                % Scheme for Reactive Solute Transport in Porous Media, veja o paragrafo da
                % equacao 5.5, dh=dt.

                errortotal2=(hnsum1+options.extras.vnsum); % artigo:Convergence analysis for a mixed finite element scheme for flow in
                %strictly unsaturated porous media, ultimo paragrafo da equacao 37, dh=dt.

                fprintf('Erro total L1 (integrado no tempo): %.6e\n', errortotal1);
                fprintf('Erro total L2  (integrado no tempo): %.6e\n', errortotal2);
            end

        end

        % ── 16. Escrita de resultados em arquivo ──────────────────
        % Salva os campos h, theta, kmap e os centroides em .txt
        % para pos-processamento externo (ex: Python, MATLAB scripts)
        function escreverResultados(obj, env, h_storage, theta_storage, ...
                kmap_storage, time_storage, centelem, extras)
            filepath  = env.mainpathfolders.path;
            tabfolder = env.mainpathfolders.tabfolder;
            fname = fullfile(filepath, tabfolder);

            writematrix(h_storage,     [fname 'h_steptime3.txt']);
            writematrix(theta_storage, [fname 'WaterContent_steptime3.txt']);
            writematrix(centelem,      [fname 'centrocell3.txt']);
            writematrix(time_storage,  [fname 'time_step3.txt']);
            writematrix(kmap_storage,  [fname 'condhydraulic_steptime3.txt']);

            writematrix(extras.exact_solution_storage,  [fname 'exact_solution_steptime3.txt']);
            writematrix(extras.theta_storage_analitica, [fname 'WaterContent_analitica_steptime3.txt']);
            writematrix(extras.kmap_storage_analitica,  [fname 'condhydraulic_analitica_steptime3.txt']);
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
