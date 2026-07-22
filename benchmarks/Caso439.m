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
        function [parms] = inicializar(obj, env, parms, time)            
        end

        % ── 12. Atualizacao dentro do loop temporal ────────────────
        % A cada passo de tempo:
        %   1. Atualiza h_old com under-relaxation fisicamente motivada
        %      (h_old = +20 na zona saturada, -30 na nao saturada)
        %   2. Chama o pos-processador para salvar VTK
        %   3. Armazena h e theta nos pontos de monitoramento
        function [parms] = atualizarEstado(obj, env, parms, ...
                h, theta_n, time, count)

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
            postprocessor(h, theta_n, 0*flowresultZ, time, env, count, parms);
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
            end
            
            elem=env.geometry.elem;
            coord=env.geometry.coord;

            % quadrilatero
            figure(1)
            % tempo=2
            % MPFA-D e TPFA ortogonal
            A1=[300.0  65.0;
                100.0  70.0656;
                28.9426 78.4877];
            plot(A1(:,1), A1(:,2),'o')
            hold on
            % MPFA-D distorcido
            C1=[0	81.8198
                13.8135	81.1523
                27.8174	80.4857
                40.0743	78.4163
                65.6374	74.6297
                98.1984	72.2581
                144.064	68.8685
                187.473	68.2704
                300	65];
            plot(C1(:,1), C1(:,2))
            hold on
            % TPFA distorcido
            C2=[0	82.7826
                25.8440	81.3913
                51.6880	77.5652
                84.1676	73.7391
                148.079	68.8696
                207.101	67.4783
                300.000	65.0   ];
            plot(C2(:,1), C2(:,2))
            hold on
            B1=[0	80
                5.93714	79.9302
                11.5250	79.9302
                17.4622	79.5812
                22.7008	78.8831
                27.9395	78.5340
                35.2736	77.4869
                42.2584	76.4398
                49.5925	75.7417
                58.3236	74.3456
                65.6577	73.2984
                76.1350	71.9023
                100.233	69.4590
                112.806	68.0628
                129.220	67.0157
                139.697	65.9686
                150.524	65
                300.000	65 ];
            plot(B1(:,1), B1(:,2),'-')
            hold on
            % tempo=3
            % MPFA-D e TPFA
            A2=[300	65
                161.698	75.5158
                70.00	91.5765
                2.28041	    100];
            plot(A2(:,1), A2(:,2),'o')
            hold on
            % MPFA-D quadrilatero distorcido
            C3=[0	101.903
                14.5254	101.099
                34.8427	98.6191
                59.8216	93.6154
                84.7992	89.0345
                111.045	84.8801
                151.683	78.6520
                202.473	73.7211
                300.0	65.0];
            plot(C3(:,1), C3(:,2))
            hold on
            % TPFA quadrilatero distorcido
            C4=[0.0	100.524
                24.7094	99.1565
                53.0363	94.6519
                85.5607	88.7560
                129.974	81.8269
                205.859	73.1891
                300.0	65.0 ];
            plot(C4(:,1), C4(:,2))
            hold on

            B2=[0	100
                11.0337	99.9812
                23.5784	99.3283
                33.3361	97.9411
                46.2305	95.8582
                57.0343	93.7729
                68.8838	91.3391
                80.0365	88.9046
                91.5371	86.8201
                104.780	84.3880
                123.599	81.6127
                142.765	79.1875
                164.371	76.0
                196.779	73.3064
                228.141	70.5457
                266.473	67.0939
                300	65.0];
            plot(B2(:,1), B2(:,2),'-')
            hold on
            % tempo=4
            % MPFA-D e TPFA quad ortogonal
            A3=[300	65
                161.698	82.5585
                129.780	89.98
                70.00	101.622
                2.28041	110];
            plot(A3(:,1), A3(:,2),'o')
            hold on
            % MPFA-D quad distorcido

            C4=[0.0	109.185
                22.3412	108.202
                45.3857	104.159
                61.2612	101.125
                82.2579	97.0789
                112.474	90.4977
                159.584	83.9451
                200.038	77.8914
                245.610	72.8669
                300.0	65.0 ];
            plot(C4(:,1), C4(:,2),'o')
            hold on
            % TPFA quad distorcido
            C5=[0.0	108.668
                23.2068	107.400
                42.1941	104.440
                61.1814	99.7886
                90.2954	95.5603
                131.224	88.3721
                176.793	80.7611
                241.350	73.1501
                300.00	65.0];
            plot(C5(:,1), C5(:,2))
            hold on

            B3=[0	110
                9.28058	109.610
                19.3859	108.923
                28.0979	107.886
                39.2497	106.154
                53.1902	103.378
                67.1311	100.253
                82.4658	97.1292
                95.7094	94.3523
                112.438	90.8814
                126.378	88.8034
                146.242	85.6851
                183.879	80.4933
                233.364	73.9191
                300	65];
            plot(B3(:,1), B3(:,2),'-')
            hold on
            % tempo=8
            % MPFA-D e TPFA quad ortogonal
            A4=[300	65
                161.698	95.2342
                129.780	102.227
                70.0	114.354
                36.6371	119.324
                11.6036	119.705
                2.28041	120];
            plot(A4(:,1), A4(:,2),'o')
            hold on
            % MPFA-D  quad distorcido
            C6=[0.0	120.676
                18.4741	120.281
                36.2736	118.196
                53.6513	114.845
                79.0833	108.973
                119.348	101.435
                158.766	93.8956
                188.859	87.6088
                265.576	72.5268
                300.0	65.0];
            plot(C6(:,1), C6(:,2))
            hold on
            % TPFA  quad distorcido
            C7=[0.0	118.395
                20.5957	118.001
                43.4824	114.227
                74.8466	108.352
                114.688	100.797
                167.668	91.1465
                249.468	76.4625
                300.0	65.00];
            plot(C7(:,1), C7(:,2))
            hold on

            B4=[0.697674	120
                9.41860	119.930
                36.2791	118.531
                72.9070	110.839
                103.953	103.497
                134.302	97.2028
                162.907	91.6084
                215.930	81.8182
                257.093	73.4266
                300.000	65];

            plot(B4(:,1), B4(:,2),'-')

            xlabel('Aquifer Lenght')
            ylabel('Z')
            grid


            %% triangulo
            figure(2)
            % tempo=2
            % experimental
            A1=[300.0  65.0;
                100.0  70.0656;
                28.9426 78.4877];
            plot(A1(:,1), A1(:,2),'o')
            hold on
            % MPFA-D tri
            C1=[0.0	83.7209
                14.9059	83.2773
                35.1840	80.7117
                59.2637	77.2951
                94.3281	73.0173
                170.378	68.6818
                300.0	65.0];
            plot(C1(:,1), C1(:,2))
            hold on
            % TPFA tri
            C2=[0.0	70.1571
                10.4895	69.1099
                41.9580	69.4590
                103.147	68.0628
                154.895	66.3176
                242.308	65.0
                300.000	65.0
                ];
            plot(C2(:,1), C2(:,2))
            hold on
            % MPFAD tri distorcido
            B1=[0.0 	83.3005
                40.8202	79.1879
                90.3452	72.5621
                177.524	67.7338
                300.0	65.0 ];
            plot(B1(:,1), B1(:,2),'-')
            hold on
            % TPFA tri distorcido
            W1=[0.0	71.2058
                36.8553	70.9413
                74.9243	70.3318
                135.350	67.6801
                214.981	66.1201
                300.0	65.0];
            plot(W1(:,1), W1(:,2))
            hold on
            %--------------------------------------------------------------
            % tempo=3
            % experimental
            A2=[300	65
                161.698	75.5158
                70.00	91.5765
                2.28041	    100];
            plot(A2(:,1), A2(:,2),'o')
            hold on
            % MPFAD  tri
            C3=[0.0	104.001
                16.8049	101.918
                34.6039	100.259
                49.8616	97.7541
                71.4790	92.7318
                96.4861	87.7142
                133.358	82.2922
                168.110	77.7094
                213.032	73.1409
                300.00	65.00      ];
            plot(C3(:,1), C3(:,2))
            hold on
            % TPFA tri
            C4=[0.00000	92.6829
                11.1758	92.3345
                39.8137	92.6829
                61.8161	91.2892
                87.3108	88.1533
                122.584	82.9268
                163.097	77.0035
                232.247	70.3833
                300.00	65.00 ];
            plot(C4(:,1), C4(:,2))
            hold on
            % MPFA-D tri distocido
            B2=[0.00	103.832
                29.1552	101.328
                74.0825	92.5125
                115.180	85.0999
                176.834	76.9427
                239.885	70.5244
                300.00	65.00   ];
            plot(B2(:,1), B2(:,2),'-')
            hold on
            % TPFA tri distorcido
            W1=[0.241752	92.1743
                33.3812	92.2127
                55.0095	91.8900
                80.4774	89.4847
                110.134	84.6495
                152.699	78.7858
                205.032	72.9335
                300.273	65.0439];
            plot(W1(:,1), W1(:,2))
            hold on
            %--------------------------------------------------------------
            % tempo=4
            % experimental
            A3=[300	65
                161.698	82.5585
                129.780	89.98
                70.00	101.622
                2.28041	110];
            plot(A3(:,1), A3(:,2),'o')
            hold on
            % MPFA-D tri

            C4=[0.00	112.237
                9.93192	111.527
                21.1325	110.465
                34.0831	109.051
                51.5823	105.884
                68.0305	102.019
                96.0282	96.0421
                121.227	91.1176
                154.475	85.1348
                199.625	78.7884
                300.00	65.00];
            plot(C4(:,1), C4(:,2))
            hold on
            % TPFA tri
            C5=[0.00000	101.571
                8.03260	100.524
                33.5274	101.222
                52.7357	100.175
                85.9139	96.3351
                105.122	92.4956
                131.665	87.9581
                161.001	83.4206
                300.00	65.00  ];
            plot(C5(:,1), C5(:,2))
            hold on
            % MPFAD tri distorcido
            B3=[0.00000	112.391
                18.5315	111.344
                33.5664	109.599
                58.0420	103.316
                93.0070	96.6841
                128.322	90.4014
                176.573	83.0716
                227.622	74.6946
                300.00	65.00
                ];
            plot(B3(:,1), B3(:,2))
            hold on
            % TPFA tri distorcido
            W1=[0.00	101.394
                34.1860	100.000
                71.1628	97.2125
                104.651	92.3345
                151.395	84.6690
                202.326	77.3519
                258.488	70.3833
                300.00	65.00
                ];
            plot(W1(:,1), W1(:,2))
            hold on
            % tempo=8
            % experimental
            A4=[300	65
                161.698	95.2342
                129.780	102.227
                70.0	114.354
                36.6371	119.324
                11.6036	119.705
                2.28041	120];
            plot(A4(:,1), A4(:,2),'o')
            hold on
            % MPFA-D  tri
            C6=[0.00	124.297
                19.9319	122.762
                39.3526	119.693
                63.8842	114.066
                95.0596	107.417
                137.479	98.7212
                194.208	87.9795
                300.00	65.00];
            plot(C6(:,1), C6(:,2))
            hold on
            % TPFA  tri
            W1=[0.920582	110.491
                11.7605	109.804
                44.6276	109.842
                61.0616	109.512
                78.1977	106.734
                99.8809	102.564
                146.394	93.8769
                249.214	74.7659
                300.974	65.0361
                ];
            plot(W1(:,1), W1(:,2))
            hold on
            % MPFAD tri distorcido
            Z1=[0.00	122.996
                0.494247	122.996
                25.3463	121.225
                48.7957	117.365
                79.2422	110.012
                130.689	100.196
                166.386	92.1410
                220.982	81.6246
                300.00	65.00];
            plot(Z1(:,1), Z1(:,2))
            hold on
            % TPFA tri distorcido
            Z2=[0.00	108.902
                37.2862	108.247
                67.7074	106.886
                83.7937	104.810
                102.679	101.342
                129.259	95.7881
                157.587	90.5855
                204.100	82.2626
                300.00	65.00];

            plot(Z2(:,1), Z2(:,2))

            hold on
            xlabel('Aquifer Lenght')
            ylabel('Z')
            grid



            %% --------------------------------------------------------------------------
            if max(max(elem(:,4)))~=0
                figure(3)

                % Malha quadrilateral ortogonal
                %centro = obj.elementos_centroide_na_caixa(elem, coord, [-Inf 200], [20 25]);
                cacheFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'centro_x=20_quad.mat');
                if isfile(cacheFile), centro = load(cacheFile).centro;
                else
                 centro = obj.elementos_centroide_na_caixa(elem, coord, [-Inf 200], [20 25]); 
                 save(cacheFile, 'centro'); 
                end
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_quad_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_08_1WaterContent_steptime3.txt'));
                centelem = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_08_1centrocell3.txt'));
                centroY=centelem(centro,2);
                % MPFA-D
                % theta_n --> T=8
                theta_aux=theta_n(:,end);
                theta_init=theta_n(:,2);
                theta=theta_aux(centro);
                
                plot(theta, centroY)
                hold on
                % theta_n --> T=0
                plot(theta_init(centro),centroY)
                hold on

                % Malha quadrilateral distorcido
                % MPFA D
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_quad_distorcido_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_distorcido_08_1WaterContent_steptime3.txt'));
                % MPFA-D
                % theta_n --> T=8
                theta_aux=theta_n(:,end);

                theta=theta_aux(centro);

                plot(theta, centroY)
                hold on
                % TPFA
                % theta_n --> T=8
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_TPFA_quad_distorcido_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_TPFA_quad_distorcido_08_1WaterContent_steptime3.txt'));
                theta_aux=theta_n(:,end);

                theta=theta_aux(centro);

                plot(theta, centroY)
                hold on

                % theta experimental TEMPO=0
                T1=[0.309198	200-130.634
                    0.304946	200-121.653
                    0.319766	200-111.612
                    0.143656	200-100.410
                    0.142113	200-90.3880
                    0.0642259	200-81.4914
                    0.0558460	200-70.4409
                    0.0406789	200-61.1264
                    0.0473110	200-50.7499
                    0.0307680	200-40.7457
                    0.0156071	200-31.7768
                    0.00315489	200-21.7679
                    0.00980558	200-12.4283
                    ];
                plot(T1(:,1), T1(:,2),'o')
                hold on

                % theta experimental TEMPO= 8
               T2=[0.319091	200-100.900
                    0.301364	200-90.9002
                    0.283636	200-80.8976
                    0.271364	200-70.8871
                    0.271364	200-60.8625
                    0.253636	200-50.8693
                    0.249545	200-40.1613
                    0.250909	200-31.5192
                    0.250909	200-21.4883
                    0.253636	200-12.1518];

                plot(T2(:,1), T2(:,2),'o')
                xlabel('Water Content')
                ylabel('Z')
                title('x=21')
                grid
                %% ============================================================
                % x=80
                figure(4)
                
                cacheFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'centro_x=80_quad.mat');
                if isfile(cacheFile), centro_80 = load(cacheFile).centro_80;
                else
                 centro_80 = obj.elementos_centroide_na_caixa(elem, coord, [-Inf 200], [80 85]); 
                 save(cacheFile, 'centro_80'); 
                end

                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_quad_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_08_1WaterContent_steptime3.txt'));
                centelem = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_08_1centrocell3.txt'));
                centroY_80=centelem(centro_80,2);
                % MPFA-D
                % theta_n --> T=0 e 8
                theta_aux=theta_n(:,end);
                theta_init=theta_n(:,2);
                theta_80=theta_aux(centro_80);
                
                plot(theta_80, centroY_80)
                hold on
                plot(theta_init(centro_80),centroY_80)
                hold on
                % malha quadrilateral distorcido
                % MPFAD
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_quad_distorcido_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_distorcido_08_1WaterContent_steptime3.txt'));

                % theta_n --> T=0 e 8
                theta_aux=theta_n(:,end);
                theta_80=theta_aux(centro_80);
                
                plot(theta_80, centroY_80)
                hold on
                % TPFA
                % theta_n --> 8
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_TPFA_quad_distorcido_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_TPFA_quad_distorcido_08_1WaterContent_steptime3.txt'));
                theta_aux=theta_n(:,end);
                theta_80=theta_aux(centro_80);
                plot(theta_80, centroY_80)
                hold on

                % theta experimental x=80. t=0
                T3=[0.313699	200-121.724
                    0.324658	200-110.690
                    0.121918	200-100.345
                    0.153425	200-91.3793
                    0.0684932	200-81.3793
                    0.0602740	200-71.3793
                    0.0821918	200-61.3793
                    0.0917808	200-50.6897
                    0.0397260	200-41.0345
                    0.0383562	200-31.0345
                    0.0260274	200-21.7241
                    0.00273973	200-11.3793];

                plot(T3(:,1), T3(:,2),'o')
                hold on

                % theta experimental x=80. t=8
                T4=[0.319178	200-100.345
                    0.321918	200-91.0345
                    0.280822	200-81.0345
                    0.184932	200-71.7241
                    0.212329	200-60.6897
                    0.215068	200-51.3793
                    0.163014	200-41.3793
                    0.157534	200-31.3793
                    0.135616	200-21.0345
                    0.135616	200-11.3793
                    0.116438	200-3.10345];
                plot(T4(:,1), T4(:,2),'o')
                xlabel('Water Content')
                ylabel('Z')
                title('x=80')
                grid

                %% ============================================================
                figure(5)
                
                cacheFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'centro_x=140_quad.mat');
                if isfile(cacheFile), centro_140 = load(cacheFile).centro_140;
                else
                 centro_140 = obj.elementos_centroide_na_caixa(elem, coord, [-Inf 200], [140 145]); 
                 save(cacheFile, 'centro_140'); 
                end
                centroY_140=centelem(centro_140,2);

                % MPFAD
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_quad_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_08_1WaterContent_steptime3.txt'));

                theta_aux=theta_n(:,end);
                theta_140=theta_aux(centro_140);
                theta_init=theta_n(:,2);
                
                plot(theta_140, centroY_140)
                hold on
                plot(theta_init(centro_140),centroY_140)
                hold on
                % malha quadrilateral distorcido
                % MPFAD
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_quad_distorcido_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_distorcido_08_1WaterContent_steptime3'));
                theta_aux=theta_n(:,end);
                theta_140=theta_aux(centro_140);
                centroY_140=centelem(centro_140,2);
                plot(theta_140, centroY_140)
                hold on
                % TPFA
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_TPFA_quad_distorcido_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_TPFA_quad_distorcido_08_1WaterContent_steptime3.txt'));
                theta_aux=theta_n(:,end);
                theta_140=theta_aux(centro_140);
                centroY_140=centelem(centro_140,2);
                plot(theta_140, centroY_140)
                hold on

                % Theta experimental x=140. t=0
                T3=[0.271622	200-121.379
                    0.239189	200-111.379
                    0.187838	200-101.034
                    0.101351	200-91.3793
                    0.0729730	200-80.3448
                    0.0581081	200-71.3793
                    0.0567568	200-60.3448
                    0.0432432	200-51.0345
                    0.0391892	200-41.0345
                    0.0162162	200-21.7241
                    0.0148649	200-10.6897 ];

                plot(T3(:,1), T3(:,2),'o')
                hold on
                % theta experimental x=140. t=8
                T4=[0.278378	200-112.069
                    0.314865	200-101.379
                    0.275676	200-91.0345
                    0.259459	200-81.0345
                    0.228378	200-71.3793
                    0.117568	200-60.3448
                    0.0743243	200-51.0345
                    0.0527027	200-41.3793
                    0.0432432	200-31.7241
                    0.0135135	200-21.7241];
                plot(T4(:,1), T4(:,2),'o')

                title('x=140')
                grid

                %% =============================================================
                
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_quad_08';
                fname = fullfile(filepath);
                h_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_08_1h_steptime3.txt'));
                time2= readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_08_1time_step3.txt'));
                h_n(:, 1:2:end) = [];
                
               cacheFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'centros_6pontos_quadrilateral.mat');

                if isfile(cacheFile)
                    S = load(cacheFile);
                    centro1 = S.centro1; centro2 = S.centro2; centro3 = S.centro3;
                    centro4 = S.centro4; centro5 = S.centro5; centro6 = S.centro6;
                else
                    centro1 = obj.elemento_no_ponto(elem, coord, 12.5, 107.5);
                    centro2 = obj.elemento_no_ponto(elem, coord, 12.5, 132.5);
                    centro3 = obj.elemento_no_ponto(elem, coord, 12.5, 187.5);
                    centro4 = obj.elemento_no_ponto(elem, coord, 162.5, 82.5);
                    centro5 = obj.elemento_no_ponto(elem, coord, 162.5, 117.5);
                    centro6 = obj.elemento_no_ponto(elem, coord, 162.5, 157.5);
                    save(cacheFile, 'centro1', 'centro2', 'centro3', 'centro4', 'centro5', 'centro6');
                end

                h_time1 = h_n(centro1,:);  
                h_time2 = h_n(centro2,:);  
                h_time3 = h_n(centro3,:);  
                h_time4 = h_n(centro4,:);  
                h_time5 = h_n(centro5,:);  
                h_time6 = h_n(centro6,:);    
               %===========================================================
               % malha quadrilateral distorcido
               % MPFA-D
               filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_quad_distorcido_08';
                fname = fullfile(filepath);
                h_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_distorcido_08_1h_steptime3.txt'));
                h_n_MPFAD=h_n;
                h_n_MPFAD(:, 1:2:end) = [];
                 h_time1_M=h_n_MPFAD(centro1,:);
                 h_time2_M=h_n_MPFAD(centro2,:);
                 h_time3_M=h_n_MPFAD(centro3,:);
                 h_time4_M=h_n_MPFAD(centro4,:);
                 h_time5_M=h_n_MPFAD(centro5,:);
                 h_time6_M=h_n_MPFAD(centro6,:);
                % TPFA 
                 filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_TPFA_quad_distorcido_08';
                fname = fullfile(filepath);
                h_n = readmatrix(fullfile(fname, 'Tables_teste_TPFA_quad_distorcido_08_1h_steptime3.txt'));
                h_n_TPFA=h_n;
                h_n_TPFA(:, 1:2:end) = [];
                h_time1_T=h_n_TPFA(centro1,:);
                h_time2_T=h_n_TPFA(centro2,:);
                h_time3_T=h_n_TPFA(centro3,:);
                h_time4_T=h_n_TPFA(centro4,:);
                h_time5_T=h_n_TPFA(centro5,:);
                h_time6_T=h_n_TPFA(centro6,:);
                figure(6)
                plot(time2,h_time1)
                plot(time2,h_time1_M)
                plot(time2,h_time1_T)
                hold on
                plot(time2,h_time2)
                plot(time2,h_time2_M)
                plot(time2,h_time2_T)
                hold on
                plot(time2,h_time3)
                plot(time2,h_time3_M)
                plot(time2,h_time3_T)
                xlabel('Time')
                ylabel('Water content ')
                title('x=11')
                grid

                figure(7)
                plot(time2,h_time4)
                plot(time2,h_time4_M)
                plot(time2,h_time4_T)
                hold on
                plot(time2,h_time5)
                plot(time2,h_time5_M)
                plot(time2,h_time5_T)
                hold on
                plot(time2,h_time6)
                plot(time2,h_time6_M)
                plot(time2,h_time6_T)

                xlabel('Time')
                ylabel('Water Pressure ')
                title('x=161')
                grid
               
                %% ============================================================
            else
                figure(3)
                % Malha triangular
                % x=20
                cacheFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'centro_x=20_tri.mat');
                if isfile(cacheFile), centro = load(cacheFile).centro;
                else
                 centro = obj.elementos_centroide_na_caixa(elem, coord, [-Inf 200], [20 25]); 
                 save(cacheFile, 'centro'); 
                end

                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_tri_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_tri_08_1WaterContent_steptime3.txt'));
                centelem_aux = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_tri_08_1centrocell3.txt'));
                centroY=centelem_aux(centro,2);
                
                % MPFA-D
                theta_aux=theta_n(:,end);
                theta_init=theta_n(:,2);
                theta=theta_aux(centro);
                % theta_n --> T=0
                plot(theta_init(centro),centroY)
                hold on
                % theta_n --> T=8
                plot(theta, centroY)
                hold on
                
                % TPFA
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_TPFA_tri_08';
                fname = fullfile(filepath);
                theta_n_T = readmatrix(fullfile(fname, 'Tables_teste_TPFA_tri_08_1WaterContent_steptime3.txt'));
                % theta_n --> T=8
                theta_aux_T=theta_n_T(:,end);
                theta_T=theta_aux_T(centro);
                plot(theta_T, centroY)
                hold on
                % Malha triangular distorcido
                % MPFA D
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_tri_distorcido_08';
                fname = fullfile(filepath);
                theta_n_M = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_tri_distorcido_08_1WaterContent_steptime3.txt'));
                % theta_n --> T=8
                theta_aux_M=theta_n_M(:,end);
                theta_M=theta_aux_M(centro);

                plot(theta_M, centroY)
                hold on
                % TPFA
                % theta_n --> T=8
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_TPFA_tri_distorcido_08';
                fname = fullfile(filepath);
                theta_n_TP = readmatrix(fullfile(fname, 'Tables_teste_TPFA_tri_distorcido_08_1WaterContent_steptime3.txt'));
                theta_aux=theta_n_TP(:,end);
                theta_TP=theta_aux(centro);

                plot(theta_TP, centroY)
                legend('Initially solution','MPFA-D: triangle','TPFA: triangle','MPFA-D: distorted triangle','TPFA: distorted triangle')

                hold on

                % theta experimental TEMPO=0
                T1=[0.309198	200-130.634
                    0.304946	200-121.653
                    0.319766	200-111.612
                    0.143656	200-100.410
                    0.142113	200-90.3880
                    0.0642259	200-81.4914
                    0.0558460	200-70.4409
                    0.0406789	200-61.1264
                    0.0473110	200-50.7499
                    0.0307680	200-40.7457
                    0.0156071	200-31.7768
                    0.00315489	200-21.7679
                    0.00980558	200-12.4283
                    ];
                plot(T1(:,1), T1(:,2),'o')
                hold on

                % theta experimental TEMPO= 8
                T2=[0.319091	200-100.900
                    0.301364	200-90.9002
                    0.283636	200-80.8976
                    0.271364	200-70.8871
                    0.271364	200-60.8625
                    0.253636	200-50.8693
                    0.249545	200-40.1613
                    0.250909	200-31.5192
                    0.250909	200-21.4883
                    0.253636	200-12.1518];

                plot(T2(:,1), T2(:,2),'o')
                xlabel('Water Content')
                ylabel('Z')
                title('x=21')
                grid
                %% ============================================================
                figure(4)
                % x=80
              
                cacheFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'centro_x=80_tri.mat');
                if isfile(cacheFile), centro_80 = load(cacheFile).centro_80;
                else
                 centro_80 = obj.elementos_centroide_na_caixa(elem, coord, [-Inf 200], [80 85]); 
                 save(cacheFile, 'centro_80'); 
                end


                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_tri_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_tri_08_1WaterContent_steptime3.txt'));
                centelem = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_tri_08_1centrocell3.txt'));
                centroY_80=centelem(centro_80,2);
                % MPFA-D
                % theta_n --> T=0 e 8
                theta_aux=theta_n(:,end);
                theta_init=theta_n(:,2);
                theta_80=theta_aux(centro_80);
                plot(theta_80, centroY_80)
                
                hold on
                plot(theta_init(centro_80),centroY_80)
                
                hold on
                % TPFA

                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_TPFA_tri_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_TPFA_tri_08_1WaterContent_steptime3.txt'));
                % MPFA-D
                % theta_n --> T=0 e 8
                theta_aux=theta_n(:,end);
                theta_80=theta_aux(centro_80);
                plot(theta_80, centroY_80)
               
                hold on

                % malha triangular distorcido
                % MPFAD
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_tri_distorcido_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_tri_distorcido_08_1WaterContent_steptime3.txt'));

                % theta_n --> T=0 e 8
                theta_aux=theta_n(:,end);
                theta_80=theta_aux(centro_80);
                plot(theta_80, centroY_80)
                
                hold on
                % TPFA
                % theta_n --> 8
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_TPFA_tri_distorcido_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_TPFA_tri_distorcido_08_1WaterContent_steptime3.txt'));
                theta_aux=theta_n(:,end);
                theta_80=theta_aux(centro_80);
                plot(theta_80, centroY_80)
                legend('MPFA-D: triangle','Initially solution','TPFA: triangle','MPFA-D: distorted triangle','TPFA: distorted triangle')
                hold on

                % theta experimental x=80. t=0
                T3=[0.313699	200-121.724
                    0.324658	200-110.690
                    0.121918	200-100.345
                    0.153425	200-91.3793
                    0.0684932	200-81.3793
                    0.0602740	200-71.3793
                    0.0821918	200-61.3793
                    0.0917808	200-50.6897
                    0.0397260	200-41.0345
                    0.0383562	200-31.0345
                    0.0260274	200-21.7241
                    0.00273973	200-11.3793];

                plot(T3(:,1), T3(:,2),'o')
                hold on

                % theta experimental x=80. t=8
                T4=[0.319178	200-100.345
                    0.321918	200-91.0345
                    0.280822	200-81.0345
                    0.184932	200-71.7241
                    0.212329	200-60.6897
                    0.215068	200-51.3793
                    0.163014	200-41.3793
                    0.157534	200-31.3793
                    0.135616	200-21.0345
                    0.135616	200-11.3793
                    0.116438	200-3.10345];
                plot(T4(:,1), T4(:,2),'o')
                legend('MPFA-D: tri', 'TPFA:tri','MPFA-D: distorted tri', 'TPFA: distorted tri')
                xlabel('Water Content')
                ylabel('Z')
                title('x=80')
                grid
                %% ============================================================
                % X=140
                figure(5)

                cacheFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'centro_x=140_tri.mat');
                if isfile(cacheFile), centro_140 = load(cacheFile).centro_140;
                else
                 centro_140 = obj.elementos_centroide_na_caixa(elem, coord, [-Inf 200], [140 145]); 
                 save(cacheFile, 'centro_140'); 
                end

                % malha triangular
                % MPFAD
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_tri_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_tri_08_1WaterContent_steptime3.txt'));
                centelem = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_tri_08_1centrocell3.txt'));
                centroY_140=centelem(centro_140,2);

                theta_aux=theta_n(:,end);
                theta_140=theta_aux(centro_140);
                theta_init=theta_n(:,2);
                plot(theta_140, centroY_140)
                hold on
                plot(theta_init(centro_140),centroY_140)
                hold on
                % TPFA
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_TPFA_tri_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_TPFA_tri_08_1WaterContent_steptime3.txt'));
                theta_aux=theta_n(:,end);
                theta_140=theta_aux(centro_140);
                plot(theta_140, centroY_140)
                hold on
                % malha triangular distorcido
                % MPFAD
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_tri_distorcido_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_tri_distorcido_08_1WaterContent_steptime3'));
                theta_aux=theta_n(:,end);
                theta_140=theta_aux(centro_140);
                plot(theta_140, centroY_140)
                hold on
                % TPFA
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_TPFA_tri_distorcido_08';
                fname = fullfile(filepath);
                theta_n = readmatrix(fullfile(fname, 'Tables_teste_TPFA_tri_distorcido_08_1WaterContent_steptime3.txt'));
                theta_aux=theta_n(:,end);
                theta_140=theta_aux(centro_140);
                plot(theta_140, centroY_140)
                hold on

                % Theta experimental x=140. t=0
                T3=[0.271622	200-121.379
                    0.239189	200-111.379
                    0.187838	200-101.034
                    0.101351	200-91.3793
                    0.0729730	200-80.3448
                    0.0581081	200-71.3793
                    0.0567568	200-60.3448
                    0.0432432	200-51.0345
                    0.0391892	200-41.0345
                    0.0162162	200-21.7241
                    0.0148649	200-10.6897 ];

                plot(T3(:,1), T3(:,2),'o')
                hold on
                % theta experimental x=140. t=8
                T4=[0.278378	200-112.069
                    0.314865	200-101.379
                    0.275676	200-91.0345
                    0.259459	200-81.0345
                    0.228378	200-71.3793
                    0.117568	200-60.3448
                    0.0743243	200-51.0345
                    0.0527027	200-41.3793
                    0.0432432	200-31.7241
                    0.0135135	200-21.7241];
                plot(T4(:,1), T4(:,2),'o')
                xlabel(' Water content')
                ylabel('Z')
                title('x=140')
                legend('MPFA-D: tri', 'TPFA:tri','MPFA-D: distorted tri', 'TPFA: distorted tri')

                grid
                %%=============================================================
                %% =============================================================
                % malha triangular
                % MPFAD
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_tri_08';
                fname = fullfile(filepath);
                h_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_tri_08_1h_steptime3.txt'));
                time2= readmatrix(fullfile(fname, 'Tables_teste_MPFAD_tri_08_1time_step3.txt'));
                h_n_M=h_n;
                h_n_M(:, 1:2:end) = [];
                % TPFA
                filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_TPFA_tri_08';
                fname = fullfile(filepath);
                h_n = readmatrix(fullfile(fname, 'Tables_teste_TPFA_tri_08_1h_steptime3.txt'));
                h_n_T=h_n;
                h_n_T(:, 1:2:end) = [];

                cacheFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'centros_6pontos_tri.mat');

                if isfile(cacheFile)
                    S = load(cacheFile);
                    centro1 = S.centro1; centro2 = S.centro2; centro3 = S.centro3;
                    centro4 = S.centro4; centro5 = S.centro5; centro6 = S.centro6;
                else
                    centro1 = obj.elemento_no_ponto(elem, coord, 12.5, 107.5);
                    centro2 = obj.elemento_no_ponto(elem, coord, 12.5, 132.5);
                    centro3 = obj.elemento_no_ponto(elem, coord, 12.5, 187.5);
                    centro4 = obj.elemento_no_ponto(elem, coord, 162.5, 82.5);
                    centro5 = obj.elemento_no_ponto(elem, coord, 162.5, 117.5);
                    centro6 = obj.elemento_no_ponto(elem, coord, 162.5, 157.5);
                    save(cacheFile, 'centro1', 'centro2', 'centro3', 'centro4', 'centro5', 'centro6');
                end

                h_time1_M = h_n_M(centro1,:);  h_time1_T = h_n_T(centro1,:);
                h_time2_M = h_n_M(centro2,:);  h_time2_T = h_n_T(centro2,:);
                h_time3_M = h_n_M(centro3,:);  h_time3_T = h_n_T(centro3,:);
                h_time4_M = h_n_M(centro4,:);  h_time4_T = h_n_T(centro4,:);
                h_time5_M = h_n_M(centro5,:);  h_time5_T = h_n_T(centro5,:);
                h_time6_M = h_n_M(centro6,:);  h_time6_T = h_n_T(centro6,:);   
               %===========================================================
               % malha triangular distorcido
               % MPFA-D
               filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_MPFAD_tri_distorcido_08';
                fname = fullfile(filepath);
                h_n = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_tri_distorcido_08_1h_steptime3.txt'));
                h_n_MPFAD=h_n;
                h_n_MPFAD(:, 1:2:end) = [];
                 h_time1_MTD=h_n_MPFAD(centro1,:);
                 h_time2_MTD=h_n_MPFAD(centro2,:);
                 h_time3_MTD=h_n_MPFAD(centro3,:);
                 h_time4_MTD=h_n_MPFAD(centro4,:);
                 h_time5_MTD=h_n_MPFAD(centro5,:);
                 h_time6_MTD=h_n_MPFAD(centro6,:);
                % TPFA 
                 filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic409\teste_TPFA_tri_distorcido_08';
                fname = fullfile(filepath);
                h_n = readmatrix(fullfile(fname, 'Tables_teste_TPFA_tri_distorcido_08_1h_steptime3.txt'));
                h_n_TPFA=h_n;
                h_n_TPFA(:, 1:2:end) = [];
                h_time1_TTD=h_n_TPFA(centro1,:);
                h_time2_TTD=h_n_TPFA(centro2,:);
                h_time3_TTD=h_n_TPFA(centro3,:);
                h_time4_TTD=h_n_TPFA(centro4,:);
                h_time5_TTD=h_n_TPFA(centro5,:);
                h_time6_TTD=h_n_TPFA(centro6,:);
                figure(6)
                plot(time2,h_time1_M)
                 plot(time2,h_time1_T)
                plot(time2,h_time1_MTD)
                plot(time2,h_time1_TTD)
                hold on
                plot(time2,h_time2_M)
                plot(time2,h_time2_T)
                plot(time2,h_time2_MTD)
                plot(time2,h_time2_TTD)
                hold on
                plot(time2,h_time3_M)
                plot(time2,h_time3_T)
                plot(time2,h_time3_MTD)
                plot(time2,h_time3_TTD)

                legend('MPFA-D: tri', 'TPFA:tri','MPFA-D: distorted tri', 'TPFA: distorted tri')

                xlabel('Time')
                ylabel('Water content ')
                title('x=11')
                grid

                figure(7)
                plot(time2,h_time4_M)
                plot(time2,h_time4_T)
                plot(time2,h_time4_MTD)
                plot(time2,h_time4_TTD)
                hold on
                plot(time2,h_time5_M)
                plot(time2,h_time5_T)
                plot(time2,h_time5_MTD)
                plot(time2,h_time5_TTD)
                hold on
                plot(time2,h_time6_M)
                plot(time2,h_time6_T)
                plot(time2,h_time6_MTD)
                plot(time2,h_time6_TTD)

                xlabel('Time')
                ylabel('Water Pressure ')
                title('x=161')
                grid
            end

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
