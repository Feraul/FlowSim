classdef Caso436 < SimulacaoBase
    %======================================================================
    % CASO436 — Teste de convergencia da equacao de Richards com
    % condutividade hidraulica nao-linear (dependente de h).
    %
    % Ideia geral do caso:
    %   Este e um teste do tipo "solucao manufaturada" (Method of
    %   Manufactured Solutions, MMS). Em vez de comparar com uma solucao
    %   analitica "natural" do problema fisico, escolhe-se A PRIORI uma
    %   funcao suave qualquer,
    %
    %       u_exact(x,y,t) = -3*t*x*(1-x)*y*(1-y) - 1
    %
    %   definida no dominio unitario [0,1] x [0,1], e usa-se essa funcao
    %   para:
    %     (a) impor a condicao de contorno de Dirichlet (configurarContorno),
    %     (b) calcular o erro numerico em cada passo de tempo, comparando
    %         a pressao h obtida pelo solver com u_exact avaliada nos
    %         mesmos pontos (atualizarEstado).
    %
    %   Como x*(1-x) >= 0 e y*(1-y) >= 0 dentro do dominio, o termo
    %   -3*t*x*(1-x)*y*(1-y) e sempre <= 0, e subtraindo mais 1 garante
    %   que u_exact seja sempre negativo (h < 0) em todo o dominio e em
    %   todo instante t > 0. Fisicamente isso significa que o dominio
    %   inteiro permanece na zona NAO SATURADA, o que faz sentido para
    %   testar justamente o termo nao-linear kr(h) da condutividade
    %   hidraulica relativa (o objetivo do caso e verificar a ordem de
    %   convergencia do esquema numerico quando kr(h) e nao-linear).
    %
    %   O modelo de retencao de agua no solo usado aqui, apesar de os
    %   comentarios originais dizerem "Brooks-Corey", e na verdade o
    %   modelo de van Genuchten-Mualem:
    %
    %       Se(h) = theta(h) = (1 + (alpha*|h|)^n)^(-(n-1)/n)      se h<0
    %       kr(h) = sqrt(theta) * (1 - (1 - theta^(n/(n-1)))^((n-1)/n))^2
    %
    %   com theta_s = 1 e theta_r = 0 (forma simplificada/normalizada).
    %   Essa e a formula classica de van Genuchten (1980) combinada com
    %   o modelo de condutividade relativa de Mualem (1976).
    %
    %   A malha e resolvida com um esquema de volumes finitos (MPFA/LPEW2,
    %   a julgar pelos nomes bedge/inedge/elemarea/normals e pela funcao
    %   calcularTermoNeumannVet), e a nao-linearidade de Richards e tratada
    %   por iteracao de Picard (h_old sendo o "chute" da iteracao anterior).
    %======================================================================
    properties
        Nome   = 'Convergence test: Non-Linear Hydraulic Conductivity'   % nome legivel usado nos logs/relatorios da simulacao
        TipoID = 436                     % identificador numerico do caso, deve bater com "numcase" no arquivo Start.dat que seleciona este caso
    end

    methods
        % ── 1. Permeabilidade relativa — (rotulado Brooks-Corey, mas a formula e van Genuchten-Mualem) ─────────
        % Calcula, para cada elemento da malha, a condutividade hidraulica
        % nao-linear K(h) = Ksat * kr(h) a partir da carga de pressao
        % "old" (h_old, valor da iteracao de Picard anterior), e monta o
        % tensor de permeabilidade isotropico usado na montagem do
        % sistema linear (MPFA).
        function [env, parms] = configurarPermeabilidade(obj, env, parms, time)
            elem   = env.geometry.elem;      % tabela de conectividade dos elementos da malha
            nelem  = size(elem,1);           % numero total de elementos (celulas) da malha
            theta  = ones(nelem,1);          % conteudo volumetrico de agua; inicializa em 1 (= saturado, theta_s) para todos os elementos
            kr     = ones(nelem,1);          % condutividade relativa; inicializa em 1 (K = Ksat) para todos os elementos
            n      = parms.nvg;              % parametro "n" do modelo de van Genuchten
            alpha  = parms.alpha;            % parametro alpha (inverso do "entry pressure") do modelo de van Genuchten
            h      = parms.h_old;            % carga de pressao da iteracao anterior de Picard (chute atual da nao-linearidade)
            neg    = h <= 0;                 % mascara logica: elementos na zona nao saturada (h <= 0)
            % Nos elementos nao saturados, substitui theta=1 e kr=1 pelas
            % formulas nao-lineares de van Genuchten-Mualem:
            theta(neg) = (1 + (-alpha*h(neg)).^n).^(-(n-1)/n);
            kr(neg)    = sqrt(theta(neg)) .* (1 - (1 - theta(neg).^(n/(n-1))).^((n-1)/n)).^2;
            % Condutividade hidraulica final = Ksat (env.config.perm(1,1)) vezes kr(h):
            coef                    = env.config.perm(1,1) .* kr;
            env.config.kmap         = obj.iso(env, coef);      % monta o tensor isotropico K = coef * Identidade, formato exigido pela montagem MPFA
            parms.auxperm           = env.config.kmap;         % guarda uma copia do tensor atual em parms, para uso posterior (ex.: ajuste de contorno)
            env.config.auxkmap      = obj.isoConst(env);        % tensor de referencia com a permeabilidade SATURADA constante (Ksat), sem depender de h
            env.geometry.elem(:,5)  = env.utils.idx;             % atualiza a coluna de indices dos elementos (associa cada elemento ao seu indice de material/regiao)
        end

        % Indica ao solver que a permeabilidade deve ser recalculada a
        % cada iteracao de Picard / passo de tempo, pois ela depende de h
        % (problema nao-linear).
        function flag = precisaAtualizarPermeabilidade(obj)
            flag = true;
        end

        % ── 2. Condicoes de contorno de Dirichlet ─────────────────
        % Atribui o valor de h nas faces/vertices de contorno.
        % Observacao: apesar do comentario original no cabecalho da
        % classe descrever um esquema com 3 flags distintos (entrada,
        % nivel piezometrico linear h=65-z, dreno h=0), a IMPLEMENTACAO
        % efetiva deste metodo usa, para TODO ponto de contorno, o valor
        % da solucao manufaturada u_exact(x,y,t) avaliada no ponto medio
        % da aresta/no de contorno. Ou seja, neste caso (436) a condicao
        % de contorno de Dirichlet e simplesmente "h = solucao exata",
        % tipico de testes de convergencia por solucao manufaturada.
        function bcattrib = configurarContorno(obj, vertices, flagptr, time, env, pR)
            % calcula o ponto medio de cada aresta (coordmid)
            if size(vertices,2) > 1
                % vertices tem 2 colunas → representa uma aresta (v1, v2);
                % calcula o ponto medio da aresta a partir das coordenadas
                % dos dois vertices que a formam.
                coordmid         = (env.geometry.coord(vertices(:,1),1:2) + ...
                    env.geometry.coord(vertices(:,2),1:2)) ./ 2;
                auxvertices(:,1) = 1:size(flagptr,1);  % reindexa sequencialmente (1..N) para casar com o vetor de saida bcattrib
                vertices         = auxvertices;
            else
                % vertices tem 1 coluna → lista de indices de nos (vertices)
                % de contorno; usa diretamente as coordenadas dos nos.
                coordmid = env.geometry.coord(:, 1:2);
            end
            % Avalia a solucao manufaturada exata no ponto medio (x,y) e no
            % instante "time" atual, e usa esse valor como carga de pressao
            % prescrita no contorno:
            %   u_exact(x,y,t) = -3*t*x*(1-x)*y*(1-y) - 1
            bcattrib= -3.*time.*coordmid (:,1).*(1-coordmid (:,1)).*coordmid (:,2).*(1-coordmid (:,2))-1;
        end

        % ── 3. Flags de contorno (nflag e nflagface) ──────────────
        % Monta os vetores de flags para vertices e faces de contorno.
        % nflag(:,1)    → tipo de BC (Dirichlet < 200, Neumann >= 200)
        % nflag(:,2)    → valor prescrito de h no vertice
        % nflagface(:,1) → tipo de BC na face
        % nflagface(:,2) → valor prescrito de h na face
        function [nflag, nflagface] = configurarFlags(obj, env, pR, time)
            nelem_nodes = size(env.geometry.coord, 1);   % numero total de nos (vertices) da malha
            nelem_faces = size(env.geometry.bedge, 1);    % numero total de faces de contorno (bedge = boundary edges)

            % inicializa com valor sentinela (5000 = interior, ou seja,
            % "nenhuma condicao de contorno atribuida ainda")
            nflag     = 5000 * ones(nelem_nodes, 2);
            nflagface = zeros(nelem_faces, 2);

            % mapeia bcflag → vertices e faces de contorno: obtem, a
            % partir da geometria/flags de entrada (arquivo de malha),
            % quais nos e faces pertencem ao contorno e qual o codigo de
            % condicao de contorno (bcflag) associado a cada um.
            [vertex_idx, face_idx, bcflag_vertex, bcflag_face, ...
                bc_row_vertex, bc_row_face] = obj.prepararIndices(env);

            % preenche flags e valores nos vertices de contorno
            nflag(vertex_idx,1) = bcflag_vertex(:,1);          % tipo de condicao de contorno para cada no de contorno
            mmmm                = obj.configurarContorno(vertex_idx, bc_row_vertex', time, env, pR); % avalia u_exact em cada no
            nflag(vertex_idx,2) = mmmm(vertex_idx);            % armazena o valor prescrito (Dirichlet) para cada no de contorno

            % preenche flags e valores nas faces de contorno
            nflagface(:,1) = bcflag_face(:,1);                 % tipo de condicao de contorno para cada face de contorno
            nflagface(:,2) = obj.configurarContorno(...
                [face_idx(:,1) face_idx(:,2)], bc_row_face', time, env, pR); % avalia u_exact no ponto medio de cada face de contorno
        end

        % ── 4. Pre-processamento fisico do caso ───────────────────
        % Define os parametros fisicos e as condicoes iniciais antes do
        % inicio do loop temporal.
        function [parms, env] = preprocessar(obj, env, parms)
            parms.alpha=0.1844;     % parametro alpha do modelo de van Genuchten (cm^-1 ou unidade equivalente do dominio)
            parms.nvg=3;            % parametro n do modelo de van Genuchten (controla a forma/curvatura da curva de retencao)
            % Carga de pressao inicial "verdadeira" (usada apenas como
            % valor de referencia / placeholder aqui — o valor real de
            % h_init passa a ser a solucao exata avaliada em t, atualizada
            % dentro de atualizarEstado a cada passo):
            parms.h_init=-1*ones(size(env.geometry.elem,1),1);
            % Chute inicial (h_old) para a primeira iteracao de Picard, um
            % valor uniforme e ligeiramente negativo (zona nao saturada),
            % coerente com o fato de u_exact ser sempre negativo neste caso:
            parms.h_old=-5*ones(size(env.geometry.elem,1),1);
            % passo de tempo fixo, dt = 1/128 (dias, ou unidade de tempo do
            % problema) — valor pequeno tipico de testes de convergencia,
            % onde se quer isolar o erro espacial/temporal do esquema.
            parms.dt= (1/128);
        end

        % ── 5. Fontes e pocos ─────────────────────────────────────
        % Delega a definicao de pocos injetores/produtores (termos fonte)
        % para a funcao padrao defineWells. Este caso nao possui pocos
        % proprios; usa a implementacao generica da base.
        function source_wells = definirFontes(obj, env, parms,time)
            source_wells = defineWells(env, parms,time);
        end

        % ── 6a. Modelo de retencao hidrica — van Genuchten ─────────
        % Calcula o conteudo volumetrico de agua theta(h):
        %   theta = theta_s = 1                             se h >= 0  (saturado)
        %   theta = (1 + (alpha*|h|)^n)^(-(n-1)/n)           se h <  0  (nao saturado, van Genuchten)
        % (theta_r = 0 e theta_s = 1 nesta forma normalizada/simplificada)
        function theta = calcularTheta(obj, h, parms)
            alpha=parms.alpha;
            nvg=parms.nvg;
            theta = zeros(size(h));  % inicializa o vetor (obs.: caso h>=0 nao seja tratado, o valor fica 0; ver observacao abaixo)
            idx_neg = h < 0;          % indices onde h < 0 (zona nao saturada)
            % h < 0: aplica a curva de retencao de van Genuchten
            theta(idx_neg) = (1 + (-alpha*h(idx_neg)).^nvg).^(-(nvg-1)/nvg);
            % OBS: diferente de configurarPermeabilidade/exataAuxiliares
            % (que inicializam theta=1 para h>=0), aqui theta permanece 0
            % para h>=0 pois nao ha um ramo "else" explicito — atencao a
            % essa inconsistencia caso h>=0 ocorra em algum elemento.
        end

        % ── 6b. Capacidade hidrica especifica — dtheta/dh ─────────
        % Derivada analitica de theta em relacao a h (van Genuchten):
        %   C(h) = dtheta/dh
        % Usada na montagem da matriz de massa/capacidade em
        % soil_properties, para o metodo de Picard aplicado a forma
        % nao-linear da equacao de Richards (o termo de acumulacao
        % d(theta)/dt e linearizado via C(h)*dh/dt).
        function dthetadh = calcularCapacidade(obj, h, parms)
            alpha=parms.alpha;
            nvg=parms.nvg;
            % calculo da capacidade hidrica especifica (derivada de theta)
            dthetadh = zeros(size(h));      % inicializa em zero (na zona saturada h>0 a capacidade e nula, pois theta=theta_s=constante)
            idx = (h <= 0);                  % zona nao saturada
            dthetadh(idx)  = alpha*(nvg-1) .* ((-alpha*h(idx)).^(nvg-1)) .* ...
                (((-alpha*h(idx)).^nvg + 1).^((1/nvg) - 2));
        end

        % ── 7. Flowrate boundary ──────────────────────────────────
        % Indica que este caso NAO possui contribuicao gravitacional
        % adicional nas faces de contorno (diferente do caso 439, citado
        % no comentario original, que teria essa contribuicao).
        function flag = temFlowrateBoundary(obj)
            flag = false;
        end

        % Como temFlowrateBoundary = false, este metodo nao precisa alterar
        % nada: mantem o flowrateZ (vazao na direcao z / gravitacional)
        % inalterado, equivalente ao comportamento padrao (fallback) da
        % classe base.
        function flowrateZ = ajustarFlowrate(obj, flowrateZ, bedge)
            % nao faz nada — fallback da base seria identico
        end

        % ── 8. Permeabilidade na fronteira — (metodo nao utilizado neste caso) ─────────
        % Assinatura prevista para ajustar K11/K22 nas faces de contorno
        % com Dirichlet usando a permeabilidade relativa avaliada em
        % h_contorno (ex.: coef = Ksat * kr(h_contorno)). Neste caso o
        % corpo esta vazio, ou seja, a permeabilidade de contorno usada e
        % simplesmente a calculada globalmente em configurarPermeabilidade
        % (nenhum ajuste especial e feito face a face).
        function [K11, K12, K21, K22] = ajustarKContorno(obj, env, parms, ...
                auxkmap, matid, h_contorno, maskT)
        end

        % ── 9. Termo temporal de Richards ─────────────────────────
        % Adiciona a matriz de capacidade hidrica (proporcional a C(h) =
        % dtheta/dh) e o vetor de acumulacao ao sistema linear global —
        % corresponde a discretizacao do termo d(theta)/dt na equacao de
        % Richards. A montagem propriamente dita (matriz esparsa em bloco
        % diagonal, dividida por dt) e delegada a funcao externa
        % soil_properties, que tem acesso a malha completa (env) e aos
        % parametros fisicos (parms).
        function [M,I] = adicionarTermoTemporal(obj, M, I, parms, flowresultZ, env)
            [M,I] = soil_properties(M, I, parms, flowresultZ, env);
        end

        % ── 10. Interpolacao de Neumann — LPEW2 ───────────────────
        % Calcularia o termo "s" usado nos nos de contorno com condicao de
        % Neumann, para a interpolacao de pesos LPEW2 (Local Point-based
        % Extended Weighted Interpolation, usada no esquema MPFA-D para
        % reconstruir valores nodais a partir dos valores nos centroides).
        % Como este caso so possui contorno de Dirichlet (ver
        % configurarContorno), a rotina retorna simplesmente "false",
        % desativando esse calculo.
        function s = calcularTermoNeumannVet(obj, r, sum_lambda, N, env)
            s=false;
        end
    end

    methods
        % ── 11. Inicializacao antes do loop temporal ───────────────
        % Prepara as estruturas de acumulacao de erro (hnsum, vnsum) e
        % inicializa as series de armazenamento (solucao exata, theta e
        % kmap "analiticos") no instante inicial "time", avaliando a
        % solucao manufaturada e as propriedades hidraulicas dela
        % decorrentes em cada centroide da malha.
        function [parms, extras] = inicializar(obj, env, parms, time)
            extras.hnsum = 0;   % acumulador do erro de pressao (norma L2 integrada no tempo) — comeca zerado
            extras.vnsum = 0;   % acumulador do erro de velocidade/fluxo (norma L2 integrada no tempo) — comeca zerado
            extras.MBE      = 0;  % erro de balanco de massa por passo (indice 1 = t=0, nao usado)
            extras.MBE_rel  = 0;  % MBE normalizado por passo
            extras.maxRloc  = 0;  % maior residuo local por elemento, por passo
            extras.L1_MBE   = 0;  % norma L1 intensiva do residuo, por passo
            extras.L2_MBE   = 0;  % norma L2 intensiva do residuo, por passo
            centelem = env.geometry.centelem;   % coordenadas dos centroides de todos os elementos da malha
            x = centelem(:,1); y = centelem(:,2); t = time;
            % Avalia a solucao exata manufaturada no instante inicial em
            % cada centroide:
            u_exact0 = -3*t .* x.*(1-x) .* y.*(1-y) - 1;
            % Calcula theta e kmap (condutividade) "analiticos" a partir
            % de u_exact0, usando o mesmo modelo de van Genuchten-Mualem
            % (funcao auxiliar exataAuxiliares, item 12b abaixo):
            [theta0, kmap0] = obj.exataAuxiliares(u_exact0, parms, env);
            % Guarda, para cada centroide, o par (coordenada y, valor) das
            % grandezas exatas — usado depois para exportar/plotar perfis
            % ao longo de uma coluna (ex.: x fixo, variando y):
            extras.exact_solution_storage  = [centelem(:,2), u_exact0];
            extras.theta_storage_analitica = [centelem(:,2), theta0];
            extras.kmap_storage_analitica  = [centelem(:,2), kmap0];
        end

        % ── 12. Atualizacao dentro do loop temporal ────────────────
        % Executada a cada passo de tempo, depois que o solver obtem a
        % nova pressao numerica h. Responsavel por:
        %   1. Atualizar parms.p_old (indicador de zona saturada/nao
        %      saturada) e parms.h_init (usado como referencia da solucao
        %      exata no passo atual).
        %   2. Chamar o pos-processador para salvar os resultados em VTK.
        %   3. Calcular a solucao exata e as propriedades hidraulicas
        %      exatas (theta, K) nos centroides E nas faces de contorno e
        %      internas (bedge/inedge), e a partir delas o fluxo/velocidade
        %      exata em cada face.
        %   4. Comparar essas grandezas exatas com o resultado numerico do
        %      solver (h e flowrate) para acumular o erro L2 de pressao
        %      (hnsum) e de velocidade (vnsum), integrados no tempo (soma
        %      ponderada por dt a cada passo — um metodo tipo
        %      "quadratura retangular no tempo").
        function [parms,extras] = atualizarEstado(obj, env, parms, ...
                h, theta_n, time, count,flowrate,extras)
            centelem=env.geometry.centelem;
            x = centelem(:,1);
            y = centelem(:,2);
            t = time;
            nvG=parms.nvg;

            % ── MBE: captura theta_old ANTES de parms.h_init ser sobrescrito ──
            % parms.h_init aqui ainda e o h_n que soil_properties usou como
            % nivel de tempo anterior PARA ESTE passo (so e sobrescrito mais
            % abaixo, com o valor exato do novo tempo).
            theta_old_MBE = obj.calcularTheta(parms.h_init, parms);

            % Solucao exata manufaturada no instante atual, avaliada nos
            % centroides dos elementos:
            u_exact = -3*t .* x.*(1-x) .* y.*(1-y) - 1;
            theta_exact = thetafunction(u_exact, parms, env);
            k_theta_exact= obj.calcularKappa(theta_exact, 0.03, 1, nvG);
            k_theta= obj.calcularKappa(theta_n, 0.03, 1, nvG);

            extras.u_exact=u_exact;
            extras.theta_exact=theta_exact;
            extras.k_theta_exact=k_theta_exact;
            extras.theta_n=theta_n;
            extras.k_theta=k_theta;
            extras.h=h;
            %% calcula theta e salva VTK
            % Chama o pos-processador generico, passando a pressao
            % numerica h, o theta numerico (theta_n) e a solucao exata
            % (exact_sol), para escrita de arquivos de visualizacao VTK.
            postprocessor(env, count, time, pressure=h, theta_n=theta_n,...
                exact_sol=u_exact, k_theta=k_theta,k_theta_exact=k_theta_exact,...
                theta_exact=theta_exact);

            % Classifica os elementos pela pressao NUMERICA h obtida pelo
            % solver (nao pela exata): p_oldaux1 = zona nao saturada
            % (h<=0), p_oldaux2 = zona saturada (h>0). parms.p_old recebe
            % um "rotulo" numerico (10 para saturado, -5 para nao
            % saturado) — provavelmente usado apenas para diagnostico ou
            % como novo chute em outra parte do solver, ja que h_old em si
            % e atualizado em configurarPermeabilidade a partir de h em si.
            p_oldaux1 = (h <= 0);
            p_oldaux2 = (h > 0);
            parms.p_old     = 10*p_oldaux2 - 5*p_oldaux1;
            parms.h_init=u_exact;   % atualiza a "referencia exata" do passo atual, usada em outros pontos do fluxo (ex.: relatorios)
            [env, parms] = PLUG_kfunction(env, parms, time);
            
            %==============================================================
            % A partir daqui: calculo do fluxo/velocidade EXATA (dada pela
            % lei de Darcy aplicada a solucao manufaturada) em cada face da
            % malha, para depois comparar com o fluxo NUMERICO (flowrate)
            % retornado pelo solver, e assim estimar o erro de velocidade.
            centelem = env.geometry.centelem;
            bedge    = env.geometry.bedge;    % faces de contorno (boundary edges)
            inedge   = env.geometry.inedge;   % faces internas (interior edges)
            coord    = env.geometry.coord;    % coordenadas dos nos da malha
            normals  = env.geometry.normals;   % vetores normais de cada face (confirme o nome exato desse campo no seu env)
            elemarea = env.geometry.elemarea;  % area de cada elemento da malha
            alpha = parms.alpha;
            nvg   = parms.nvg;
            dt    = parms.dt;
            t     = time;
            sizebedge       = size(bedge,1);              % numero de faces de contorno
            sizebedgeinedge = sizebedge + size(inedge,1);  % numero total de faces (contorno + internas)

            %------------------------------------------------------------------
            % Bloco das faces de CONTORNO (bedge)
            v1 = bedge(:,1); v2 = bedge(:,2);   % indices dos dois nos que formam cada face de contorno
            m_bedge = 0.5*(coord(v1,1:2)+coord(v2,1:2));   % ponto medio de cada face de contorno
            areanormal = zeros(sizebedgeinedge,1);
            % comprimento de cada face de contorno (usado para converter
            % vazao total na face em velocidade normal media):
            areanormal(1:sizebedge,1) = norm(coord(v1,1:2)-coord(v2,1:2));
            xbedge = m_bedge(:,1); ybedge = m_bedge(:,2);
            % Solucao exata avaliada no ponto medio de cada face de contorno:
            u_exact_bedge = -3*t .* xbedge.*(1-xbedge) .* ybedge.*(1-ybedge) - 1;
            % Gradiente ANALITICO de u_exact (derivadas parciais exatas em
            % relacao a x e y), obtido derivando a expressao polinomial da
            % solucao manufaturada:
            dpdx = -3*t .* (1 - 2*xbedge) .* ybedge.*(1-ybedge);
            dpdy = -3*t .* xbedge.*(1-xbedge) .* (1 - 2*ybedge);
            % Condutividade relativa exata kr(u_exact) nas faces de contorno
            % (mesmo modelo de van Genuchten-Mualem, ramo simplificado
            % kr=theta, usado aqui como aproximacao para o calculo do
            % fluxo de Darcy exato):
            ktheta = ones(size(u_exact_bedge));
            mask = u_exact_bedge <= 0;
            ktheta(mask) = (1 + (-alpha*u_exact_bedge(mask)).^nvg).^(-(nvg-1)/nvg);
            % Velocidade de Darcy exata: v = -K(h) * grad(h)  (sem termo
            % gravitacional, pois temFlowrateBoundary = false neste caso):
            vx1 = -ktheta .* dpdx;  vy1 = -ktheta .* dpdy;
            vn = zeros(sizebedgeinedge,1);
            % Projeta a velocidade exata na direcao normal da face, para
            % obter o fluxo normal exato "vn" (comparavel ao fluxo
            % numerico do solver):
            vn(1:sizebedge,1) = vx1 .* normals(1:sizebedge,1) + vy1 .* normals(1:sizebedge,2);

            %------------------------------------------------------------------
            % Bloco das faces INTERNAS (inedge) — mesma logica do bloco
            % anterior, agora aplicada as faces internas da malha.
            v11 = inedge(:,1); v21 = inedge(:,2);
            m_inedge = 0.5*(coord(v11,1:2)+coord(v21,1:2));   % ponto medio de cada face interna
            areanormal(sizebedge+1:sizebedgeinedge,1) = norm(coord(v11,1:2)-coord(v21,1:2));   % comprimento de cada face interna
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
            % Fluxo/velocidade NUMERICO: o solver retorna "flowrate" como
            % vazao total atraves de cada face; dividindo pelo comprimento
            % (areanormal) obtem-se a velocidade normal media numerica,
            % diretamente comparavel a "vn" (exata) calculada acima:
            vn_numerico = flowrate ./ areanormal;
            %----------------------------------------------------------------
            % Erro de PRESSAO: soma ponderada pela area de cada elemento
            % (media espacial do erro quadratico), acumulada no tempo
            % multiplicando por dt a cada passo (integracao temporal do
            % erro L2 espacial, ao estilo de uma quadratura retangular):
            wA = elemarea;
            diffp = u_exact - h;
            extras.hnsum = extras.hnsum + dt*(sum(wA .* diffp.^2) / sum(wA));

            % Peso "Q" de cada face para o calculo do erro de velocidade:
            % para faces de contorno, usa a area do elemento vizinho unico
            % (bedge(:,3)); para faces internas, usa a media das areas dos
            % dois elementos vizinhos (inedge(:,3) e inedge(:,4)):
            Q = zeros(size(inedge,1) + size(bedge,1), 1);
            Q(1:size(bedge,1)) = elemarea(bedge(:,3));
            Q(size(bedge,1)+1:end) = 0.5.*(elemarea(inedge(:,3)) + elemarea(inedge(:,4)));

            % Erro de VELOCIDADE: diferenca entre velocidade normal exata
            % (vn) e numerica (vn_numerico), ponderada por Q e acumulada
            % no tempo da mesma forma que o erro de pressao:
            e  = vn - vn_numerico;
            extras.vnsum = extras.vnsum + dt*(Q'*e.^2)/sum(Q');

            %% ── MBE (balanco de massa) deste passo ────────────────────
            % sourcevector deve ser avaliado no MESMO tempo usado no solve
            % que gerou "h" (h_new): como source_wells e recalculado em
            % hydraulic_RE APOS o incremento de "time", a fonte realmente
            % usada neste passo foi avaliada em (time - dt), nao em "time".
            time_old_MBE = time - dt;
            [flowrate_MBE, flowresult_MBE, ~, ~] = env.metodo.calcularFlowrate(h, env, parms);
            sourcevector_MBE = PLUG_sourcefunction(h, env, time_old_MBE, parms);

            [MBEval, MBE_rel, Rloc, L1_MBE, L2_MBE] = ferncodes_MBE(theta_old_MBE, theta_n, ...
                flowrate_MBE, flowresult_MBE, sourcevector_MBE, elemarea, dt, ...
                sizebedge, +1);

            extras.MBE(count)          = MBEval;
            extras.MBE_rel(count)      = MBE_rel;
            extras.maxRloc(count)      = max(abs(Rloc));
            extras.L1_MBE(count)       = L1_MBE;   % norma L1 intensiva (para grafico log-log)
            extras.L2_MBE(count)       = L2_MBE;   % norma L2 intensiva (para grafico log-log)
            extras.time_storage(count) = time;   % para plotar MBE(t) em finalizar

            % Guarda os campos do passo atual para uso posterior (ex.: em
            % finalizar, para o calculo do erro final consolidado):
            extras.u_exact=u_exact;
            extras.h=h;

            % Recalcula theta e kmap "exatos" (a partir de u_exact) neste
            % passo de tempo, e acumula nas series historicas de
            % armazenamento (uma coluna de coordenada y + uma coluna de
            % valor, por passo de tempo "count"):
            [theta_a, kmap_a] = obj.exataAuxiliares(u_exact, parms, env);
            extras.exact_solution_storage(:, 2*count-1:2*count)  = [centelem(:,2), u_exact];
            extras.theta_storage_analitica(:, 2*count-1:2*count) = [centelem(:,2), theta_a];
            extras.kmap_storage_analitica(:, 2*count-1:2*count)  = [centelem(:,2), kmap_a];
        end

        % ── 12b. Funcao auxiliar: propriedades exatas a partir de u_exact ─────
        % Calcula theta(u_exact) e kmap(u_exact) usando o modelo de van
        % Genuchten-Mualem — EXATAMENTE a mesma formula usada em
        % configurarPermeabilidade (item 1), mas aqui aplicada a solucao
        % exata manufaturada em vez de a pressao numerica h_old. Serve
        % como "gabarito" para comparar com os valores numericos theta_n e
        % com a permeabilidade calculada pelo solver.
        function [theta, kmap] = exataAuxiliares(obj, u_exact, parms, env)
            alpha = parms.alpha; nvg = parms.nvg;
            theta = ones(size(u_exact));   % zona saturada (default): theta = theta_s = 1
            kr    = ones(size(u_exact));   % zona saturada (default): kr = 1 (K = Ksat)
            neg   = u_exact <= 0;          % mascara: pontos na zona nao saturada
            theta(neg) = (1 + (-alpha*u_exact(neg)).^nvg).^(-(nvg-1)/nvg);
            kr(neg)    = sqrt(theta(neg)) .* (1 - (1 - theta(neg).^(nvg/(nvg-1))).^((nvg-1)/nvg)).^2;
            kmap = env.config.perm(1,1) .* kr;   % K = Ksat * kr(u_exact)
        end

        % ── 13. Criterio de parada ────────────────────────────────
        % Este caso nao usa criterio de parada antecipada: a simulacao so
        % termina quando o criterio padrao "stopcriteria >= 100" (tempo
        % final atingido) e satisfeito em outra parte do codigo. Por isso
        % este metodo sempre retorna false.
        function parar = deveParar(obj, parms, premethod, stopcriteria)
            parar = false;
        end

        % ── 14. Atualiza flags no loop ────────────────────────────
        % As condicoes de contorno deste caso, embora dependam do tempo
        % (u_exact varia com t), sao recalculadas em outro ponto do fluxo
        % (via configurarContorno/configurarFlags chamadas diretamente a
        % cada passo pelo solver). Este flag apenas informa que NAO e
        % necessario reconstruir a estrutura de flags de contorno em si
        % (topologia de quais nos/faces sao Dirichlet/Neumann), que
        % permanece fixa ao longo da simulacao — apenas os VALORES
        % prescritos mudam a cada passo.
        function flag = precisaAtualizarFlags(obj, time)
            flag = false;
        end

        % ── 15. Finalizacao — relatorio de erro ────────────────────
        % Executado uma unica vez, ao final do loop temporal. Calcula e
        % imprime duas metricas de erro total (integradas no tempo),
        % combinando o erro de pressao com o erro de velocidade segundo
        % duas convencoes diferentes encontradas na literatura (ver
        % referencias nos comentarios do codigo original).
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

            % u_exact=options.extras.u_exact;
            % theta_exact=options.extras.theta_exact;
            % k_theta_exact=options.extras.k_theta_exact;
            % theta_n=options.extras.theta_n;
            % k_theta=options.extras.k_theta;
            % h=options.extras.h;
            % filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic406';
            % fname = fullfile(filepath);
            % celulas = readmatrix(fullfile(fname, 'elem_id_diagonal_edwards_quad_distorted_32.xlsx'));
            % 
            % 
            % filepath='C:\Users\flc59\Documents\Benchmark_Cases\BenchHydraulic406\teste_MPFAD_quad_32_k_theta';
            % fname = fullfile(filepath);
            % k_theta_exact = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_32_k_theta_1condhydraulic_analitica_steptime3.txt'));
            % k_theta = readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_32_k_theta_1condhydraulic_steptime3.txt'));
            % theta_n= readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_32_k_theta_1WaterContent_steptime3.txt'));
            % theta_exact= readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_32_k_theta_1WaterContent_analitica_steptime3.txt'));
            % u_exact=readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_32_k_theta_1exact_solution_steptime3.txt'));
            % h=readmatrix(fullfile(fname, 'Tables_teste_MPFAD_quad_32_k_theta_1h_steptime3.txt'));
            
            % 
            % figure(1)
            % plot(theta_exact(celulas, end),k_theta_exact(celulas,end));
            % hold on
            % plot(theta_n(celulas,end),k_theta(celulas,end));
            % 
            % figure(2)
            % plot(u_exact(celulas,end),theta_exact(celulas,end));
            % hold on
            % plot(h(celulas,end),theta_n(celulas,end));



            if ~isempty(options.extras)
                %fprintf('Erro L2 pressao (integrado no tempo): %.6e\n', sqrt(options.extras.hnsum));
                %fprintf('Erro L2 velocidade (integrado no tempo): %.6e\n', sqrt(options.extras.vnsum));
                elemarea=env.geometry.elemarea;
                wA = elemarea;                         % pesos por elemento (area)
                % Erro de pressao do ULTIMO passo de tempo apenas (nao
                % integrado no tempo, diferente de extras.hnsum):
                diffp = options.extras.u_exact-options.extras.h;                   % diferenca analitica - numerica no ultimo passo
                hnsum1=(sum(wA .* diffp.^2) / sum(wA) );
                % errortotal1: soma do erro de pressao JA integrado no
                % tempo (extras.hnsum, acumulado passo a passo) com o erro
                % de velocidade integrado no tempo (extras.vnsum).
                % Convencao de erro "L1" segundo o artigo citado:
                % "Analysis of an Euler Implicit-Mixed Finite Element
                % Scheme for Reactive Solute Transport in Porous Media"
                % (equacao 5.5, considerando dh=dt):
                errortotal1=(options.extras.hnsum+options.extras.vnsum);
                % errortotal2: usa em vez disso apenas o erro de pressao do
                % ULTIMO passo (hnsum1, sem integracao no tempo) somado ao
                % erro de velocidade integrado no tempo. Convencao de erro
                % "L2" segundo o artigo: "Convergence analysis for a mixed
                % finite element scheme for flow in strictly unsaturated
                % porous media" (ultimo paragrafo da equacao 37, dh=dt):
                errortotal2=(hnsum1+options.extras.vnsum);
                fprintf('Erro total L1 (integrado no tempo): %.6e\n', errortotal1);
                fprintf('Erro total L2  (integrado no tempo): %.6e\n', errortotal2);
            end
            if isfield(options.extras,'MBE')
                fprintf('MBE global (soma sobre os passos): %.6e\n', sum(options.extras.MBE));
                fprintf('MBE_rel maximo entre os passos:    %.6e\n', max(options.extras.MBE_rel));
                fprintf('max|Rloc| (pior elemento, todos os passos): %.6e\n', max(options.extras.maxRloc));
                fprintf('L1_MBE maximo entre os passos (intensivo):  %.6e\n', max(options.extras.L1_MBE));
                fprintf('L2_MBE maximo entre os passos (intensivo):  %.6e\n', max(options.extras.L2_MBE));

                % ── Figura para o artigo: MBE ao longo do tempo ──────────
                % indice 1 = t=0 (nao usado, ver "inicializar"); descarta na plotagem
                tt   = options.extras.time_storage(2:end);
                mbe  = abs(options.extras.MBE(2:end));
                mrel = options.extras.MBE_rel(2:end);
                rloc = options.extras.maxRloc(2:end);

                figure('Name','Balanco de massa (Caso436)');
                subplot(2,1,1)
                semilogy(tt, mbe, '-o', 'LineWidth', 1.2); hold on
                semilogy(tt, mrel, '-s', 'LineWidth', 1.2);
                xlabel('Tempo'); ylabel('Erro de balanco de massa');
                legend('|MBE| (absoluto)', 'MBE_{rel}', 'Location', 'best');
                grid on; title('MBE global por passo de tempo');

                subplot(2,1,2)
                semilogy(tt, rloc, '-^', 'LineWidth', 1.2, 'Color', [0.6 0 0]);
                xlabel('Tempo'); ylabel('max_i |R_i|');
                grid on; title('Pior residuo local por elemento (consistencia MPFA-D)');

                % ── salva a figura (mesma pasta/prefixo das tabelas) ──────
                filepath  = env.mainpathfolders.path;
                tabfolder = env.mainpathfolders.tabfolder;
                fname     = fullfile(filepath, tabfolder);
                savefig(gcf,  [fname 'MBE_figure.fig']);
                saveas(gcf,   [fname 'MBE_figure.png']);
            end
        end

        % ── 16. Escrita de resultados em arquivo ──────────────────
        % Salva os campos h, theta, kmap (numericos e seus equivalentes
        % analiticos/exatos) e os centroides em arquivos .txt (formato
        % CSV via writematrix), para pos-processamento externo (ex.:
        % scripts em Python/MATLAB para gerar graficos de convergencia).
        function escreverResultados(obj, env, h_storage, theta_storage, ...
                kmap_storage, time_storage, centelem, extras)
            filepath  = env.mainpathfolders.path;       % pasta raiz de saida da simulacao
            tabfolder = env.mainpathfolders.tabfolder;   % subpasta destinada a tabelas/arquivos de resultado
            fname = fullfile(filepath, tabfolder);
            writematrix(h_storage,     [fname 'h_steptime3.txt']);                         % pressao numerica ao longo do tempo
            writematrix(theta_storage, [fname 'WaterContent_steptime3.txt']);               % conteudo de agua numerico ao longo do tempo
            writematrix(centelem,      [fname 'centrocell3.txt']);                          % coordenadas dos centroides da malha
            writematrix(time_storage,  [fname 'time_step3.txt']);                           % instantes de tempo salvos
            writematrix(kmap_storage,  [fname 'condhydraulic_steptime3.txt']);              % condutividade hidraulica numerica ao longo do tempo
            writematrix(extras.exact_solution_storage,  [fname 'exact_solution_steptime3.txt']);        % solucao exata (manufaturada) ao longo do tempo
            writematrix(extras.theta_storage_analitica, [fname 'WaterContent_analitica_steptime3.txt']); % theta exato ao longo do tempo
            writematrix(extras.kmap_storage_analitica,  [fname 'condhydraulic_analitica_steptime3.txt']); % condutividade hidraulica exata ao longo do tempo

            % ── MBE por passo de tempo (para reuso/artigo, com cabecalho) ──
            % indice 1 = t=0 (nao usado, ver "inicializar") e descartado aqui
            if isfield(extras,'MBE')
                n = numel(extras.MBE);
                mbeTable = table( extras.time_storage(2:n)', extras.MBE(2:n)', ...
                    extras.MBE_rel(2:n)', extras.maxRloc(2:n)', ...
                    extras.L1_MBE(2:n)', extras.L2_MBE(2:n)', ...
                    'VariableNames', {'time','MBE','MBE_rel','maxRloc','L1_MBE','L2_MBE'});
                writetable(mbeTable, [fname 'MBE_steptime3.csv']);
            end
        end
    end

    methods(Static)

        function kappa = calcularKappa(theta, kappa_abs, mu, nvG)
            % CALCULARKAPPA  Condutividade hidraulica relativa via van Genuchten-Mualem
            %   kappa(theta) = (kappa_abs/mu) * sqrt(theta) .* ...
            %                  (1 - (1 - theta.^(nvG/(nvG-1))).^((nvG-1)/nvG)).^2
            %
            %   theta      : vetor (Nx1 ou 1xN) com o conteudo volumetrico de agua
            %                efetivo, valores em [0,1]
            %   kappa_abs  : permeabilidade absoluta (escalar)
            %   mu         : viscosidade do fluido (escalar)
            %   nvG        : parametro n do modelo de van Genuchten (escalar)

            m = (nvG - 1) / nvG;                       % expoente m = (nvG-1)/nvG
            expoente = nvG / (nvG - 1);                 % expoente nvG/(nvG-1)

            kappa = (kappa_abs / mu) .* sqrt(theta) .* (1 - (1 - theta.^expoente).^m).^2;

        end
    end
end