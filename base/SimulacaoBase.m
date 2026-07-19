%--------------------------------------------------------------------------
% SimulacaoBase — Classe abstrata base para todos os benchmarks
%
% Esta classe define o CONTRATO que todo benchmark deve cumprir.
% Cada numcase do Start.dat corresponde a uma subclasse concreta
% (ex: Caso439, Caso437, Caso341) que herda desta classe.
%
% ARQUITETURA:
%   O simulador tem tres eixos independentes de variacao:
%     env.benchmark = CasoXXX()    ← fisica do problema  (este contrato)
%     env.metodo    = MetodoXXX()  ← metodo numerico      (MetodoBase)
%     env.sim       = SimXXX()     ← tipo de simulacao    (SimulacaoBase)
%
% COMO CRIAR UM NOVO BENCHMARK:
%   1. Crie o arquivo benchmarks/CasoNNN.m
%   2. Declare: classdef CasoNNN < SimulacaoBase
%   3. Implemente TODOS os metodos Abstract listados abaixo
%   4. Adicione uma linha em createBenchmark.m: case NNN, bench = CasoNNN();
%   5. Os metodos concretos (calcularTheta, iso, etc.) sao herdados
%      automaticamente — sobrescreva apenas os que precisarem de
%      comportamento especifico para o seu caso
%
% METODOS ABSTRACT (obrigatorios em toda subclasse):
%   preprocessar, definirFontes, configurarPermeabilidade,
%   configurarContorno, configurarFlags, inicializar,
%   atualizarEstado, finalizar, escreverResultados
%
% METODOS CONCRETOS (herdados, sobrescreva se necessario):
%   calcularTheta, calcularCapacidade, iso, isoConst,
%   ajustarKContorno, temFlowrateBoundary, ajustarFlowrate,
%   calcularNeumannBoundary, adicionarTermoTemporal,
%   calcularTermoNeumannVet, prepararIndices, initParms
%--------------------------------------------------------------------------
classdef SimulacaoBase < handle

    properties (Abstract)
        Nome    % string descritiva do caso (ex: 'Processo de Recarga')
        TipoID  % numcase numerico (ex: 439) — deve coincidir com Start.dat
    end

    %----------------------------------------------------------------------
    % BLOCO ABSTRACT — toda subclasse DEVE implementar estes metodos
    % Se algum estiver faltando, o MATLAB lanca erro ao instanciar a classe
    % Dica de diagnostico: meta.class.fromName('CasoNNN').MethodList
    %----------------------------------------------------------------------
    methods (Abstract)

        %% ── Fluxo principal ──────────────────────────────────────

        % preprocessar — inicializa os parametros fisicos do caso
        %
        % Chamado por preRE() antes do loop temporal.
        % Deve preencher parms com todos os campos necessarios:
        %   parms.theta_s, theta_r  → parametros de retencao hidrica
        %   parms.alpha, pp, q, nvg → parametros de Van Genuchten / Brooks-Corey
        %   parms.h_init            → condicao inicial de h (vetor nelem x 1)
        %   parms.h_old             → chute inicial para Picard
        %   parms.dt                → passo de tempo
        %
        % Entrada:  env   (geometria e config ja carregadas)
        %           parms (struct vazio criado por initParms)
        % Saida:    parms preenchido, env (pode ser modificado se necessario)
        [parms, env] = preprocessar(obj, env, parms)

        % definirFontes — define pocos injetores e produtores
        %
        % Normalmente delega para defineWells(env, parms).
        % Sobrescreva apenas se o caso tiver logica especial de pocos.
        wells = definirFontes(obj, env, parms)

        %% ── Permeabilidade e contorno ────────────────────────────

        % configurarPermeabilidade — calcula kmap a cada passo de tempo
        %
        % Chamado por PLUG_kfunction() antes do loop e a cada iteracao.
        % Deve calcular o tensor de permeabilidade kmap usando h_old:
        %   kr = kr(h_old)             (permeabilidade relativa)
        %   kmap = Ks * kr             (tensor isotropico: use obj.iso())
        %   kmap = [K11 K12 K21 K22]   (tensor anisotropico: construa manualmente)
        %
        % Deve atualizar:
        %   env.config.kmap    → tensor atual (usado por Kde_Ded_Kt_Kn)
        %   parms.auxperm      → copia de kmap (usada por soil_properties)
        %   env.config.auxkmap → tensor saturado de referencia
        %   env.geometry.elem(:,5) → id de material por elemento
        [env, parms] = configurarPermeabilidade(obj, env, parms, time)

        % configurarContorno — calcula o valor de h prescrito nas faces/vertices
        %
        % Chamado por configurarFlags() e por PLUG_bcfunction().
        % Recebe os indices de vertices/faces de contorno e retorna
        % o valor prescrito de h (ou concentracao) em cada um.
        %
        % Argumentos:
        %   vertices  → indices dos vertices ou [v1 v2] das arestas
        %   flagptr   → flags de BC (ex: 1=Dirichlet entrada, 2=Dirichlet saida)
        %   time      → tempo atual (para BC dependentes do tempo)
        %   env, parms → estado atual da simulacao
        %
        % O calculo do ponto medio (coordmid) deve ser feito internamente:
        %   if size(vertices,2) > 1  → aresta: coordmid = media dos dois vertices
        %   else                     → vertice: coordmid = coord do vertice
        bcattrib = configurarContorno(obj, vertices, flagptr, time, env, parms)

        % configurarFlags — monta nflag e nflagface
        %
        % Chamado por ferncodes_calflag() uma vez antes do loop temporal
        % (e a cada passo se precisaAtualizarFlags() retornar true).
        %
        % Deve preencher:
        %   nflag(:,1)     → tipo de BC no vertice (< 200: Dirichlet, >= 200: Neumann)
        %   nflag(:,2)     → valor prescrito de h no vertice
        %   nflagface(:,1) → tipo de BC na face
        %   nflagface(:,2) → valor prescrito de h na face
        %
        % Use obj.prepararIndices(env) para obter os mapeamentos
        % e obj.configurarContorno() para calcular os valores
        [nflag, nflagface] = configurarFlags(obj, env, parms, time)

        %% ── Loop temporal ────────────────────────────────────────

        % inicializar — executado UMA vez antes do while em hydraulic_RE
        %
        % Use para:
        %   - Localizar elementos de monitoramento (find com mascaras)
        %   - Salvar indices em .mat para reutilizar em proximas simulacoes
        %   - Inicializar series temporais (h_time, theta_time, MBE...)
        %   - Chamar postprocessor para o instante t=0
        %   - Calcular theta/massa inicial para calculo de erros
        %
        % Saida: extras — struct livre para armazenar qualquer dado
        %        que precise persistir entre os passos de tempo
        [parms, extras] = inicializar(obj, env, parms, time)

        % atualizarEstado — executado a CADA passo de tempo no while
        %
        % Use para:
        %   - Atualizar parms.h_old (chute inicial para proximo Picard)
        %   - Atualizar parms.h_init (condicao inicial para proximo dt)
        %   - Chamar postprocessor para salvar VTK
        %   - Armazenar h e theta nos pontos de monitoramento
        %   - Calcular erros L2, MBE, fluxos de contorno
        %   - Atualizar source terms dependentes de h ou t (caso 436)
        %
        % ATENCAO: theta_n NAO e passado como argumento aqui.
        % Se precisar de theta, chame: thetafunction(h, parms, env)
        % ou env.benchmark.calcularTheta(h, parms)
        [parms, extras] = atualizarEstado(obj, env, parms, extras, ...
            h, flowrate, time, count)

        % finalizar — executado UMA vez apos o while em hydraulic_RE
        %
        % Use para:
        %   - Plotar graficos (figure, plot, xlabel, grid...)
        %   - Calcular e exibir erros finais (L2, H1, MBE)
        %   - Comparar com solucao analitica
        %
        % ATENCAO: a assinatura atual e (obj, env, extras, theta_n)
        % Se precisar de outros campos, adicione ao struct extras em inicializar()
        finalizar(obj, env, options)

        % escreverResultados — executado apos finalizar() em hydraulic_RE
        %
        % Use para salvar os campos em arquivo .txt ou .mat:
        %   h_storage     → campo de h em todos os passos de tempo
        %   theta_storage → campo de theta em todos os passos de tempo
        %   kmap_storage  → condutividade hidraulica em todos os passos
        %   time_storage  → vetor de instantes de tempo
        %
        % Exemplo de uso tipico:
        %   writematrix(h_storage, [fname 'h_steptime.txt']);
        escreverResultados(obj, env, h_storage, theta_storage, kmap_storage, ...
            time_storage, centelem, extras)

    end

    %----------------------------------------------------------------------
    % BLOCO CONCRETO — metodos compartilhados por todos os benchmarks
    % Nao e necessario implementar nas subclasses, a menos que o
    % comportamento especifico do caso seja diferente do padrao aqui definido
    %----------------------------------------------------------------------
    methods

        % exibir — imprime identificacao do benchmark no console
        % Chamado no main apos createBenchmark()
        function exibir(obj)
            fprintf('[Benchmark %g] %s\n', obj.TipoID, obj.Nome);
        end

        %% ── Modelos de retencao hidrica ──────────────────────────

        % calcularTheta — conteudo volumetrico de agua theta(h)
        %
        % FALLBACK: Van Genuchten standard (casos 431, 432, 433, 434, 435)
        %   theta = theta_s                                      se h >= 0
        %   theta = theta_r + (theta_s-theta_r)/(1+(-ah)^n)^m   se h <  0
        %
        % Sobrescreva em subclasses para outros modelos:
        %   Caso436 → Van Genuchten com nvg (expoente diferente)
        %   Caso437 → Gardner: theta = theta_r + (theta_s-theta_r)*exp(alpha*h)
        %   Caso438 → Cubico:  theta = (2-h)^(-1/3)
        %   Caso439 → Brooks-Corey: Se = c/(c+|h|^D)
        function theta = calcularTheta(obj, h, parms)
            theta_s = parms.theta_s;
            theta_r = parms.theta_r;
            alpha   = parms.alpha;
            pp      = parms.pp;
            q       = parms.q;
            theta   = theta_s * ones(size(h));
            idx_neg = h < 0;
            theta(idx_neg) = theta_r + (theta_s - theta_r) ./ ...
                (1 + (-alpha*h(idx_neg)).^pp).^q;
        end

        %% ── Helpers de permeabilidade ────────────────────────────

        % iso — monta tensor de permeabilidade isotropico [idx, K, 0, 0, K]
        %
        % Formato de kmap: [idx | K11 | K12 | K21 | K22]
        % Para caso isotropico: K11 = K22 = coef, K12 = K21 = 0
        %
        % Uso tipico em configurarPermeabilidade:
        %   coef = Ks .* kr;
        %   env.config.kmap = obj.iso(env, coef);
        function kmap = iso(obj, env, coef)
            idx   = env.utils.idx;
            nelem = env.utils.nelem;
            kmap  = [idx, coef, zeros(nelem,1), zeros(nelem,1), coef];
        end

        % calcSe — saturacao efetiva de Van Genuchten
        %   Se = 1 / (1 + (-alpha*h)^pp)^q   para h <= 0
        %   Se = 1                             para h >  0
        % Usado internamente por calcularTheta e calcularCapacidade
        function Se = calcSe(obj, h_old, alpha, pp, q, nelem)
            Se      = ones(nelem,1);
            neg     = h_old <= 0;
            Se(neg) = 1 ./ (1 + (-alpha*h_old(neg)).^pp).^q;
        end

        % calcTheta — theta de Van Genuchten com parametro nvg (caso 436)
        %   Theta = (1 + (-alpha*h)^n)^(-(n-1)/n)   para h <= 0
        %   Theta = 1                                  para h >  0
        function Theta = calcTheta(obj, h_old, alpha, n, nelem)
            Theta      = ones(nelem,1);
            neg        = h_old <= 0;
            Theta(neg) = (1 + (-alpha*h_old(neg)).^n).^(-(n-1)/n);
        end

        % isoConst — tensor isotropico com permeabilidade saturada constante
        %
        % Retorna kmap com K = perm(1,1) (condutividade hidraulica saturada)
        % em todos os elementos. Usado como referencia em auxkmap.
        function auxkmap = isoConst(obj, env)
            auxkmap = obj.iso(env, env.config.perm(1,1).*ones(env.utils.nelem,1));
        end

        % bcflagDefault — valor de BC diretamente de bcflag
        %
        % Fallback para casos com BC simples (valor fixo lido do Start.dat).
        % Chamado em configurarContorno quando nao ha formula especial.
        function bcattrib = bcflagDefault(obj, flagptr, env)
            bcattrib = env.config.bcflag(flagptr, 2);
        end

        %% ── Helpers de contorno ──────────────────────────────────

        % prepararIndices — mapeia bcflag para vertices e faces de contorno
        %
        % Retorna os indices e valores de BC para todos os vertices e faces
        % do contorno, mapeados a partir de env.config.bcflag.
        % Usado em configurarFlags() de todas as subclasses.
        %
        % Saidas:
        %   vertex_idx    → indices dos vertices de contorno em coord
        %   face_idx      → indices [v1 v2] das faces de contorno em bedge
        %   bcflag_vertex → linhas de bcflag correspondentes aos vertices
        %   bcflag_face   → linhas de bcflag correspondentes as faces
        %   bc_row_vertex → linha de bcflag para cada vertice (mapeamento)
        %   bc_row_face   → linha de bcflag para cada face (mapeamento)
        function [vertex_idx, face_idx, bcflag_vertex, bcflag_face, ...
                bc_row_vertex, bc_row_face] = prepararIndices(obj, env)
            vertex_flag = env.geometry.bedge(:,4);
            face_flag   = env.geometry.bedge(:,5);
            vertex_idx  = env.geometry.bedge(:,1);
            face_idx    = env.geometry.bedge(:,1:2);

            % para cada vertice/face, encontra a linha correspondente em bcflag
            x_vertex = env.config.bcflag(:,1) == vertex_flag';
            x_face   = env.config.bcflag(:,1) == face_flag';
            [~, bc_row_vertex] = max(x_vertex, [], 1);
            [~, bc_row_face]   = max(x_face,   [], 1);

            bcflag_vertex = env.config.bcflag(bc_row_vertex, :);
            bcflag_face   = env.config.bcflag(bc_row_face,   :);
        end

        % initParms — cria struct de parametros fisicos com valores padrao zero
        %
        % Chamado no main antes de sim.preprocessar().
        % Garante que todos os campos existem mesmo que o benchmark
        % nao os preencha (evita erros de campo inexistente).
        %
        % Campos:
        %   theta_s, theta_r → umidade de saturacao e residual
        %   alpha, pp, q     → parametros de Van Genuchten
        %   nvg              → expoente de Van Genuchten (caso 436)
        %   hs               → pressao de entrada (Brooks-Corey)
        %   valuecontor      → valor no contorno (caso 437)
        %   h_init, h_old    → condicao inicial e chute para Picard
        %   dt               → passo de tempo
        function parms = initParms(obj)
            parms.theta_s     = 0;
            parms.theta_r     = 0;
            parms.alpha       = 0;
            parms.pp          = 0;
            parms.q           = 0;
            parms.nvg         = 0;
            parms.hs          = 0;
            parms.valuecontor = 0;
            parms.h_init      = [];
            parms.h_old       = [];
            parms.dt          = 0;
        end

        %% ── Permeabilidade na fronteira ──────────────────────────

        % ajustarKContorno — ajusta K11..K22 nas faces de contorno Dirichlet
        %
        % FALLBACK: retorna K sem ajuste (permeabilidade base do elemento)
        % Usado em ferncodes_Kde_Ded_Kt_Kn e ferncodes_Kde_Ded_Kt_Kn_TPFA.
        %
        % Sobrescreva para modelos nao-lineares onde K depende de h_contorno:
        %   Caso437 → Gardner:      K = Ks * exp(alpha*h_contorno)
        %   Caso438 → Cubico:       K = Ks * (2-h)^(-1)
        %   Caso439 → Brooks-Corey: K = 35 * kr(h_contorno)
        %
        % Argumentos:
        %   auxkmap    → tensor de permeabilidade dos elementos (nelem x 5)
        %   matid      → id de material de cada face de contorno
        %   h_contorno → valor de h prescrito em cada face (nflagface(:,2))
        %   maskT      → mascara: true onde bedge(:,5) < 200 (Dirichlet)
        function [K11, K12, K21, K22] = ajustarKContorno(obj, env, parms, ...
                auxkmap, matid, h_contorno, maskT)
            K11 = auxkmap(matid,2);
            K12 = auxkmap(matid,3);
            K21 = auxkmap(matid,4);
            K22 = auxkmap(matid,5);
        end

        % temFlowrateBoundary — controla se o bloco de fluxo gravitacional
        % nas faces de contorno e calculado em Kde_Ded_Kt_Kn
        %
        % FALLBACK: false — nao calcula (casos lineares e estacionarios)
        %
        % Sobrescreva retornando true para casos Richards (431..439)
        % onde o fluxo gravitacional nas faces de contorno e necessario
        % para a montagem correta do vetor b (soil_properties)
        function flag = temFlowrateBoundary(obj)
            flag = false;
        end

        % ajustarFlowrate — pos-ajuste do flowrateZ apos calculo
        %
        % FALLBACK: nao faz nada
        %
        % Sobrescreva para casos que precisam de ajuste especial:
        %   Caso435 → inverte sinal em faces com flag 101 (fluxo de saida)
        function flowrateZ = ajustarFlowrate(obj, flowrateZ, bedge)
            % nao faz nada — fallback identico para maioria dos casos
        end

        %% ── Contorno de Neumann ──────────────────────────────────

        % calcularNeumannBoundary — calcula contribuicao de Neumann no vetor b
        %
        % FALLBACK: formula geral — fluxo prescrito * norma da aresta
        %   valsI_N = norN .* bcflag(loc,2) + flowrateZ(mask222)
        %
        % Sobrescreva para casos com Neumann especial:
        %   Caso341 → usa ferncodes_K(x,y) para permeabilidade variavel
        %
        % Chamado em ferncodes_globalmatrix_MPFAD e ferncodes_globalmatrix_TPFA
        function valsI_N = calcularNeumannBoundary(obj, isNeu, bedge, ...
                bcflag, nflagface, flowrateZ, nor_b, normals, env)
            flagN   = bedge(isNeu,5);
            norN    = nor_b(isNeu);
            [~,loc] = ismember(flagN, bcflag(:,1));
            mask222 = bedge(:,5) > 200;
            valsI_N = norN.*bcflag(loc,2) + flowrateZ(find(mask222),1);
        end

        % adicionarTermoTemporal — adiciona matriz de massa e vetor de acumulacao
        %
        % FALLBACK: nao faz nada (problemas estacionarios)
        %
        % Chamado em ferncodes_globalmatrix_MPFAD e ferncodes_globalmatrix_TPFA
        % apos montagem da parte difusiva.
        %
        % Sobrescreva para:
        %   Richards (431..439) → chama soil_properties (capacidade hidrica)
        %   Groundwater transiente → chama ferncodes_implicitandcranknicolson
        function [M,I] = adicionarTermoTemporal(obj, M, I, parms, flowresultZ, env)
            % nao faz nada por padrao
        end

        %% ── Interpolacao LPEW2 ───────────────────────────────────

        % calcularTermoNeumannVet — calcula o termo "s" para nos de contorno
        % com condicao de Neumann, usado na interpolacao LPEW2
        %
        % FALLBACK: formula geral (eq. 26 do artigo LPEW2 — Agelas et al. 2010)
        %   s(No) = -(1/sum_lambda) * (r1*flux1 + r2*flux2)
        %
        % Chamado UMA vez fora do loop sobre nos em ferncodes_Pre_LPEW_2_vect,
        % apos o calculo dos pesos weight para todos os nos.
        % Esta abordagem e mais eficiente que chamar dentro do loop
        % (elimina overhead de despacho OOP por no).
        %
        % Sobrescreva para casos com Neumann especial:
        %   Caso341 → usa ferncodes_K(x,y) para K variavel na fronteira
        %
        % Argumentos:
        %   r          → fatores geometricos de contorno (nNodes x 2)
        %   sum_lambda → soma dos lambdas por no (nNodes x 1)
        %   N          → mapeamento face→no (da ferncodes_elementface)
        %   env        → estado atual (geometria, config, bcflag...)
        function s = calcularTermoNeumannVet(obj, r, sum_lambda, N, env)
            nNodes = size(env.geometry.coord, 1);
            s      = zeros(nNodes, 1);
            bedge  = env.geometry.bedge;
            bcflag = env.config.bcflag;
            nb     = size(bedge, 1);
            ns2    = env.geometry.nsurn2;   % estrutura CSR: vizinhos por no

            for No = 1:nNodes
                nns   = ns2(No+1) - ns2(No);   % numero de nos vizinhos
                face1 = N(No,1);               % primeira face de contorno do no
                face2 = N(No,nns);             % ultima  face de contorno do no

                % no nao toca o contorno — pula
                if face1 > nb || face2 > nb, continue; end

                % verifica se o no pertence ao contorno (bedge(:,1) == No)
                MM  = bedge(:,1) == No;
                MMM = find(MM, 1);
                if isempty(MMM), continue; end

                % verifica se tem condicao de Neumann (flag > 200)
                if ~(200 < bedge(MMM,4)), continue; end

                % localiza o fluxo prescrito correspondente em bcflag
                a  = bcflag(:,1) == bedge(face1,5);
                b  = bcflag(:,1) == bedge(face2,5);
                s1 = find(a, 1);
                s2 = find(b, 1);
                if isempty(s1) || isempty(s2), continue; end

                % eq. 26 do artigo LPEW2 (Agelas et al., 2010)
                s(No) = -(1/sum_lambda(No)) * ...
                    (r(No,1)*bcflag(s1,2) + r(No,2)*bcflag(s2,2));
            end
        end

    end
end