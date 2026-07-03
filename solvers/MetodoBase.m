%--------------------------------------------------------------------------
% MetodoBase — Classe abstrata base para todos os metodos numericos
%
% Esta classe define o CONTRATO que todo metodo numerico deve cumprir.
% Cada valor de pmethod no Start.dat corresponde a uma subclasse concreta:
%
%   pmethod = 'tpfa'   → MetodoTPFA
%   pmethod = 'mpfad'  → MetodoMPFAD
%   pmethod = 'mpfah'  → MetodoMPFAH
%   pmethod = 'nlfvpp' → MetodoNLFVPP
%   pmethod = 'mpfaql' → MetodoMPFAQL
%
% ARQUITETURA:
%   O metodo numerico e responsavel pela discretizacao espacial do fluxo.
%   Ele NAO sabe qual caso fisico esta sendo simulado (isso e responsabilidade
%   do benchmark). A separacao e clara:
%
%     env.benchmark → fisica do problema  (K, theta, BC, parametros)
%     env.metodo    → metodo numerico     (transmissibilidades, montagem, solver)
%
% COMO CRIAR UM NOVO METODO NUMERICO:
%   1. Crie o arquivo metodos/MetodoXXX.m
%   2. Declare: classdef MetodoXXX < MetodoBase
%   3. Implemente TODOS os metodos Abstract listados abaixo
%   4. Adicione uma linha em createMetodo.m: case 'xxx', metodo = MetodoXXX();
%   5. O metodo concreto calcGravidade() e herdado automaticamente
%
% METODOS ABSTRACT (obrigatorios em toda subclasse):
%   preprocessar, atualizarPremethod, montarSistema, resolver, calcularFlowrate
%
% METODOS CONCRETOS (herdados, geralmente nao precisam ser sobrescritos):
%   calcGravidade, exibir
%--------------------------------------------------------------------------
classdef MetodoBase < handle

    properties (Abstract)
        Nome      % string descritiva do metodo (ex: 'MPFA-D', 'TPFA')
                  % exibida no console durante a execucao
        MetodoID  % identificador string lido do Start.dat (ex: 'mpfad', 'tpfa')
                  % deve coincidir EXATAMENTE com o valor em pmethod
    end

    %----------------------------------------------------------------------
    % BLOCO ABSTRACT — toda subclasse DEVE implementar estes metodos
    % Se algum estiver faltando, o MATLAB lanca erro ao instanciar a classe
    % Dica de diagnostico: meta.class.fromName('MetodoXXX').MethodList
    %----------------------------------------------------------------------
    methods (Abstract)

        % preprocessar — calcula todos os parametros geometrico-fisicos do metodo
        %
        % Chamado UMA vez em preprocessmethod(), apos ferncodes_calflag().
        % E o momento de calcular tudo que depende da geometria e do kmap atual
        % mas NAO muda a cada passo de tempo (ex: V, N, F, Kde, pesos LPEW2).
        %
        % Deve popular env.premethod.XXXX.* com os campos necessarios.
        % Exemplo para MPFA-D (env.premethod.MPFAD.*):
        %   V, N, F      → mapeamentos geometricos (ferncodes_elementface)
        %   Hesq         → altura ortogonal centroide-aresta
        %   Kde          → transmissibilidade normal das arestas internas
        %   Ded          → correcao de anisotropia (desvio do gradiente)
        %   Kt           → permeabilidade tangencial (contorno)
        %   Kn           → permeabilidade normal (contorno)
        %   weight       → pesos de interpolacao LPEW2 (nos internos)
        %   s            → termos de correcao Neumann (nos de contorno)
        %   flowrateZ    → fluxo gravitacional por face
        %   flowresultZ  → fluxo gravitacional acumulado por elemento
        %
        % Tambem deve calcular env.preGravity (via obj.calcGravidade)
        % e, se necessario, pre-processar a concentracao acoplada.
        %
        % Entradas:
        %   env   → estado atual (geometria, config, kmap, nflag, nflagface)
        %   parms → parametros fisicos (h_old, auxperm, alpha, pp, q...)
        % Saida:
        %   env   → com env.premethod.XXXX.* e env.preGravity populados
        [env, parms] = preprocessar(obj, env, parms)

        % atualizarPremethod — recalcula parametros quando kmap muda no loop
        %
        % Chamado em hydraulic_RE a CADA passo de tempo, logo apos
        % PLUG_kfunction() atualizar kmap com o novo h convergido.
        %
        % Para Richards nao-linear, K(h) muda a cada passo de tempo,
        % entao as transmissibilidades (Kde, Ded, Kn, Kt) e os pesos
        % de interpolacao precisam ser recalculados com o novo kmap.
        %
        % IMPORTANTE: NAO recalcule V, N, F (geometria pura, invariante)
        % Apenas recalcule o que depende de kmap:
        %   MPFA-D → ferncodes_Kde_Ded_Kt_Kn + ferncodes_Pre_LPEW_2_vect
        %   TPFA   → ferncodes_Kde_Ded_Kt_Kn_TPFA (apenas Kde)
        %
        % Entradas:
        %   env   → estado atual com kmap ja atualizado
        %   parms → parametros fisicos com h_old atualizado (novo chute)
        % Saida:
        %   env   → com env.premethod.XXXX.* recalculado
        [env] = atualizarPremethod(obj, env, parms)

        % montarSistema — monta a matriz global A e o vetor independente b
        %
        % Chamado dentro de ferncodes_solver(), antes do iterador de Picard.
        % Deve montar o sistema linear A*h = b usando os parametros pre-calculados
        % em env.premethod.XXXX.* e as condicoes de contorno em env.config.nflag.
        %
        % A escolha da rotina de montagem depende do tipo de problema:
        %   Richards (400-500) → usa globalmatrix_XXXX com soil_properties
        %   Outros             → usa ferncodes_globalmatrix sem termo temporal
        %
        % Argumentos:
        %   env   → estado atual (premethod, nflag, nflagface, benchmark...)
        %   parms → parametros fisicos (h_old, auxperm, dt...)
        %   dt    → passo de tempo (usado apenas em problemas transientes)
        % Saidas:
        %   M     → matriz esparsa do sistema (nelem x nelem)
        %   I     → vetor independente (nelem x 1)
        [M, I] = montarSistema(obj, env, parms, dt)

        % resolver — resolve o sistema linear A*h = b com iteracao nao-linear
        %
        % Chamado dentro de ferncodes_solver(), apos montarSistema().
        % Deve escolher o iterador de acordo com env.config.acel:
        %
        %   'FPI'     → Picard classico (ferncodes_iterpicard)
        %               convergencia linear, robusto para Richards
        %               recomendado para a maioria dos casos
        %
        %   'AA'      → Picard com Aceleracao de Anderson
        %               convergencia superlinear, menos iteracoes de Picard
        %               indicado quando FPI converge lentamente
        %
        %   'LSCHEME' → L-scheme (regularizacao diagonal)
        %               incondicionalamente estavel mas mais iteracoes
        %               indicado para solos muito secos (Richards extremo)
        %
        % Dentro de cada iterador, a cada iteracao de Picard:
        %   1. Resolve h = A\b
        %   2. Atualiza kmap via PLUG_kfunction
        %   3. Recalcula transmissibilidades via atualizarPremethod
        %   4. Remonta A e b via montarSistema
        %   5. Verifica convergencia pelo residuo relativo ||A*h-b||/||r0||
        %
        % Apos convergencia, calcula flowrate via calcularFlowrate.
        %
        % Argumentos:
        %   M, I          → sistema linear pre-montado
        %   parms         → parametros fisicos (h_old, dt...)
        %   env           → estado atual
        %   tempo         → tempo atual (para BC dependentes do tempo)
        %   dt            → passo de tempo
        %   source_wells  → pocos injetores/produtores
        % Saidas:
        %   p             → campo de pressao/carga hidraulica convergido
        %   flowrate      → fluxo normal em cada face (nb+ni x 1)
        %   flowresult    → fluxo acumulado por elemento (nelem x 1)
        %   flowratedif   → fluxo difusivo (sem contribuicao gravitacional)
        %   faceaux       → indices auxiliares de faces
        %   parms         → atualizado (h_old = p convergido)
        %   env           → atualizado (premethod com novo kmap)
        [p, flowrate, flowresult, flowratedif, faceaux, ...
            parms, env] = resolver(obj, M, I, parms, env, tempo, dt, source_wells)

        % calcularFlowrate — calcula o fluxo normal nas faces apos convergencia
        %
        % Chamado dentro do iterador (ferncodes_iterpicard) com o campo
        % de pressao p ja convergido.
        %
        % A formula depende do metodo:
        %   MPFA-D → interpola p nos vertices (LPEW2) e usa formula diamante
        %            q = Kde*(p_R - p_L - Ded*(z2-z1)) + contribuicao vertices
        %   TPFA   → formula de dois pontos direta
        %            q = Kde*(h_R - h_L)
        %   NLFVPP → interpola com pesos nao-lineares
        %
        % Saidas:
        %   flowrate    → fluxo normal por face (nb+ni x 1) [m³/s ou m/s]
        %   flowresult  → flowrate acumulado por elemento (conservacao de massa)
        %   flowratedif → contribuicao difusiva isolada
        %   faceaux     → indices auxiliares
        [flowrate, flowresult, flowratedif, faceaux] = ...
            calcularFlowrate(obj, p, env, parms)

    end

    %----------------------------------------------------------------------
    % BLOCO CONCRETO — metodos compartilhados por todos os metodos numericos
    % Geralmente nao precisam ser sobrescritos nas subclasses
    %----------------------------------------------------------------------
    methods

        % exibir — imprime identificacao do metodo no console
        % Chamado no main apos createMetodo()
        function exibir(obj)
            fprintf('[Metodo] %s\n', obj.Nome);
        end

        % calcGravidade — calcula a taxa gravitacional para cada face
        %
        % Chamado no final de preprocessar() de cada metodo.
        % Necessario para simulacoes com gradiente gravitacional
        % (escoamento vertical, Richards, groundwater com gravidade).
        %
        % Retorna preGravity = [] se keygravity = 'n' (sem gravidade).
        % Retorna preGravity.gravrate se keygravity = 'y'.
        %
        % Atualmente ativo apenas para numcase < 300 (monofasico e concentracao).
        % Para Richards (400-500), a gravidade e tratada diretamente em
        % flowrateZ dentro de ferncodes_Kde_Ded_Kt_Kn.
        %
        % Se precisar ativar para outros intervalos de numcase,
        % adicione as condicoes no bloco if abaixo.
        function preGravity = calcGravidade(obj, env)
            preGravity = [];   % padrao: sem gravidade

            if strcmp(env.config.keygravity,'y')
                if env.config.numcase < 300
                    % calcula vetores de gravidade nos elementos e faces
                    [vec_gravelem, vec_gravface] = PLUG_Gfunction;

                    % calcula a taxa gravitacional por face usando kmap atual
                    % gravrate(i) = K(i) * g * cos(theta(i))
                    % onde theta e o angulo entre a face e a vertical
                    [gravrate] = gravitation(env.config.kmap, vec_gravelem, vec_gravface);

                    preGravity.gravrate = gravrate;
                end
            end
        end

    end
end