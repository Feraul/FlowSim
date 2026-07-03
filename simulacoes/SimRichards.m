%--------------------------------------------------------------------------
% SimRichards — Tipo de simulacao: Equacao de Richards (phasekey = 6)
%
% Esta classe representa o TIPO de simulacao Richards, ou seja, fluxo
% em meio poroso parcialmente saturado (zona vadosa e zona saturada).
%
% PAPEL DESTA CLASSE:
%   SimRichards define o FLUXO de execucao da simulacao Richards:
%     - Como pre-processar (chama preRE)
%     - Como definir os pocos (chama defineWells)
%
%   A FISICA ESPECIFICA (K(h), theta(h), h_init, parametros) e
%   responsabilidade do benchmark (Caso439, Caso437, etc.).
%
% SEPARACAO DE RESPONSABILIDADES:
%
%   SimRichards (esta classe)        CasoXXX (benchmark)
%   ─────────────────────────────    ──────────────────────────────────
%   preprocessar → chama preRE       preprocessar → preenche theta_s, h_init, dt
%   definirFontes → chama defineWells configurarPermeabilidade → K(h) especifico
%   configurarFlags → inicializa vazio calcularTheta → modelo hidrico especifico
%   inicializar → struct vazio        inicializar → pontos de monitoramento
%   atualizarEstado → vazio           atualizarEstado → postprocessor, h_time
%   finalizar → vazio                 finalizar → graficos, erros L2
%   escreverResultados → vazio        escreverResultados → writematrix
%
% POR QUE OS METODOS ABAIXO SAO VAZIOS:
%   SimulacaoBase declara todos estes metodos como Abstract,
%   obrigando TODA subclasse a implementa-los.
%   Como a logica real esta no benchmark (CasoXXX), SimRichards
%   implementa os metodos com corpo vazio apenas para satisfazer
%   o contrato — o benchmark sobrescreve o que precisa.
%
% FLUXO DE EXECUCAO (main → hydraulic_RE):
%   main:
%     sim = createSimulacao(6)          → instancia SimRichards
%     [env, parms] = sim.preprocessar() → chama preRE → benchmark.preprocessar
%     wells = sim.definirFontes()       → chama defineWells
%   hydraulic_RE:
%     benchmark.inicializar()           → benchmark faz o que precisa
%     while loop:
%       benchmark.atualizarEstado()     → benchmark faz o que precisa
%     benchmark.finalizar()             → benchmark faz o que precisa
%     benchmark.escreverResultados()    → benchmark faz o que precisa
%--------------------------------------------------------------------------
classdef SimRichards < SimulacaoBase

    properties
        Nome   = 'Richards Equation'   % exibido no console pelo main
        TipoID = 6                      % corresponde ao phasekey no Start.dat
    end

    methods

        %% ── Responsabilidade REAL de SimRichards ─────────────────

        % preprocessar — inicializa o pre-processamento da equacao de Richards
        %
        % Chama preRE(env), que por sua vez:
        %   1. Cria parms vazio via benchmark.initParms()
        %   2. Chama benchmark.preprocessar(env, parms) que preenche:
        %      parms.theta_s, theta_r, alpha, pp, q (ou nvg)
        %      parms.h_init  → condicao inicial de h
        %      parms.h_old   → chute inicial para Picard
        %      parms.dt      → passo de tempo
        %
        % NOTA: a ordem de retorno e [env, parms] aqui,
        % mas preRE retorna [parms, env] internamente — verifique a assinatura
        % de preRE se houver erros de atribuicao
        function [env, parms] = preprocessar(obj, env, parms)
            [parms, env] = preRE(env);
        end

        % definirFontes — define pocos injetores e produtores
        %
        % Delega para defineWells(env, parms), que le as posicoes e
        % vazoes dos pocos do Start.dat e monta a estrutura source_wells.
        % O resultado e passado para setmethod → hydraulic_RE → ferncodes_solver.
        function wells = definirFontes(obj, env, parms)
            wells = defineWells(env, parms);
        end

        %% ── Metodos vazios — benchmark concreto sobrescreve ──────
        %
        % Os metodos abaixo existem APENAS para satisfazer o contrato
        % Abstract de SimulacaoBase. A logica real esta nos benchmarks.
        %
        % Se adicionar logica aqui, ela sera executada PARA TODOS os
        % casos Richards (431..439). Use apenas para comportamento
        % verdadeiramente comum a todos — caso contrario, implemente
        % no benchmark especifico.

        % configurarPermeabilidade — vazio
        % A permeabilidade e calculada por PLUG_kfunction via
        % benchmark.configurarPermeabilidade(). SimRichards nao interfere.
        function [env, parms] = configurarPermeabilidade(obj, env, parms, time)
            % vazio — PLUG_kfunction delega ao benchmark (CasoXXX)
        end

        % configurarContorno — fallback simples
        % Retorna o valor fixo do bcflag para cada flag de contorno.
        % Os benchmarks que tem BC especiais (caso 439: h=65-z, caso 436: BC(t))
        % sobrescrevem este metodo com a formula correta.
        function bcattrib = configurarContorno(obj, vertices, flagptr, time, env, parms)
            % fallback: valor direto de bcflag — benchmarks sobrescrevem
            bcattrib = env.config.bcflag(flagptr, 2);
        end

        % configurarFlags — inicializa com valores sentinela
        % O preenchimento real de nflag e nflagface e feito por
        % ferncodes_calflag(), que chama benchmark.configurarFlags().
        % SimRichards apenas retorna valores iniciais para evitar
        % erros de campo inexistente antes de ferncodes_calflag ser chamado.
        %
        % Valores sentinela:
        %   nflag = 5000     → indica "no interior" (nao processado ainda)
        %   nflagface = 0    → sem flag de contorno nas faces
        function [nflag, nflagface] = configurarFlags(obj, env, parms, time)
            % inicializa com sentinelas — ferncodes_calflag preenche os reais
            nflag     = 5000 * ones(size(env.geometry.coord,1), 2);
            nflagface = zeros(size(env.geometry.bedge,1), 2);
        end

        % inicializar — vazio
        % Chamado em hydraulic_RE antes do loop temporal.
        % O benchmark concreto sobrescreve para:
        %   - localizar elementos de monitoramento
        %   - inicializar series temporais
        %   - salvar indices em .mat
        %   - postprocessor t=0
        % SimRichards retorna extras vazio para que hydraulic_RE
        % possa passar extras para os metodos subsequentes.
        function [parms, extras] = inicializar(obj, env, parms, time)
            extras = struct();   % benchmark concreto preenche o que precisar
        end

        % atualizarEstado — vazio
        % Chamado em hydraulic_RE a cada passo de tempo.
        % O benchmark concreto sobrescreve para:
        %   - atualizar parms.h_old (chute para proximo Picard)
        %   - chamar postprocessor (salvar VTK)
        %   - armazenar h e theta nos pontos de monitoramento
        %   - calcular erros L2 ou MBE
        function [parms, extras] = atualizarEstado(obj, env, parms, extras, ...
                h, flowrate, time, count)
            % vazio — benchmark concreto (CasoXXX) implementa a logica
        end

        % finalizar — vazio
        % Chamado em hydraulic_RE apos o loop temporal.
        % O benchmark concreto sobrescreve para:
        %   - plotar graficos (h(t), theta(z), frentes de umidade)
        %   - calcular e exibir erros finais (L2, H1, MBE)
        %   - comparar com solucao analitica
        function finalizar(obj, env, parms, extras, h, theta_n, time)
            % vazio — benchmark concreto (CasoXXX) implementa a logica
        end

        % escreverResultados — vazio
        % Chamado em hydraulic_RE apos finalizar().
        % O benchmark concreto sobrescreve para salvar em arquivo:
        %   h_storage, theta_storage, kmap_storage, time_storage
        % Formato tipico: writematrix(dados, [fname 'arquivo.txt'])
        function escreverResultados(obj, env, h_storage, theta_storage, ...
                kmap_storage, time_storage, extras)
            % vazio — benchmark concreto (CasoXXX) implementa a logica
        end

    end
end