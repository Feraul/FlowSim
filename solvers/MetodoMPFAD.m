%--------------------------------------------------------------------------
% MetodoMPFAD — Metodo de Volumes Finitos com Patch Diamante (MPFA-D)
%
% Implementa o metodo MPFA-D (Multi-Point Flux Approximation — Diamond)
% para aproximacao do fluxo difusivo em malhas nao-estruturadas com
% tensores de permeabilidade heterogeneos e anisotropicos.
%
% Referencias:
%   Contreras et al. (2016, 2019, 2021) — MPFA-D e NL-TPFA
%   Agelas et al. (2010) — LPEW2 (pesos de interpolacao)
%
% Responsabilidades:
%   preprocessar()       → calcula V, N, F, Kde, Ded, Kt, Kn, pesos LPEW2
%   atualizarPremethod() → recalcula Kde e pesos quando kmap muda (Richards)
%   montarSistema()      → monta a matriz global A e vetor b
%   resolver()           → escolhe o iterador (Picard, AA, L-scheme) e resolve
%   calcularFlowrate()   → interpola pressao e calcula fluxo nas faces
%--------------------------------------------------------------------------
classdef MetodoMPFAD < MetodoBase

    properties
        Nome     = 'MPFA-D'    % nome para exibicao no console
        MetodoID = 'mpfad'     % identificador string (lido do Start.dat)
    end

    methods

        %% ── 1. Pre-processamento geometrico e fisico ─────────────
        % Chamado UMA vez antes do loop temporal (em preprocessmethod).
        % Calcula e armazena em env.premethod.MPFAD.*:
        %
        %   V, N, F  → mapeamentos geometricos de faces, vertices e elementos
        %              (usados na interpolacao LPEW2 e montagem da matriz)
        %
        %   Kde      → transmissibilidade normal da aresta interna
        %   Ded      → correcao de anisotropia (desvio do gradiente)
        %   Kt       → permeabilidade tangencial (contorno)
        %   Kn       → permeabilidade normal (contorno)
        %
        %   weight   → pesos de interpolacao LPEW2 nos vertices internos
        %   s        → termos de correcao Neumann nos vertices de contorno
        %
        %   flowrateZ   → fluxo gravitacional por face (Richards)
        %   flowresultZ → fluxo gravitacional acumulado por elemento
        function [env, parms] = preprocessar(obj, env, parms)

            % ── geometria das faces, vertices e elementos ─────────
            % V: vertices por face | N: faces por vertice | F: faces por elemento
            [V, N, F] = ferncodes_elementface(env);
            env.premethod.MPFAD.V = V;
            env.premethod.MPFAD.N = N;
            env.premethod.MPFAD.F = F;

            % ── transmissibilidades Kde, Ded, Kt, Kn ──────────────
            % Calcula os parametros de transmissibilidade para todas as
            % arestas (internas e de contorno) com o kmap atual
            [env] = ferncodes_Kde_Ded_Kt_Kn(env, parms);

            % ── pesos de interpolacao LPEW2 ───────────────────────
            % Calcula os pesos de interpolacao da pressao nos vertices
            % internos (weight) e os termos de correcao Neumann (s)
            % Metodo: LPEW2 (Linearity Preserving — Edge Weighting 2)
            [env, ~, ~] = ferncodes_Pre_LPEW_2_vect(env, parms);

            % ── pre-processamento de concentracao acoplada ─────────
            % Apenas para simulacoes de transporte (casos 200-300 e 350-400)
            % Calcula Kdec, Knc, pesos e flags para a equacao de concentracao
            nc = env.config.numcase;
            if (200 < nc && nc < 300) || (350 < nc && nc < 400)
                env.premethod.conMPFAD = obj.preprocessConcentracao(env, N);
            end

            % ── gravidade ─────────────────────────────────────────
            % Calcula gravrate (taxa gravitacional por face) se keygravity='y'
            env.preGravity = obj.calcGravidade(env);
        end

        %% ── 2. Atualiza premethod no loop temporal ───────────────
        % Chamado em hydraulic_RE a CADA passo de tempo, apos PLUG_kfunction
        % atualizar kmap com o novo h.
        %
        % Necessario para Richards nao-linear: K(h) muda a cada iteracao,
        % entao Kde, Ded, Kn, Kt e os pesos LPEW2 precisam ser recalculados.
        %
        % NAO recalcula V, N, F (geometria pura — invariante no tempo)
        function [env] = atualizarPremethod(obj, env, parms)
            [env]       = ferncodes_Kde_Ded_Kt_Kn(env, parms);
            %[env, ~, ~] = ferncodes_Pre_LPEW_2_vect(env, parms);
            [env,~,~] = ferncodes_Pre_LPEW_2_vect_antigo(env,parms);

        end

        %% ── 3. Monta a matriz global do sistema linear ───────────
        % Seleciona a rotina de montagem adequada ao tipo de problema:
        %
        %   numcase == 331 ou 400-500 (Richards):
        %     → ferncodes_globalmatrix_MPFAD: monta A com termo temporal
        %       (capacidade hidrica via soil_properties)
        %
        %   outros (estacionario ou groundwater):
        %     → ferncodes_globalmatrix: monta A sem termo temporal
        %       (usa premethod.MPFAD diretamente)
        function [M, I] = montarSistema(obj, env, parms, dt)
            if env.config.numcase == 331 || ...
                    (400 < env.config.numcase && env.config.numcase < 500)
                % Richards / fluxo nao-saturado — inclui termo dtheta/dt
                [M, I] = ferncodes_globalmatrix_MPFAD(env, parms);
            else
                % estacionario ou groundwater — sem termo temporal
                [M, I, ~] = ferncodes_globalmatrix(env, env.premethod.MPFAD, parms, dt);
            end
        end

        %% ── 4. Resolve o sistema linear A*h = b ──────────────────
        % Seleciona o iterador de acordo com env.config.acel:
        %
        %   'FPI'     → Picard classico (Fixed Point Iteration)
        %               convergencia linear, robusto para Richards
        %
        %   'AA'      → Picard com Aceleracao de Anderson
        %               convergencia superlinear, menos iteracoes
        %
        %   'LSCHEME' → L-scheme (estabilizacao por termo de regularizacao)
        %               incondicionalalmente estavel, mais iteracoes
        %               indicado para solos muito secos (alta nao-linearidade)
        function [p, flowrate, flowresult, flowratedif, faceaux, parms, env] = ...
                resolver(obj, M, I, parms, env, tempo, dt, source_wells)

            switch upper(env.config.acel)

                case 'FPI'
                    % Picard classico — mais comum para Richards
                    [p, flowrate, flowresult, flowratedif, faceaux, parms, env] = ...
                        ferncodes_iterpicard(M, I, parms, env, tempo, dt, source_wells);

                case 'AA'
                    % Aceleracao de Anderson — reduz numero de iteracoes de Picard
                    [p, flowrate, flowresult, flowratedif, faceaux, parms] = ...
                        ferncodes_iterpicardANLFVPP2(M, I, env.premethod.MPFAD, ...
                        parms, env, tempo, dt, source_wells);

                case 'LSCHEME'
                    % L-scheme — adiciona L*I a diagonal para garantir estabilidade
                    [p, flowrate, flowresult, flowratedif, faceaux, parms, ...
                        env.premethod.MPFAD] = ...
                        L_scheme(M, I, env.premethod.MPFAD, parms, env, ...
                        tempo, dt, source_wells);
                case 'LNEAR'
                    p = M \ I;
                    % auxiliary variables interpolation
                    [pinterp,~]=ferncodes_pressureinterpNLFVPP(p,env);
                    %Get the flow rate (Diamond)
                    [flowrate,flowresult,flowratedif,faceaux] = ferncodes_flowrate(p,pinterp,...
                        env);

            end
        end

        %% ── 5. Calcula o flowrate (fluxo normal nas faces) ───────
        % Apos convergencia do Picard, calcula o fluxo normal discreto:
        %   1. Interpola a pressao nos vertices (ferncodes_pressureinterpNLFVPP)
        %      usando os pesos LPEW2 pre-calculados
        %   2. Calcula o fluxo normal em cada face (ferncodes_flowrate)
        %      usando a formula MPFA-D: q = Kde*(p_R - p_L - Ded*(z2-z1))
        function [flowrate, flowresult, flowratedif, faceaux] = ...
                calcularFlowrate(obj, p, env, parms)
            [pinterp, ~] = ferncodes_pressureinterpNLFVPP(p, env);
            [flowrate, flowresult, flowratedif, faceaux] = ...
                ferncodes_flowrate(p, pinterp, env);
        end

    end

    methods (Access = private)

        %% ── Pre-processamento de concentracao acoplada ───────────
        % Chamado apenas para simulacoes de transporte (casos 200-300, 350-400)
        % Calcula os parametros necessarios para a equacao de concentracao:
        %
        %   Kdec, Knc, Ktc, Dedc → transmissibilidades para o tensor de difusao
        %   wightc, sc            → pesos LPEW2 para a concentracao
        %   weightDMPc            → pesos DMP (discrete maximum principle)
        %   dparameter            → parametro de difusao numerica
        %   nflagnoc, nflagfacec  → flags de contorno para a concentracao
        function env = preprocessConcentracao(obj, env, N)
            % condicao inicial de concentracao
            [Con, lastimelevel, lastimeval] = applyinicialcond;

            % parametros auxiliares para o tensor de difusao molecular
            [~, Kdec, Knc, Ktc, Dedc, wightc, sc, weightDMPc, dparameter] = ...
                parametersauxiliary(env.config.dmap, N);

            % flags de contorno para a equacao de concentracao
            [nflagnoc, nflagfacec] = ferncodes_calflag_con(lastimeval);

            % armazena tudo em env.conpre.*
            env.conpre.Con          = Con;
            env.conpre.lastimelevel = lastimelevel;
            env.conpre.lastimeval   = lastimeval;
            env.conpre.Kdec         = Kdec;
            env.conpre.Knc          = Knc;
            env.conpre.Ktc          = Ktc;
            env.conpre.Dedc         = Dedc;
            env.conpre.wightc       = wightc;
            env.conpre.sc           = sc;
            env.conpre.weightDMPc   = weightDMPc;
            env.conpre.dparameter   = dparameter;
            env.conpre.nflagnoc     = nflagnoc;
            env.conpre.nflagfacec   = nflagfacec;
        end

    end
end