%--------------------------------------------------------------------------
% MetodoTPFA — Metodo de Volumes Finitos de Dois Pontos (TPFA)
%
% Implementa o metodo TPFA (Two-Point Flux Approximation) para
% aproximacao do fluxo difusivo. E o metodo mais simples e eficiente,
% porem limitado a malhas ortogonais ou tensores de permeabilidade
% alinhados com a malha (sem termos de anisotropia cruzada).
%
% Diferenca em relacao ao MPFA-D:
%   TPFA  → fluxo calculado com DOIS pontos (centroide L e centroide R)
%            sem interpolacao nos vertices, sem Ded, sem pesos LPEW2
%   MPFA-D → fluxo calculado com MULTIPLOS pontos (patch diamante)
%             captura anisotropia completa do tensor K
%
% Responsabilidades:
%   preprocessar()       → calcula Hesq, Kde, Kn (sem Ded, Kt, pesos)
%   atualizarPremethod() → recalcula Kde quando kmap muda (Richards)
%   montarSistema()      → monta a matriz global A e vetor b
%   resolver()           → escolhe o iterador (Picard, AA, L-scheme) e resolve
%   calcularFlowrate()   → calcula fluxo nas faces via formula TPFA
%--------------------------------------------------------------------------
classdef MetodoTPFA < MetodoBase

    properties
        Nome     = 'TPFA'    % nome para exibicao no console
        MetodoID = 'tpfa'    % identificador string (lido do Start.dat)
    end

    methods

        %% ── 1. Pre-processamento geometrico e fisico ─────────────
        % Chamado UMA vez antes do loop temporal (em preprocessmethod).
        % Calcula e armazena em env.premethod.TPFA.*:
        %
        %   Hesq        → altura do triangulo formado pelo centroide e a aresta
        %                 (distancia ortogonal do centroide a face)
        %   Kde         → transmissibilidade normal da aresta interna
        %                 Kde = -|e| * Kn_L * Kn_R / (Kn_L*H_R + Kn_R*H_L)
        %   Kn          → permeabilidade normal nas faces de contorno
        %   flowrateZ   → fluxo gravitacional por face (Richards)
        %   flowresultZ → fluxo gravitacional acumulado por elemento
        %
        % TPFA NAO calcula Ded, Kt nem pesos LPEW2
        % (esses termos existem apenas no MPFA-D para capturar anisotropia)
        function [env, parms] = preprocessar(obj, env, parms)
            [env]          = ferncodes_Kde_Ded_Kt_Kn_TPFA(env, parms);
            env.preGravity = obj.calcGravidade(env);
        end

        %% ── 2. Atualiza premethod no loop temporal ───────────────
        % Chamado em hydraulic_RE a CADA passo de tempo, apos PLUG_kfunction
        % atualizar kmap com o novo h.
        %
        % Para TPFA apenas Kde precisa ser recalculado (depende de Kn que
        % depende de kmap). NAO ha pesos LPEW2 para recalcular.
        function [env] = atualizarPremethod(obj, env, parms)
            [env] = ferncodes_Kde_Ded_Kt_Kn_TPFA(env, parms);
        end

        %% ── 3. Monta a matriz global do sistema linear ───────────
        % Seleciona a rotina de montagem adequada ao tipo de problema:
        %
        %   numcase == 331 ou 400-500 (Richards):
        %     → ferncodes_globalmatrix_TPFA: monta A com termo temporal
        %       (capacidade hidrica via soil_properties)
        %
        %   outros (estacionario ou groundwater):
        %     → ferncodes_globalmatrix: monta A sem termo temporal
        function [M, I] = montarSistema(obj, env, parms, dt)
            if env.config.numcase == 331 || ...
                    (400 < env.config.numcase && env.config.numcase < 500)
                % Richards / fluxo nao-saturado — inclui termo dtheta/dt
                [M, I] = ferncodes_globalmatrix_TPFA(env, parms);
            else
                % estacionario ou groundwater — sem termo temporal
                [M, I, ~] = ferncodes_globalmatrix(env, env.premethod.TPFA, parms, dt);
            end
        end

        %% ── 4. Resolve o sistema linear A*h = b ──────────────────
        % Seleciona o iterador de acordo com env.config.acel:
        %
        %   'FPI'     → Picard classico (Fixed Point Iteration)
        %   'AA'      → Picard com Aceleracao de Anderson
        %   'LSCHEME' → L-scheme (estabilizacao por regularizacao diagonal)
        %
        % Identico ao MetodoMPFAD — o iterador e independente do metodo
        % de discretizacao espacial (a matriz A ja foi montada antes)
        function [p, flowrate, flowresult, flowratedif, faceaux, parms, env] = ...
                resolver(obj, M, I, parms, env, tempo, dt, source_wells)

            switch upper(env.config.acel)

                case 'FPI'
                    % Picard classico
                    [p, flowrate, flowresult, flowratedif, faceaux, parms, env] = ...
                        ferncodes_iterpicard(M, I, parms, env, tempo, dt, source_wells);

                case 'AA'
                    % Aceleracao de Anderson — reduz numero de iteracoes
                    [p, flowrate, flowresult, flowratedif, faceaux, parms] = ...
                        ferncodes_iterpicardANLFVPP2(M, I, env.premethod.TPFA, ...
                        parms, env, tempo, dt, source_wells);

                case 'LSCHEME'
                    % L-scheme — adiciona L*I a diagonal para garantir estabilidade
                    [p, flowrate, flowresult, flowratedif, faceaux, parms, ...
                        env.premethod.TPFA] = ...
                        L_scheme(M, I, env.premethod.TPFA, parms, env, ...
                        tempo, dt, source_wells);
            end
        end

        %% ── 5. Calcula o flowrate (fluxo normal nas faces) ───────
        % Formula TPFA — usa apenas os dois centroides adjacentes:
        %   q_face = Kde * (h_R - h_L)
        %
        % NAO interpola pressao nos vertices (diferente do MPFA-D)
        % NAO usa pesos LPEW2
        % Adequado para malhas ortogonais ou K alinhado com a malha
        function [flowrate, flowresult, flowratedif, faceaux] = ...
                calcularFlowrate(obj, p, env, parms)
            [flowrate, flowresult, flowratedif, faceaux] = ...
                ferncodes_flowrateTPFA(p, env);
        end

    end
end