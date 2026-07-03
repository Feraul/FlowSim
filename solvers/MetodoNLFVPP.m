%--------------------------------------------------------------------------
% MetodoNLFVPP — Non-Linear Finite Volume with Positive Preserving property
%
% Implementa o esquema NLFV-PP (Contreras et al. 2021) — metodo nao-linear
% que garante o principio do maximo discreto (sem oscilacoes espurias).
% Mais caro que MPFA-D pois requer iteracao interna para os pesos nao-lineares.
% Recomendado quando MPFA-D produz valores negativos em problemas altamente
% anisotropicos.
%
% Origem: refatorado a partir de SolverNLFVPP.m (orphan que herdava de
% SolverBase — classe base ausente). Agora herda de MetodoBase.
% Ver PR-A2 / ADR-004.
%
% Nota de escopo PR-A2:
%   Instantiavel + delega para wrappers legacy. Vetorizacao completa em PR-C5.
%--------------------------------------------------------------------------
classdef MetodoNLFVPP < MetodoBase

    properties
        Nome     = 'NL-TPFA / NLFVPP'
        MetodoID = 'nlfvpp'
    end

    methods

        %% ── 1. Pre-processamento ─────────────────────────────────
        % NLFVPP precisa de LPEW2 (weights + s) + parametros de conectividade
        % Popula env.premethod.NLTPFA.* com V, N, F, weight, s, parameter, contnorm, p_old
        function [env, parms] = preprocessar(obj, env, parms)
            [V, N, F] = ferncodes_elementface(env);
            [env, weight, s] = ferncodes_Pre_LPEW_2_vect(env, parms);
            [parameter, contnorm] = ferncodes_coefficient(env.config.kmap, env.geometry.elem);

            env.premethod.NLTPFA.V         = V;
            env.premethod.NLTPFA.N         = N;
            env.premethod.NLTPFA.F         = F;
            env.premethod.NLTPFA.weight    = weight;
            env.premethod.NLTPFA.s         = s;
            env.premethod.NLTPFA.p_old     = 1e1 * ones(size(env.geometry.elem, 1), 1);
            env.premethod.NLTPFA.parameter = parameter;
            env.premethod.NLTPFA.contnorm  = contnorm;

            env.preGravity = obj.calcGravidade(env);
        end

        %% ── 2. Atualiza premethod no loop temporal ───────────────
        function [env] = atualizarPremethod(obj, env, parms) %#ok<INUSD>
            [env, weight, s] = ferncodes_Pre_LPEW_2_vect(env, parms);
            [parameter, contnorm] = ferncodes_coefficient(env.config.kmap, env.geometry.elem);
            env.premethod.NLTPFA.weight    = weight;
            env.premethod.NLTPFA.s         = s;
            env.premethod.NLTPFA.parameter = parameter;
            env.premethod.NLTPFA.contnorm  = contnorm;
        end

        %% ── 3. Monta a matriz global ─────────────────────────────
        % Precisa de pinterp que depende de p_old — assim, montarSistema
        % assume p_old em env.premethod.NLTPFA.p_old (atualizado por resolver).
        function [M, I] = montarSistema(obj, env, parms, dt) %#ok<INUSD>
            pre = env.premethod.NLTPFA;

            % interpolacao inicial das pressoes (usa iterate corrente)
            % Nota: NLFVPP requer atualizacao iterativa dentro do resolver
            [pinterp, ~] = ferncodes_pressureinterpNLFVPP(pre.p_old, env);

            viscosity = getOr(env.config, 'viscosity', 1);
            SS        = getOr(env.config, 'SS',        0);
            h         = getOr(parms,      'h',         []);
            MM        = getOr(env.config, 'MM',        0);
            gravrate  = getOr(env.preGravity, 'gravrate', 0);
            nflag     = env.config.nflag;

            [M, I] = ferncodes_assemblematrixNLFVPP( ...
                pinterp, pre.parameter, viscosity, pre.contnorm, ...
                SS, dt, h, MM, gravrate, nflag);
        end

        %% ── 4. Resolve com iterador Picard/AA ────────────────────
        function [p, flowrate, flowresult, flowratedif, faceaux, parms, env] = ...
                resolver(obj, M, I, parms, env, tempo, dt, source_wells) %#ok<INUSD>
            % Delega ao iterador Picard/AA — mesma logica de MetodoMPFAD,
            % adaptada para os campos NLTPFA em env.premethod.
            pre = env.premethod.NLTPFA;
            faceaux = [];

            switch upper(env.config.acel)
                case 'FPI'
                    [p, flowrate, flowresult, flowratedif, faceaux, parms, env] = ...
                        ferncodes_iterpicard(M, I, parms, env, tempo, dt, source_wells);
                case 'AA'
                    [p, flowrate, flowresult, flowratedif, faceaux, parms] = ...
                        ferncodes_iterpicardANLFVPP2(M, I, pre, parms, env, tempo, dt, source_wells);
                otherwise
                    % fallback: solve direto (sem Picard)
                    p = solver(M, I);
                    [flowrate, flowresult, flowratedif, faceaux] = ...
                        obj.calcularFlowrate(p, env, parms);
            end
        end

        %% ── 5. Calcula flowrate ──────────────────────────────────
        function [flowrate, flowresult, flowratedif, faceaux] = ...
                calcularFlowrate(obj, p, env, parms) %#ok<INUSD>
            [pinterp, ~] = ferncodes_pressureinterpNLFVPP(p, env);
            [flowrate, flowresult, flowratedif, faceaux] = ferncodes_flowrate(p, pinterp, env);
        end

    end
end

function v = getOr(s, f, d)
    if isstruct(s) && isfield(s, f), v = s.(f); else, v = d; end
end
