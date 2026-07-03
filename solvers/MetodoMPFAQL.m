%--------------------------------------------------------------------------
% MetodoMPFAQL — Multi-Point Flux Approximation — Quasi-Linear scheme
%
% Variante do MPFA que usa pesos baseados em quadrilateros locais.
% Referencia: Contreras et al. 2019.
%
% Origem: classe NOVA criada em PR-A2 (nao existia arquivo previo apesar
% de createMetodo.m dispatcher expor `case 'mpfaql': metodo = MetodoMPFAQL()`).
% Ver ADR-004.
%
% Nota de escopo PR-A2:
%   Instantiavel + delega para wrappers legacy. Vetorizacao em PR-C5.
%--------------------------------------------------------------------------
classdef MetodoMPFAQL < MetodoBase

    properties
        Nome     = 'MPFA-QL'
        MetodoID = 'mpfaql'
    end

    methods

        %% ── 1. Pre-processamento ─────────────────────────────────
        % MPFA-QL usa LPEW2 (weight/s) + parametros MPFA-QL specificos + weightDMP
        function [env, parms] = preprocessar(obj, env, parms)
            [V, N, F] = ferncodes_elementface(env);
            [env, weight, s] = ferncodes_Pre_LPEW_2_vect(env, parms);
            [parameter, ~]   = ferncodes_coefficient(env.config.kmap, env.geometry.elem);
            [weightDMP]      = ferncodes_weightnlfvDMP(env.config.kmap, env.geometry.elem);

            env.premethod.MPFAQL.V         = V;
            env.premethod.MPFAQL.N         = N;
            env.premethod.MPFAQL.F         = F;
            env.premethod.MPFAQL.weight    = weight;
            env.premethod.MPFAQL.s         = s;
            env.premethod.MPFAQL.parameter = parameter;
            env.premethod.MPFAQL.weightDMP = weightDMP;

            env.preGravity = obj.calcGravidade(env);
        end

        %% ── 2. Atualiza premethod no loop temporal ───────────────
        function [env] = atualizarPremethod(obj, env, parms) %#ok<INUSD>
            [env, weight, s] = ferncodes_Pre_LPEW_2_vect(env, parms);
            [parameter, ~]   = ferncodes_coefficient(env.config.kmap, env.geometry.elem);
            [weightDMP]      = ferncodes_weightnlfvDMP(env.config.kmap, env.geometry.elem);
            env.premethod.MPFAQL.weight    = weight;
            env.premethod.MPFAQL.s         = s;
            env.premethod.MPFAQL.parameter = parameter;
            env.premethod.MPFAQL.weightDMP = weightDMP;
        end

        %% ── 3. Monta a matriz global ─────────────────────────────
        function [M, I] = montarSistema(obj, env, parms, dt) %#ok<INUSD>
            pre = env.premethod.MPFAQL;
            nflag = env.config.nflag;
            mobility = getOr(env.config, 'mobility', ones(size(env.geometry.inedge, 1), 1));

            [M, I] = ferncodes_assemblematrixMPFAQL( ...
                pre.parameter, pre.weight, pre.s, nflag, pre.weightDMP, mobility);
        end

        %% ── 4. Resolve o sistema linear ──────────────────────────
        % MPFA-QL para problemas LINEARES — solve direto, sem Picard.
        function [p, flowrate, flowresult, flowratedif, faceaux, parms, env] = ...
                resolver(obj, M, I, parms, env, tempo, dt, source_wells) %#ok<INUSD>
            p = solver(M, I);
            [flowrate, flowresult, flowratedif, faceaux] = obj.calcularFlowrate(p, env, parms);
        end

        %% ── 5. Calcula flowrate ──────────────────────────────────
        function [flowrate, flowresult, flowratedif, faceaux] = ...
                calcularFlowrate(obj, p, env, parms) %#ok<INUSD>
            pre = env.premethod.MPFAQL;
            nflag = env.config.nflag;
            mobility = getOr(env.config, 'mobility', ones(size(env.geometry.inedge, 1), 1));

            [pinterp] = ferncodes_pressureinterpMPFAQL(p, nflag, pre.weight, pre.s);
            [flowrate, flowresult] = ferncodes_flowratelfvMPFAQL( ...
                pre.parameter, pre.weightDMP, mobility, pinterp, p);

            flowratedif = flowrate;
            faceaux     = [];
        end

    end
end

function v = getOr(s, f, d)
    if isstruct(s) && isfield(s, f), v = s.(f); else, v = d; end
end
