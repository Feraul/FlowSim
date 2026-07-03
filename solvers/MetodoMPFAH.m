%--------------------------------------------------------------------------
% MetodoMPFAH — Metodo de Volumes Finitos com Pontos Harmonicos (MPFA-H)
%
% Implementa o metodo MPFA-H (Multi-Point Flux Approximation — Harmonic points)
% para aproximacao do fluxo difusivo. Alternativa ao MPFA-D que usa pontos
% harmonicos em vez de patches diamante para interpolacao nas faces.
%
% Origem: refatorado a partir de SolverMPFAH.m (orphan que herdava de
% SolverBase — classe base ausente do repo). Agora herda de MetodoBase para
% unificar a hierarquia de metodos numericos. Ver PR-A2 / ADR-004.
%
% Referencias:
%   MPFA-H — Hybrid MPFA with harmonic points (Contreras et al.)
%
% Nota de escopo PR-A2:
%   Esta classe INSTANTIA (satisfaz o contrato MetodoBase) e delega o
%   pre-processamento real. A integracao completa com o loop hydraulic_RE
%   (montarSistema, resolver, calcularFlowrate cheios) e feita em PR-C4
%   quando reescrevemos ferncodes_assemblematrixMPFAH em forma triplet.
%--------------------------------------------------------------------------
classdef MetodoMPFAH < MetodoBase

    properties
        Nome     = 'MPFA-H'
        MetodoID = 'mpfah'
    end

    methods

        %% ── 1. Pre-processamento geometrico e fisico ─────────────
        % Chamado UMA vez em preprocessmethod(), apos ferncodes_calflag.
        % Popula env.premethod.MPFAH.* com:
        %   facelement   → mapeamento elemento-face (ferncodes_elementfacempfaH)
        %   pointarmonic → pontos harmonicos (ferncodes_harmonicopoint)
        %   parameter    → coeficientes de transmissibilidade MPFA-H
        %   weightDMP    → pesos DMP nas faces internas
        function [env, parms] = preprocessar(obj, env, parms)
            kmap = env.config.kmap;

            [facelement]   = ferncodes_elementfacempfaH;
            [pointarmonic] = ferncodes_harmonicopoint(kmap);
            [parameter, ~] = ferncodes_coefficientmpfaH(facelement, pointarmonic, kmap);
            [weightDMP]    = ferncodes_weightnlfvDMP(kmap, env.geometry.elem);

            env.premethod.MPFAH.facelement   = facelement;
            env.premethod.MPFAH.pointarmonic = pointarmonic;
            env.premethod.MPFAH.parameter    = parameter;
            env.premethod.MPFAH.weightDMP    = weightDMP;

            env.preGravity = obj.calcGravidade(env);
        end

        %% ── 2. Atualiza premethod no loop temporal ───────────────
        % MPFA-H tipicamente e usado para problemas lineares (steady-state).
        % Se kmap mudar (ex: Richards), recalcula parameter + weightDMP.
        function [env] = atualizarPremethod(obj, env, parms) %#ok<INUSD>
            kmap = env.config.kmap;
            fe   = env.premethod.MPFAH.facelement;
            ph   = env.premethod.MPFAH.pointarmonic;

            [parameter, ~] = ferncodes_coefficientmpfaH(fe, ph, kmap);
            [weightDMP]    = ferncodes_weightnlfvDMP(kmap, env.geometry.elem);

            env.premethod.MPFAH.parameter = parameter;
            env.premethod.MPFAH.weightDMP = weightDMP;
        end

        %% ── 3. Monta a matriz global A e vetor b ─────────────────
        % Delega ao ferncodes_assemblematrixMPFAH legacy (820 linhas).
        % PR-C4 substitui por versao vetorizada (triplet form) sob flag
        % env.config.useVectAssembly.
        function [M, I] = montarSistema(obj, env, parms, dt) %#ok<INUSD>
            % Argumentos herdados do wrapper legacy ferncodes_solverpressureMPFAH:
            parameter = env.premethod.MPFAH.parameter;
            weightDMP = env.premethod.MPFAH.weightDMP;
            nflagface = env.config.nflagface;
            SS        = getOr(env.config, 'SS',        0);
            h         = getOr(parms,      'h',         []);
            MM        = getOr(env.config, 'MM',        0);
            gravrate  = getOr(env.preGravity, 'gravrate', 0);
            viscosity = getOr(env.config, 'viscosity', 1);

            [M, I, ~] = ferncodes_assemblematrixMPFAH( ...
                parameter, nflagface, weightDMP, SS, dt, h, MM, gravrate, viscosity);
        end

        %% ── 4. Resolve o sistema linear A*h = b ──────────────────
        % MPFA-H para problemas LINEARES: solve direto, sem Picard.
        % Casos nao-lineares (se surgirem) devem usar MetodoMPFAD ou MetodoNLFVPP.
        function [p, flowrate, flowresult, flowratedif, faceaux, parms, env] = ...
                resolver(obj, M, I, parms, env, tempo, dt, source_wells) %#ok<INUSD>
            % wells + source ja aplicados em montarSistema (contrato MetodoBase).
            % Se ainda precisar adicionar aqui, descomente:
            % [M, I] = addsource(sparse(M), I, source_wells);
            p = solver(M, I);

            % Calcula flowrate via calcularFlowrate (contrato MetodoBase)
            [flowrate, flowresult, flowratedif, faceaux] = obj.calcularFlowrate(p, env, parms);
        end

        %% ── 5. Calcula flowrate apos convergencia ────────────────
        function [flowrate, flowresult, flowratedif, faceaux] = ...
                calcularFlowrate(obj, p, env, parms) %#ok<INUSD>
            parameter = env.premethod.MPFAH.parameter;
            weightDMP = env.premethod.MPFAH.weightDMP;
            nflagface = env.config.nflagface;
            viscosity = getOr(env.config, 'viscosity', 1);

            [pinterp] = ferncodes_pressureinterpHP(p, nflagface, parameter, weightDMP, 0, 0, 0, 0);
            [flowrate, flowresult] = ferncodes_flowratelfvHP( ...
                parameter, weightDMP, pinterp, p, viscosity);

            % Contrato MetodoBase pede 4 saidas; MPFA-H nao separa difusivo/faceaux
            flowratedif = flowrate;
            faceaux     = [];
        end

    end
end

% ── helper local: getOr(struct, field, default) ────────────────────────────
function v = getOr(s, f, d)
    if isstruct(s) && isfield(s, f), v = s.(f); else, v = d; end
end
