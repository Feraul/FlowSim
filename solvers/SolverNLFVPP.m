% solvers/SolverNLFVPP.m
classdef SolverNLFVPP < SolverBase
    properties
        Nome     = 'NL-TPFA / NLFVPP'
        MetodoID = 'nlfvpp'
    end
    methods
        function [premethod, preGravity] = preprocessar(obj, env, parmRichardEq)
            premethod = struct();

            [V, N, F] = ferncodes_elementface(env);
            [premethod, weight, s] = ferncodes_Pre_LPEW_2_vect(premethod, parmRichardEq, env);
            [parameter, contnorm]  = ferncodes_coefficient(env.config.kmap, env.geometry.elem);

            premethod.NLTPFA.V        = V;
            premethod.NLTPFA.N        = N;
            premethod.NLTPFA.weight   = weight;
            premethod.NLTPFA.s        = s;
            premethod.NLTPFA.p_old    = 1e1*ones(size(env.geometry.elem,1),1);
            premethod.NLTPFA.parameter = parameter;
            premethod.NLTPFA.contnorm  = contnorm;

            preGravity = obj.calcGravidade(env);
        end

        function [h_new, flowrate] = resolver(obj, env, premethod, parmRichardEq, h, dt, time)
            [h_new, flowrate] = ferncodes_solverpressureNLFVPP(...
                premethod, env, parmRichardEq, h, dt, time);
        end
    end
end