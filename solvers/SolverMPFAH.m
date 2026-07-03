% solvers/SolverMPFAH.m
classdef SolverMPFAH < SolverBase
    properties
        Nome     = 'MPFA-H'
        MetodoID = 'mpfah'
    end
    methods
        function [premethod, preGravity] = preprocessar(obj, env, parmRichardEq)
            premethod = struct();
            kmap = env.config.kmap;

            [facelement]               = ferncodes_elementfacempfaH;
            [pointarmonic]             = ferncodes_harmonicopoint(kmap);
            [parameter, ~]             = ferncodes_coefficientmpfaH(facelement, pointarmonic, kmap);
            [weightDMP]                = ferncodes_weightnlfvDMP(kmap, env.geometry.elem);

            premethod.MPFAH.facelement   = facelement;
            premethod.MPFAH.pointarmonic = pointarmonic;
            premethod.MPFAH.parameter    = parameter;
            premethod.MPFAH.weightDMP    = weightDMP;

            preGravity = obj.calcGravidade(env);
        end

        function [h_new, flowrate] = resolver(obj, env, premethod, parmRichardEq, h, dt, time)
            [h_new, flowrate] = ferncodes_solverpressureMPFAH(...
                premethod, env, parmRichardEq, h, dt, time);
        end
    end
end