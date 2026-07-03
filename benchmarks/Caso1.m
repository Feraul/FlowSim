classdef Caso1 < BenchmarkBase

    properties
        Nome = 'Caso1'
    end

    methods
        function K = calcularconductividade(obj)
            K = eye(2);
        end

        function [nflag, nflagface] = condicoesContorno(obj)
            nflag    = [];
            nflagface = [];
        end
    end

end