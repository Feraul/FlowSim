% simulacoes/SimGroundwater.m   (tipo 4)
classdef SimGroundwater < SimulacaoBase

    properties
        Nome   = 'Groundwater (Hydraulic Head)'
        TipoID = 4
    end

    methods
        function [env, parms] = preprocessar(obj, env)
            source_wells = obj.definirFontes(env, []);
            [parms, source_wells] = prehydraulic(env, source_wells);
        end

        function K = calcularconductividade(obj, env)
            K = env.config.perm;   % sua logica aqui
        end

        function [nflag, nflagface] = condicoesContorno(obj, env)
            nflag     = [];   % sua logica aqui
            nflagface = [];
        end

        function wells = definirFontes(obj, env, parms)
            wells = defineWells(env, parms);
        end
    end

end