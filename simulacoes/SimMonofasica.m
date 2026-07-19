classdef SimMonofasica < SimulacaoBase

    properties
        Nome   = 'Simulacao Monofasica Estacionaria'
        TipoID = 1
    end

    methods
        function [env, parms] = configurarPermeabilidade(obj, env, parms, time)
            [env, parms] = env.benchmark.configurarPermeabilidade(env, parms, time);
        end

        function bcattrib = configurarContorno(obj, vertices, flagptr, time, env, pR)
            bcattrib = env.benchmark.configurarContorno(vertices, flagptr, time, env, pR);
        end

        function [nflag, nflagface] = configurarFlags(obj, env, pR, time)
            [nflag, nflagface] = env.benchmark.configurarFlags(env, pR, time);
        end

        function [parms, extras] = inicializar(obj, env, parms, time)
            [parms, extras] = env.benchmark.inicializar(env, parms, time);
        end

        function [parms, extras] = atualizarEstado(obj, env, parms, extras, h, flowrate, time, count)
            [parms, extras] = env.benchmark.atualizarEstado(env, parms, extras, h, flowrate, time, count);
        end

        function finalizar(obj, env, extras, varargin)
            env.benchmark.finalizar(env, extras, varargin{:});
        end

        function escreverResultados(obj, env, varargin)
            env.benchmark.escreverResultados(env, varargin{:});
        end

        function [env, parms] = preprocessar(obj, env, parms)
            [parms, env] = env.benchmark.preprocessar(env, parms);
        end

        function wells = definirFontes(obj, env, pR)
            wells = env.benchmark.definirFontes(env, pR);
        end
    end

end