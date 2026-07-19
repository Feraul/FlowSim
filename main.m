%--------------------------------------------------------------------------
% MAIN — Simulador de fluxo em meios porosos
% Programer: Fernando Contreras
%--------------------------------------------------------------------------
% Fluxo de execucao:
%   1. Leitura do Start.dat e pre-processamento geometrico
%   2. Instanciacao dos tres objetos principais (benchmark, metodo, simulacao)
%   3. Pre-processamento fisico (h_init, parms, kmap, flags, premethod)
%   4. Chamada ao solver transiente ou estacionario
%--------------------------------------------------------------------------

% ── Path — rodar savepath() somente na primeira execucao ─────────
% Apos salvar o path com savepath(), o MATLAB encontra as pastas
% automaticamente em todas as execucoes seguintes.


%startup(); savepath();   % ← descomente apenas na primeira vez

clc;
clear classes   % limpa o cache de classes OOP do MATLAB
                % obrigatorio sempre que um arquivo .m de classe for editado
tic

%% ── 1. Pre-processamento geral ───────────────────────────────────
% Le o arquivo Start.dat e popula:
%   env.config.*    → numcase, pmethod, phasekey, perm, bcflag, totaltime...
%   env.geometry.*  → coord, elem, bedge, inedge, centelem, elemarea...
env = preprocessormod(1);

%% ── 2. Instancia os tres objetos a partir do Start.dat ──────────

% benchmark: encapsula a fisica especifica do numcase
%   → permeabilidade, BC, theta, capacidade hidrica, loop temporal
env.benchmark = createBenchmark(env.config.numcase);

% metodo: encapsula o metodo numerico (TPFA, MPFA-D, MPFA-H, NL-TPFA)
%   → pre-processamento geometrico, montagem de matriz, solver linear
env.metodo = createMetodo(env.config.pmethod);

% sim: encapsula o tipo de simulacao (phasekey 0..6)
%   → pre-processamento fisico do tipo (Richards, groundwater, etc.)
%   → definicao de fontes e pocos
sim = createSimulacao(env.config.phasekey);

%% ── 3. Pipeline de pre-processamento ────────────────────────────

% 3a. Pre-processamento fisico do tipo de simulacao
%     Preenche parms com: h_init, h_old, dt, theta_s, theta_r, alpha...
%     Cada benchmark implementa seu proprio preprocessar()
parms        = env.benchmark.initParms();    % struct vazio com campos padrao
[env, parms] = sim.preprocessar(env, parms); % preRE, prehydraulic, preRichards...

% 3b. Permeabilidade
%     Calcula kmap (tensor de permeabilidade por elemento) via benchmark
%     Popula: env.config.kmap, env.config.auxkmap, env.utils.*
[env, parms] = PLUG_kfunction(env, parms, 0);

% 3c. Flags de contorno
%     Monta nflag (por vertice) e nflagface (por face) via benchmark
%     Popula: env.config.nflag, env.config.nflagface
%     DEVE rodar antes de preprocessmethod (que usa nflagface)
[env] = ferncodes_calflag(env, parms, 0);

% 3d. Pre-processamento do metodo numerico
%     Calcula os parametros geometrico-fisicos especificos do metodo:
%       TPFA  → Hesq, Kde, Kn, flowrateZ
%       MPFA-D → Hesq, Kde, Kn, Kt, Ded, pesos LPEW2, s
%     Popula: env.premethod.TPFA.* ou env.premethod.MPFAD.*
%     Popula: env.preGravity.*
[env, parms] = preprocessmethod(env, parms);

%% ── 4. Solver ────────────────────────────────────────────────────
% Escolhe e executa o solver de acordo com phasekey:
%   case 1 → One-phase flow (estacionario)
%   case 4 → Groundwater / hydraulic head
%   case 5 → Contaminant + hydraulic head
%   case 6 → Richards (transiente nao-linear)
%
% source_wells: pocos injetores/produtores definidos pelo benchmark
setmethod(sim, sim.definirFontes(env, parms), 'i', 8, env, parms);

toc