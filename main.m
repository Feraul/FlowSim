
%Subject: numerical code used to simulate fluid flow in porous media. That
%routine calls several others which defines how the equation will be solved
%Type of file: MAIN
%Programer: PhD. Fernando Contreras,
%--------------------------------------------------------------------------
%Goals: Do the manegement of simulator. This is a MAIN program.

%--------------------------------------------------------------------------
% In this numerical routine the flow may be simulated with one or two phase
% or contaminants or groundwater (hydraulic head). The functions below are
% organized in order give flexibility to software resourcers.
% For example: the saturation and pressure fields are calculated by IMPES,
% but it may be calculated also by a fully implicit scheme. This change is
% done in the present rountine just call each function.

%|--------------------------------------|
%| read the instructions very carefully |
%| Warning: See preMPFA.m line 53 and 54|
%|--------------------------------------|

%Clear the screem
clc;
%Clear all the memory of the matlab
clear all;
%Define the format of out data
format short;
%It begins the time counter and "profile".
%% pre-processador
tic
env = preprocessormod(1);
parameters = initParams();
parmconcentra = [];
parmgroundwater = [];
parmRichardEq = [];
%% define os pocos e fonte ou sumidoro para cada caso

%% Inicializações seguras

bedgesize  = size(env.geometry.bedge,1);
inedgesize = size(env.geometry.inedge,1);

kmapaux = env.config.perm(1,1);
env.config.perm=[1 kmapaux 0  0 kmapaux];


% ============================================================
% CASOS 200–300 → Concentração
% ============================================================
if 200 < env.config.numcase && env.config.numcase < 300
    source_wells= defineWells(env,parmRichardEq);
    [env, parmconcentra] = preconcentration(env, source_wells.wells);

    if ismember(env.config.numcase, [247 249 250])
        load('Perm_Var0p1.mat')
        perm = flipud(perm);
        auxperm2 = perm(1:125,1:125);
        kmap = auxperm2';
        kmap = kmap(:);
    elseif env.config.numcase == 251
        kmap = kmap;
    end
% ============================================================
% CASOS 300–350 → Hidráulica
% ============================================================
elseif 300 < env.config.numcase && env.config.numcase < 350
    source_wells= defineWells(env,parmRichardEq);
    [parmgroundwater, source_wells] = prehydraulic;

    if ismember(env.config.numcase, [341 380.1 341.1])
        Nmod = 100;
        varK = 0.1;
    end
    if 350 < env.config.numcase && env.config.numcase < 400
        [env, parmconcentra] = preconcentration(env, source_wells);
    end
% ============================================================
% CASOS 400–500 → Richards
% ============================================================
elseif 400 < env.config.numcase && env.config.numcase < 500

    [parmRichardEq,env] = preRE(env);
    source_wells= defineWells(env,parmRichardEq);

end

% ============================================================
% Define the norm of permeability or conductivity hidraulic 
% tensor ("normk")
[env,parmRichardEq] = calcnormk(env, parmRichardEq,0);
%=============================================================
% flag boundary condition
[env] = ferncodes_calflag(env,parmRichardEq,0);
% ============================================================
% CHAMADA FINAL
% ============================================================
setmethod(source_wells, 'i', 8, env, parmconcentra, parmgroundwater,...
          parmRichardEq);
toc