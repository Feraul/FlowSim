%--------------------------------------------------------------------------
%Function "applyinicialcond"
%--------------------------------------------------------------------------

function [Sw,lastimelevel,lastimeval] = applyinicialcond
%Define global parameters
global filepath benchkey

%Attribute "zero" for "lasttimelevel" and "lastimeval"
lastimelevel = 0;
lastimeval = 0;
%Verify if there exists a restart condition
command = [char(filepath) '\' 'Results_teste_ConReport.dat'];
restartkey = exist(command,'file');

%There exists a restart condition
if restartkey ~= 0 && strcmp(benchkey,'r')
    %Get a initial condition based in a last saturation field.
    [Sw,lastimelevel,lastimeval] = setrestartinicond;
    %There exists a restart condition, BUT the user does NOT set this option.
elseif restartkey ~= 0 && strcmp(benchkey,'r') == 0
    %It deletes the "restart.dat" file
    command = ['del ' char(filepath) '\' 'Results_teste_ConReport.dat'];
    %It calls system
    system(command);
    %Attribute INITIAL CONDITION
    [Sw,] = attribinitialcond;
    %There is NO a restart condition.
else
    %Attribute INITIAL CONDITION
    [Sw,] = attribinitialcond;
end%End of IF (restart condition)
end