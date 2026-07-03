% startup.m — rode uma vez ou coloque no inicio do main
function startup()
    base = fileparts(mfilename('fullpath'));   % pasta raiz do projeto

    addpath(fullfile(base, 'base'));
    addpath(fullfile(base, 'simulacoes'));
    addpath(fullfile(base, 'benchmarks'));
    addpath(fullfile(base, 'solvers'));
    addpath(fullfile(base, 'factories'));
    %addpath(fullfile(base, 'core'));

    fprintf('Paths configurados com sucesso.\n');
end