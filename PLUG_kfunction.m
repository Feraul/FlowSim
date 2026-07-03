
%--------------------------------------------------------------------------
%Subject: numerical routine to solve flux flow in porous media
%Type of file: FUNCTION
%Programer: Fernando R. L. Contreras
%--------------------------------------------------------------------------
%Goals: this FUNCTION gets the conductivity hidraulic.
%--------------------------------------------------------------------------

%Fill the matrix of permeability as a function of element to left and right
%of half edge evaluated. This function receives "kmap" and a feature of
%element which wants to know the permeability ("kfeature").
function [env,parms] = PLUG_kfunction(env, parms, time)

% ── Utilitários compartilhados (ficam aqui, disponiveis a todos) ──
    nelem              = size(env.geometry.centelem, 1);
    idx                = (1:nelem)';
    env.utils.nelem    = nelem;
    env.utils.idx      = idx;
    env.utils.iso      = @(coef) [idx, coef, zeros(nelem,1), zeros(nelem,1), coef];

    % ── Delega ao benchmark ───────────────────────────────────────────
    if isfield(env, 'benchmark') && ~isempty(env.benchmark)
        [env, parms] = env.benchmark.configurarPermeabilidade(env, parms, time);
        return
    end
    error('PLUG_kfunction: env.benchmark nao definido. Chame createBenchmark() no main.');
end