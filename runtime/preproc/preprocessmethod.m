
%--------------------------------------------------------------------------
%Subject: numerical routine to solve flux flow in porous media
%Type of file: FUNCTION
%Programer: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%This FUNCTION calculate the

%--------------------------------------------------------------------------

function [env,parms] = preprocessmethod(env,parms)
disp('>> Preprocessing Pressure or Hydraulic head Equation...');

% ── Delega ao solver ─────────────────────────────────────────
[env,parms] = env.metodo.preprocessar(env, parms);

disp('>> "preprocessmethod" finalizado com sucesso!');

end
