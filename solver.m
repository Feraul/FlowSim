
%--------------------------------------------------------------------------
%Subject: numerical routine to solve flux flow in porous media
%Type of file: FUNCTION
%--------------------------------------------------------------------------
%Goals: %This function solves the global algebric system.
%"globalM" is a global matrix

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [pressure] = solver(globalmatrix,vector)
%Solver the algebric system using matlab's routines (black box!)
global numcase

if numcase==347
    % preconditioner
    %globalmatrix=globalmatrix+1*eye(size(globalmatrix,1));

    [L,U] = ilu(sparse(globalmatrix),struct('type','ilutp','droptol',1e-6));

    %[p_old,fl1,rr1,it1,rv1]=bicgstab(M_old,RHS_old,1e-10,1000,L,U);
    [pressure,fl1,rr1,it1,rv1]=gmres(globalmatrix,vector,10,1e-4,1000,L,U);
else
    pressure = globalmatrix\vector;
end
end