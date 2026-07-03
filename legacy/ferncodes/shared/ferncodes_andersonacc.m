%% Artur: Para uma boa convergencia tem mexer em "atol" e "rtol" a eficiencia do acelerador muito depende eles
% eu nos meus testes utilizei 1.e^-6, eu acho que vocÃª tera que fazer algum estudo parametrico para achar um 
% um adequado.



function [x,iter,res_hist] = ferncodes_andersonacc(x,...
    parameter,w,s,nflagno,mobility,wells,R0)
global pmethod
% This performs fixed-point iteration with or without Anderson
% acceleration for a given fixed-point map g and initial
% approximate solution x.
%%
% https://users.wpi.edu/~walker/Papers/anderson_accn_algs_imps.pdf
%%
% Required inputs:
%        g = fixed-point map (function handle); form gval = g(x).
%        x = initial approximate solution (column vector).
%
% Optional inputs:
%       mMax = maximum number of stored residuals (non-negative integer).
% NOTE: mMax = 0 => no acceleration.
%       itmax = maximum allowable number of iterations.
%       atol = absolute error tolerance.
%       rtol = relative error tolerance.
%       droptol = tolerance for dropping stored residual vectors to improve
% conditioning: If droptol > 0, drop residuals if the
% condition number exceeds droptol; if droptol <= 0,
% do not drop residuals.
%       beta = damping factor: If beta > 0 (and beta ~= 1), then the step is
% damped by beta; otherwise, the step is not damped.
% NOTE: beta can be a function handle; form beta(iter), where iter is
% the iteration number and 0 < beta(iter) <= 1.
%       AAstart = acceleration delay factor: If AAstart > 0, start acceleration
%       when iter = AAstart.
%
% Output:
%       x = final approximate solution.
%       iter = final iteration number.
%       res_hist = residual history matrix (iteration numbers and residual norms).
%
% Homer Walker (walker@wpi.edu), 10/14/2011.
% Set the method parameters. 
%######################################################################
if nargin < 2 
    error('AndAcc requires at least two arguments.'); 
end
%if nargin < 3 
    mMax = min(10,size(x,1)); 
%end
%if nargin < 4 
    itmax = 1000; 
%end
%if nargin < 5 
    atol = 1.e-6; 
%end
%if nargin < 6 
    rtol = 1.e-6; 
%end
%if nargin < 7 
    droptol = 1.e10; 
%end
%if nargin < 8
    beta = 1; 
%end
%if nargin < 9
    AAstart = 20; 
%end
% Initialize the storage arrays.
res_hist = []; % Storage of residual history.
DG = []; % Storage of g-value differences.
% Initialize printing.
% if mMax == 0
%     fprintf('\n No acceleration.');
% elseif mMax > 0
%     fprintf('\n Anderson acceleration, mMax = %d \n',mMax);
% else
%     error('AndAcc.m: mMax must be non-negative.');
% end
% fprintf('\n iter res_norm \n');
% Initialize the number of stored residuals.
mAA = 0;
% Top of the iteration loop.
for iter = 0:itmax
    % Apply g and compute the current residual norm.
    % interpolação nos nós ou faces
    if strcmp(pmethod,'nlfvpp')
        [pinterp_new]=ferncodes_pressureinterpNLFVPP(x,nflagno,w,s);
        [M,I]=ferncodes_assemblematrixNLFVPP(pinterp_new,parameter,mobility);
        %Often it may change the global matrix "M"
        [M_new,RHS_new] = addsource(sparse(M),I,wells);
    elseif strcmp(pmethod,'nlfvpp')
        
        
    elseif strcmp(pmethod,'nlfvdmp')
        % interpolação nos nós ou faces
        [pinterp_new]=ferncodes_pressureinterpHP(x,nflagface,parameter,weightDMP,mobility);
        % Montagem da matriz global
        
        [M,I]=ferncodes_assemblematrixDMP(x,pinterp_new,0,parameter,weightDMP,mobility);
        %--------------------------------------------------------------------------
        %Add a source therm to independent vector "mvector"
        
        %Often it may change the global matrix "M"
        [M_new,RHS_new] = addsource(sparse(M),I,wells);
    else
        [pinterp_new]=ferncodes_pressureinterpHP(x,nflagface,parameter,weightDMP,mobility);
        %% Calculo da matriz global
        [M,I]=ferncodes_assemblematrixNLFVH(pinterp_new,parameter,mobility);
        %--------------------------------------------------------------------------
        %Add a source therm to independent vector "mvector"
        
        %Often it may change the global matrix "M"
        [M_new,RHS_new] = addsource(sparse(M),I,wells);
    end
    RR = norm(M_new*x - RHS_new);
    
    if (R0 ~= 0.0)
        er = abs(RR/R0)
    else
        er = 0.0; %exact
    end
    
     %% InterpolaÃ§Ã£o das pressÃµes na arestas (faces)
%     [pinterp_new]=pressureinterp(x,nflagface,nflagno,w,s,auxflag,metodoP,...
%         parameter,weightDMP);
%                                       
%     %% Calculo da matriz global
%     [M_new,RHS_new]=globalmatrix(x,pinterp_new,gamma,nflagface,nflagno...
%         ,parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,...
%         mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
    if er<1e-8
        break
    else
        %[L,U] = ilu(M_new,struct('type','ilutp','droptol',1e-6));
       
        %[gval,fl1,rr1,it1,rv1]=bicgstab(M_new,RHS_new,1e-9,1000,L,U);
        gval=M_new\RHS_new;
    end
    %gval=M_new\RHS_new;
    fval = gval - x;
    res_norm = norm(fval);
    fprintf('%d %e \n', iter, res_norm);
    res_hist = [res_hist;[iter,res_norm]];
    % Set the residual tolerance on the initial iteration.
    if iter == 0,
        tol = max(atol,rtol*res_norm);
    end
    % Test for stopping.
    if res_norm <= tol || er<1e-8
        %fprintf('Terminate with residual norm = %e \n\n', res_norm);
        break;
    end
    if mMax == 0 || iter < AAstart,
        % Without acceleration, update x <- g(x) to obtain the next
        % approximate solution.
        x = gval;
    else
        % Apply Anderson acceleration.
        % Update the df vector and the DG array.
        if iter > AAstart,
            df = fval-f_old;
            if mAA < mMax,
                DG = [DG gval-g_old];
            else
                DG = [DG(:,2:mAA) gval-g_old];
            end
            mAA = mAA + 1;
        end
        f_old = fval;
        g_old = gval;
        if mAA == 0
            % If mAA == 0, update x <- g(x) to obtain the next approximate solution.
            x = gval;
        else
            % If mAA > 0, solve the least-squares problem and update the
            % solution.
            if mAA == 1
                % If mAA == 1, form the initial QR decomposition.
                R(1,1) = norm(df);
                Q = R(1,1)\df;
            else
                % If mAA > 1, update the QR decomposition.
                if mAA > mMax
                    % If the column dimension of Q is mMax, delete the first column and
                    % update the decomposition.
                    [Q,R] = qrdelete(Q,R,1);
                    mAA = mAA - 1;
                    % The following treats the qrdelete quirk described below.
                    if size(R,1) ~= size(R,2),
                        Q = Q(:,1:mAA-1); R = R(1:mAA-1,:);
                    end
                    % Explanation: If Q is not square, then qrdelete(Q,R,1) reduces the
                    % column dimension of Q by 1 and the column and row
                    % dimensions of R by 1. But if Q *is* square, then the
                    % column dimension of Q is not reduced and only the column
                    % dimension of R is reduced by one. This is to allow for
                    % MATLABâ€™s default "thick" QR decomposition, which always
                    % produces a square Q.
                end
                % Now update the QR decomposition to incorporate the new
                % column.
                for j = 1:mAA - 1
                    R(j,mAA) = Q(:,j)'*df;
                    df = df - R(j,mAA)*Q(:,j);
                end
                R(mAA,mAA) = norm(df);
                Q = [Q,R(mAA,mAA)\df];
            end
            if droptol > 0
                % Drop residuals to improve conditioning if necessary.
                condDF = cond(R);
                while condDF > droptol && mAA > 1
                   % fprintf(' cond(D) = %e, reducing mAA to %d \n', condDF, mAA-1);
                    [Q,R] = qrdelete(Q,R,1);
                    DG = DG(:,2:mAA);
                    mAA = mAA - 1;
                    % The following treats the qrdelete quirk described above.
                    if size(R,1) ~= size(R,2),
                        Q = Q(:,1:mAA); R = R(1:mAA,:);
                    end
                    condDF = cond(R);
                end
            end
            % Solve the least-squares problem.
            gamma = R\(Q'*fval);
            % Update the approximate solution.
            x = gval - DG*gamma;
            % Apply damping if beta is a function handle or if beta > 0
            % (and beta ~= 1).
            if isa(beta,'function_handle'),
                x = x - (1-beta(iter))*(fval - Q*R*gamma);
            else
                if beta > 0 && beta ~= 1,
                    x = x - (1-beta)*(fval - Q*R*gamma);
                end
            end
        end
    end
    
    

end
 %Bottom of the iteration loop.
%  if res_norm > tol && iter == itmax,
%      fprintf('\n Terminate after itmax = %d iterations. \n', itmax);
%      fprintf(' Residual norm = %e \n\n', res_norm);
%  end