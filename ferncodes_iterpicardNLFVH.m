function [p,flowrate,flowresult]=ferncodes_iterpicardNLFVH(M_old,RHS_old,...
   parameter,p_old,nflagface,wells,mobility,weightDMP)
global pmethod
%% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);

%% inicializando dados para iteração Picard
step=0;
er=1;

while (tolpicard<er || tolpicard==er) && (step<nitpicard)
    % atualiza iterações
    step=step+1
    
    %--------------------------------------------------------------------------
    %Solve global algebric system
    
    % calculo das pressões
    p_new = solver(M_old,RHS_old);
    
    %% plotagem no visit
    S=ones(size(p_new,1),1);
    ferncodes_postprocessor(p_new,S,step)
    [pinterp_new]=ferncodes_pressureinterpHP(p_new,nflagface,parameter,weightDMP,mobility);
    %% Calculo da matriz global
    if strcmp(pmethod,'nlfvdmp')
        [M,I]=ferncodes_assemblematrixDMP(p_new,pinterp_new,0,parameter,weightDMP,mobility);
    else
        [M,I]=ferncodes_assemblematrixNLFVH(pinterp_new,parameter,mobility);
    end
    %--------------------------------------------------------------------------
    %Add a source therm to independent vector "mvector"
    
    %Often it may change the global matrix "M"
    [M_new,RHS_new] = addsource(sparse(M),I,wells);
    %% Calculo do residuo
    
    R = norm(M_new*p_new - RHS_new);
    
    if (R0 ~= 0.0)
        er = abs(R/R0)
    else
        er = 0.0; %exact
    end
    errorelativo(step)=er;
    
    % atualizar
    M_old=M_new;
    RHS_old=RHS_new;
    
end
%--------------------------------------------------------------------------
%Solve global algebric system

% calculo das pressões
p = p_new;

%Message to user:
fprintf('\n Iteration number, iterations = %d \n',step)
fprintf('\n Residual error, error = %d \n',er)
%Message to user:
disp('>> The Pressure field was calculated with success!');
[pinterp]=ferncodes_pressureinterpHP(p,nflagface,parameter,weightDMP,mobility);
%Get the flow rate (Diamond)
if strcmp(pmethod,'nlfvdmp')
    [flowrate,flowresult]=ferncodes_flowrateDMP(p, pinterp, parameter,nflagface,0,weightDMP,mobility);
else
    [flowrate,flowresult]=ferncodes_flowrateNLFVH(p, pinterp, parameter,mobility);
end
%Message to user:
disp('>> The Flow Rate field was calculated with success!');

end