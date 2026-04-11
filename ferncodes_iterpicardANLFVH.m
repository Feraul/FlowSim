function [p,flowrate,flowresult,flowratedif,flowresultc]=...
    ferncodes_iterpicardANLFVH(M_old,RHS_old,parameter,p_old,nflagface,...
    wells,viscosity,weightDMP,w,s,nflagno,contnorm,weightDMPc,Con,...
    nflagfacec,dparameter,wightc,sc,SS,dt,h,MM,gravrate)


%% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);


%% Acelerador de Anderson
%p_old=M_old\RHS_old;
 [L,U] = ilu(M_old,struct('type','ilutp','droptol',1e-6));

%[p_old,fl1,rr1,it1,rv1]=bicgstab(M_old,RHS_old,1e-10,1000,L,U);
[p_old,fl1,rr1,it1,rv1]=gmres(M_old,RHS_old,10,1e-9,1000,L,U);
%[p_new]=ferncodes_andersonacc(p_old,parameter,w,s,nflagno,mobility,wells,...
%    weightDMP,nflagface);
[p,erro,iter]=ferncodes_andersonacc2(p_old,1e-6,parameter,w,s,...
    nflagface,weightDMP,wells,viscosity,R0,contnorm,Con,wightc,sc,...
    weightDMPc,nflagfacec,dparameter,SS,dt,h,MM,gravrate);

% se o campo de ressão é negativo ele coloca zero
%Message to user:
fprintf('\n Iteration number, iterations = %d \n',iter)
fprintf('\n Residual error, error = %d \n',erro)
disp('>> The Pressure field was calculated with success!');
[pinterp,cinterp]=ferncodes_pressureinterpHP(p,nflagface,parameter,...
    weightDMP,weightDMPc,Con,nflagfacec,dparameter);
[flowrate,flowresult,flowratedif,flowresultc]=ferncodes_flowrateNLFVH(p,...
    pinterp, parameter,viscosity,Con,cinterp,dparameter);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');
end
