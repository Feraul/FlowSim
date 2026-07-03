function [p,step,errorelativo,flowrate,flowresult]=ferncodes_iterpicardANLFVPP(M_old,RHS_old,...
    nitpicard,tolpicard,parameter,w,s,p_old,nflagno,wells,mobility)

%% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);

%% inicializando dados para iteração Picard
step=0;
er=1;
while (tolpicard<er || tolpicard==er) && (step<nitpicard)
    %% atualiza iterações
    step=step+1;
    %% Acelerador de Anderson
    if step>5
        [p_new]=ferncodes_andersonacc(p_old,parameter,w,s,...
            nflagno,mobility,wells,R0);
        % se o campo de ress�o � negativo ele coloca zero
        r=p_new(:)<0;
        x=find(r==1);
        if max(x)>0
            p_new(x)=0;
        end
     else
         p_new=M_old\RHS_old;  % invers�o sem pivotamento
     end
   
    % plotagem no visit
    S=ones(size(p_new,1),1);
    ferncodes_postprocessor(p_new,S,step)
    p_max=max(p_new);
    p_min=min(p_new);
    % Interpolação das pressões na arestas (faces)
    [pinterp_new]=ferncodes_pressureinterpNLFVPP(p_new,nflagno,w,s);
    
    % Calculo da matriz global
    [M,I]=ferncodes_assemblematrixNLFVPP(pinterp_new,parameter,mobility);
    %--------------------------------------------------------------------------
    %Add a source therm to independent vector "mvector"
    
    %Often it may change the global matrix "M"
    [M_new,RHS_new] = addsource(sparse(M),I,wells);
    %errorelativo(step)=er;
    
    R = norm(M_new*p_new - RHS_new);
    
    if (R0 ~= 0.0)
        er = abs(R/R0);
        %er = R/norm(RHS_new)
    else
        er = 0.0; %exact
    end
    errorelativo(step)=er;
    
    %% atualizar
    M_old=full(M_new);
    RHS_old=RHS_new;
    p_old=M_old\RHS_old;  % inversão sem pivotamento
    
end
%--------------------------------------------------------------------------
%Solve global algebric system

% calculo das press�es
p = solver(M_old,RHS_old);

%Message to user:
disp('>> The Pressure field was calculated with success!');
[pinterp]=ferncodes_pressureinterpNLFVPP(p,nflagno,w,s);
%Get the flow rate (Diamond)
[flowrate,flowresult]=ferncodes_flowrateNLFVPP(p, pinterp, parameter,mobility);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');

end