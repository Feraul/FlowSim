%It is called by "ferncodes_solverpressure.m"

function [preMPFAD,weight,s] = ferncodes_Pre_LPEW_2(preMPFAD,parmRichardEq,env)

auxnumcase=env.config.numcase;
if auxnumcase>400
    kmap=parmRichardEq.auxperm;
else
    kmap=env.config.perm;
end
elem=env.geometry.elem;
N=preMPFAD.N;
nflag=env.config.nflag;

% Retorna todos os parâmetros necessários às expressões dos fluxos.
apw = ones(size(env.geometry.coord,1),1);
r = zeros(size(env.geometry.coord,1),2);
s=0;
for y = 1:size(env.geometry.coord,1),
    No = y;
    % calculos dos vetores O, P, T, Q
    [O,P,T,Qo] = OPT_Interp_LPEW(No,env);
    % calculo dos angulos
    [ve2,ve1,theta2,theta1] = angulos_Interp_LPEW2(O,P,T,Qo,No,env);
    % calculo dos netas
    [neta] = netas_Interp_LPEW(O,P,T,Qo,No,env);
    
    % calculo dos Ks
    [Kt1,Kt2,Kn1,Kn2] = ferncodes_Ks_Interp_LPEW2(O,T,Qo,No,zeros(size(elem,1),1),env,kmap);
    
    % calculo dos lamdas
    [lambda,r] = Lamdas_Weights_LPEW2(Kt1,Kt2,Kn1,Kn2,theta1,theta2,ve1,...
        ve2,neta,P,O,Qo,No,T,r);
    % calculo dos pesos
    for k = 0:size(O,1) - 1,
        weight(apw(No) + k) = lambda(k + 1)/sum(lambda); %Os pesos fazem sentido%%%%%%%%%%%
    end
    
    apw(No + 1) = apw(No) + size(O,1);
    if env.config.numcase==341 
        % interpolaçao das pressões nos contornos de Neumann
        vetor = env.geometry.nsurn1(env.geometry.nsurn2(No) + 1:env.geometry.nsurn2(No + 1));
        comp1 = N(No,1);
        comp2 = N(No,length(vetor));
        % verifica se o vertices pertence ao contorno de Neumann
        if 200<nflag(No,1) && nflag(No,1)<300
            % avalia se a face comp1 esta no contorno de Neumann
            if env.geometry.bedge(comp1,5)>200
                a = env.config.bcflag(:,1) == env.geometry.bedge(comp1,5);
                s1 = find(a == 1);
                % o "r" ja esta acompanhado pela norma
                %------------------------------------------------------
                aa=0.5*(env.geometry.coord(env.geometry.bedge(comp1,1),:) + env.geometry.coord(env.geometry.bedge(comp1,2),:));
                auxkmap = ferncodes_K(aa(1,1),aa(1,2));
                
                %------------------------------------------------------
                %auxkmap=kmap(bedge(comp1,3),2);
                aux1= r(No,1)*auxkmap(1)*nflagface(s1,2);
            end
            % avalia se a face comp1 esta no contorno de Neumann
            if env.geometry.bedge(comp2,5)>200
                b = env.config.bcflag(:,1) == env.geometry.bedge(comp2,5);
                s2 = find(b == 1);
                %
                %-------------------------------------------------------
                
                aaa=0.5*(env.geometry.coord(env.geometry.bedge(comp2,1),:) + env.geometry.coord(env.geometry.bedge(comp2,2),:));
                auxkmap = ferncodes_K(aaa(1,1),aaa(1,2));
                
                %------------------------------------------------------
                % auxkmap=kmap(bedge(comp2,3),2);
                
                % o "r" ja esta acompanhado pela norma
                aux2= r(No,2)*auxkmap(1)*nflagface(s2,2);
            end
            s(No,1) = -(1/sum(lambda))*(aux1+ aux2);
        end
        
    else
        % interpolaçao das pressões nos contornos de Neumann
        vetor = env.geometry.nsurn1(env.geometry.nsurn2(No) + 1:env.geometry.nsurn2(No + 1));
        comp1 = N(No,1);
        comp2 = N(No,length(vetor));
        MM=env.geometry.bedge(:,1)==No;
        MMM= find(MM == 1);
        if comp1<= size(env.geometry.bedge,1) && comp2 <=size(env.geometry.bedge,1) && 200<env.geometry.bedge(MMM,4)
            a = env.config.bcflag(:,1) == env.geometry.bedge(comp1,5);
            s1 = find(a == 1);
            b = env.config.bcflag(:,1) == env.geometry.bedge(comp2,5);
            s2 = find(b == 1);
            
            s(No,1) = -(1/sum(lambda))*(r(No,1)*env.config.bcflag(s1,2) + ...
                r(No,2)*env.config.bcflag(s2,2));
        end  %End of IF
    end
end

  preMPFAD.weight=weight;
  preMPFAD.s=s;
end
%End of FOR




