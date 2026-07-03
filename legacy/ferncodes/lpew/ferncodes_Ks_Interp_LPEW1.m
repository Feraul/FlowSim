%It is called by the function "ferncodes_Pre_LPEW_1.m"

function [Kt1,Kt2,Kn1,Kn2] = ferncodes_Ks_Interp_LPEW1(O,T,Qo,kmap,No)
%Retorna os K(n ou t) necessários para a obtenção dos weights. kmap é a
%matriz de permeabilidade; Ex: Kt1->linhaN=Kt1(cellN);

global bedge inedge esurn2 esurn1 phasekey visc elem numcase

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Verifica quantos são os elementos em torno do nó "ni".%
N_element_No=esurn2(No+1)-esurn2(No);

%Prealocação das matrizes.%

Kt1=zeros(N_element_No,2); %As duas colunas correspondem a i=1 e i=2.
Kt2=zeros(N_element_No,2);
Kn1=zeros(N_element_No,2);
Kn2=zeros(N_element_No,2);
K=zeros(3);
K1=zeros(3);

R=[0 1 0; -1 0 0;0 0 0];
%--------------------------------------------------------------------------

%permeability construction

for k=1:N_element_No
    
    j=esurn1(esurn2(No)+k);
    for icont=1:2
        if (size(T,1)==size(O,1))&&(k==N_element_No)&&(icont==2)
                        
            K(1,1)=kmap(elem(j,5),2);
            K(1,2)=kmap(elem(j,5),3);
            K(2,1)=kmap(elem(j,5),4);
            K(2,2)=kmap(elem(j,5),5);
            
            Kn1(k,icont)=((R*(T(1,:)-Qo)')'*K*(R*(T(1,:)-Qo)'))/(norm(T(1,:)-Qo)^2);
            Kt1(k,icont)=((R*(T(1,:)-Qo)')'*K*((T(1,:)-Qo)'))/(norm(T(1,:)-Qo)^2);
        else
            K(1,1)=kmap(elem(j,5),2);
            K(1,2)=kmap(elem(j,5),3);
            K(2,1)=kmap(elem(j,5),4);
            K(2,2)=kmap(elem(j,5),5);
                     
            Kn1(k,icont)=((R*(T(k+icont-1,:)-Qo)')'*K*(R*(T(k+icont-1,:)-Qo)'))/(norm(T(k+icont-1,:)-Qo)^2);
            Kt1(k,icont)=((R*(T(k+icont-1,:)-Qo)')'*K*((T(k+icont-1,:)-Qo)'))/(norm(T(k+icont-1,:)-Qo)^2);
        end
    end
    
        
    %------------------------- Tensores ----------------------------------%
    K1(1,1)=kmap(elem(j,5),2);
    K1(1,2)=kmap(elem(j,5),3);
    K1(2,1)=kmap(elem(j,5),4);
    K1(2,2)=kmap(elem(j,5),5);
    
    % calculo dos outros K(n ou t) no paper é denotado com um "~" na parte
    % inferior de K
    for icont=0:1
        if (size(T,1)==size(O,1))&&(k==N_element_No)&&(icont==1)
            
            Kn2(k,icont+1)=((R*(O(k,:)-T(1,:))')'*K1*(R*(O(k,:)-T(1,:))'))/(norm(O(k,:)-T(1,:))^2);
            Kt2(k,icont+1)=((R*(O(k,:)-T(1,:))')'*K1*((O(k,:)-T(1,:))'))/(norm(O(k,:)-T(1,:))^2);
            
        else
            
            Kn2(k,icont+1)=((R*(O(k,:)-T(k+icont,:))')'*K1*(R*(O(k,:)-T(k+icont,:))'))/(norm(O(k,:)-T(k+icont,:))^2);
            Kt2(k,icont+1)=((R*(O(k,:)-T(k+icont,:))')'*K1*((O(k,:)-T(k+icont,:))'))/(norm(O(k,:)-T(k+icont,:))^2);
        end
    end
    
end

end

