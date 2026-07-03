function [M,I]=ferncodes_assemblematrixNLFVPP(pinterp,parameter,viscosity,...
    contnorm,SS,dt,h,MM,gravrate,nflag)
global inedge coord bedge bcflag elem numcase keygravity ...
    methodhydro dens modflowcompared elemarea;
%-----------------------inicio da rOtina ----------------------------------%
%Constrói a matriz global.

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
coeficiente=dt^-1*MM*SS;
%Initialize "M" (global matrix) and "I" (known vector)
M = sparse(size(elem,1),size(elem,1)); %Prealocação de M.
I = zeros(size(elem,1),1);
valuemin=1e-16;
jj=1;
m=0;
for ifacont=1:bedgesize
    
    if 200<numcase && numcase<300
        % equacao de concentracao
        if numcase == 246 || numcase == 245 || numcase==247 || ...
                numcase==248 || numcase==249 ||numcase==251
            % vicosity on the boundary edge
            visonface = viscosity(ifacont,:);
            %It is a Two-phase flow
        else
            visonface = 1;
        end  %End of IF
    elseif numcase<200
        % equacao de saturacao "viscosity=mobility"
        visonface=sum(viscosity(ifacont,:));
    else
        visonface=1;
    end
    
    lef=bedge(ifacont,3);
    
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    
    if bedge(ifacont,5)>200
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        %I(lef)=I(lef)- normcont*bcflag(r,2); % testes feitos em todos os
        %problemas monofásico
        I(lef)=I(lef)+ normcont*bcflag(r,2); % problema de buckley leverett Bastian
    else
        if strcmp(modflowcompared,'y')
            elembedge(jj,1)=bedge(ifacont,3);
            elembedge(jj,2)=nflag(bedge(ifacont,2),2);
            jj=jj+1;
        end
        %% calculo da contribuição do contorno, veja Eq. 2.17 (resp. eq. 24) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
        % contribuicao do termo gravitacional
        if strcmp(keygravity,'y')
            if numcase<200
                % escoamento bifasico oleo-agua
                averagedensity=(viscosity(ifacont,:)*dens')/visonface;
                m=averagedensity*gravrate(ifacont);
            else
                % concentracao soluto-solvente
                m=dens(1,1)*gravrate(ifacont)/visonface;
            end
        end
        
        alef= visonface*normcont*(parameter(1,1,ifacont)*pinterp(parameter(1,3,ifacont))+...
            parameter(1,2,ifacont)*pinterp(parameter(1,4,ifacont)));
        
        Alef=visonface*normcont*(parameter(1,1,ifacont)+parameter(1,2,ifacont));
        
        %% implementação da matriz global no contorno
        M(lef,lef)=M(lef,lef)+ Alef;
        I(lef,1)=I(lef,1)+alef+ visonface*m;
    end
end

%% Montagem da matriz global

for iface=1:inedgesize
    if 200<numcase && numcase<300
        % equacao de concentracao
        if numcase == 246 || numcase == 245 || numcase==247 ||...
                numcase==248 || numcase==249 || numcase==251
            % vicosity on the boundary edge
            visonface = viscosity(bedgesize + iface,:);
            %It is a Two-phase flow
        else
            visonface = 1;
        end  %End of IF
    elseif numcase<200
        % equacao de saturacao "viscosity=mobility"
        visonface=sum(viscosity(bedgesize + iface,:));
    else
        visonface=1;
    end
    
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma= sqrt(vd1(1,1)^2+vd1(1,2)^2);
    ifactual=iface+size(bedge,1);
    
    % calculo do a Eq. 2.7 (resp. eq. 16) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    % esquerda
   
    alef=parameter(1,1,ifactual)*pinterp(parameter(1,3,ifactual))+...
        parameter(1,2,ifactual)*pinterp(parameter(1,4,ifactual));
    
    % direita
   
    arel= parameter(2,1,ifactual)*pinterp(parameter(2,3,ifactual))+...
        parameter(2,2,ifactual)*pinterp(parameter(2,4,ifactual));
    
    mulef=(abs(arel)+1e-16)/(abs(alef)+abs(arel)+2*1e-16);
    
    murel=(abs(alef)+1e-16)/(abs(alef)+abs(arel)+2*1e-16);
    %      mulef=(abs(arel)+contnorm^2)/(abs(alef)+abs(arel)+2*contnorm^2);
    
    %      murel=(abs(alef)+contnorm^2)/(abs(alef)+abs(arel)+2*contnorm^2);
    
    % calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    
    %% calculo da contribuição, Eq. 2.10 (resp. eq. 19) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    Bsigma=murel*arel-mulef*alef;
    
    Bmas=(abs(Bsigma)+Bsigma)/2;
    Bmenos=(abs(Bsigma)-Bsigma)/2;
    
    ALL=norma*mulef*(parameter(1,1,ifactual)+parameter(1,2,ifactual));
    
    ARR=norma*murel*(parameter(2,1,ifactual)+parameter(2,2,ifactual));
    % implementação da matriz global
    % contribuição da transmisibilidade no elemento esquerda
    M(lef,lef)=M(lef,lef)+ visonface*ALL;
    M(lef,rel)=M(lef,rel)- visonface*ARR;
    % contribuição da transmisibilidade no elemento direita
    M(rel,rel)=M(rel,rel)+ visonface*ARR;
    M(rel,lef)=M(rel,lef)- visonface*ALL;
    
    
   % I(lef, 1)=I(lef,1)- (norma*Bmas*valuemin/(p(lef)+valuemin)-norma*Bmenos*valuemin/(p(rel)+valuemin));
   % I(rel, 1)=I(rel,1)+ (norma*Bmas*valuemin/(p(lef)+valuemin)-norma*Bmenos*valuemin/(p(rel)+valuemin));
    
    
    % contribuicao do termo gravitacional
    if strcmp(keygravity,'y')
         if numcase<200
            % escoamento bifasico oleo-agua
            averagedensity=(viscosity(bedgesize + iface,:)*dens')/visonface;
            m=averagedensity*gravrate(bedgesize + iface,1);
        else
            % concentracao soluto-solvente
            m=dens(1,1)*gravrate(bedgesize + iface,1)/visonface;
        end
        
        I(lef)=I(lef)+visonface*m;
        I(rel)=I(rel)-visonface*m;
    end
end
%==========================================================================
% utilizase somente quando o teste vai ser comparado com resultados do
% modflow
if strcmp(modflowcompared,'y')
    for iw = 1:size(elembedge,1)
        M(elembedge(iw,1),:)=0*M(elembedge(iw,1),:);
        M(elembedge(iw,1),elembedge(iw,1))=1;
        I(elembedge(iw,1))=elembedge(iw,2);
    end
end
% para calcular a carga hidraulica
% calcula um problema transiente
[M,I]=ferncodes_implicitandcranknicolson(M,I,SS,dt,MM,h);

end