function [flowrate,flowresult,flowratec,flowresultc]=ferncodes_flowrateNLFVPP(p, ...
    pinterp, parameter,viscosity,Con,nflagc,wightc,sc,dparameter,cinterp,gravrate,weightDMPc)
global inedge coord bedge bcflag centelem numcase bcflagc keygravity

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);


%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(centelem,1),1);
flowratec = zeros(bedgesize + inedgesize,1);
flowresultc = zeros(size(centelem,1),1);
valuemin=1e-16;
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
        %flowrate(ifacont,1)= normcont*bcflag(r,2);% testes feitos em todos os
        %problemas monofásico
        flowrate(ifacont,1)= -normcont*bcflag(r,2);% problema de buckley leverett Bastian
    else
        % contribuicao do termo gravitacional
        if strcmp(keygravity,'y')
            m=gravrate(ifacont);
        else
            m=0;
        end
        flowrate(ifacont,1)=visonface*normcont*(parameter(1,1,ifacont)*(p(lef)-pinterp(parameter(1,3,ifacont)))+...
            parameter(1,2,ifacont)*(p(lef)-pinterp(parameter(1,4,ifacont))))-visonface*m;
    end
    %Attribute the flow rate to "flowresult"
    %On the left:
    
    
    flowresult(lef) = flowresult(lef) + flowrate(ifacont);
    %% ================================================================
    if 200<numcase && numcase<300
        if bedge(ifacont,7)>200
            % condicao de contorno de Neumann
            x=bcflagc(:,1)==bedge(ifacont,7);
            r=find(x==1);
            %flowrate(ifacont,1)= normcont*bcflag(r,2);% testes feitos em todos os
            %problemas monofásico
            flowratec(ifacont,1)= -normcont*bcflagc(r,2);% problema de buckley leverett Bastian
        else
            % condicao de contorno de Dirichlet
            flowratec(ifacont,1)=normcont*(dparameter(1,1,ifacont)*(Con(lef)-cinterp(dparameter(1,3,ifacont)))+...
                dparameter(1,2,ifacont)*(Con(lef)-cinterp(dparameter(1,4,ifacont))));
        end
        %Attribute the flow rate to "flowresult"
        %On the left:
        flowresultc(lef) = flowresultc(lef) + flowratec(ifacont);
    end
end

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
    norma=norm(vd1);
    ifactual=iface+size(bedge,1);
    
    %% calculo do a Eq. 2.7 (resp. eq. 16) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    % esquerda
    alef=(parameter(1,1,ifactual)*pinterp(parameter(1,3,ifactual))+...
        parameter(1,2,ifactual)*pinterp(parameter(1,4,ifactual)));
    % direita
    
    arel= (parameter(2,1,ifactual)*pinterp(parameter(2,3,ifactual))+...
        parameter(2,2,ifactual)*pinterp(parameter(2,4,ifactual)));
    %% calculo dos "mu", Eq. 2.8 (resp. eq. 18) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    if alef==0 && arel==0
        mulef= 0.5;
        murel=1-mulef;
    else
        mulef=abs(arel)/(abs(alef)+abs(arel));
        murel=1-mulef;
    end
    %% calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    ALL=norma*mulef*(parameter(1,1,ifactual)+parameter(1,2,ifactual));
    ALR=norma*murel*(parameter(2,1,ifactual)+parameter(2,2,ifactual));
    
    % contribuicao do termo gravitacional
    if strcmp(keygravity,'y')
        
        m=gravrate(iface+size(bedge,1));
    else
        m=0;
        
    end
    flowrate(iface+size(bedge,1),1)=visonface*(ALL*p(lef)-ALR*p(rel))-visonface*m;
    
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);
    
    %% ===================================================================
    if 200<numcase && numcase<300
        % calculo do fluxo para o campo de concentracoes
        % calculo do a Eq. 2.7 (resp. eq. 16) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
        % esquerda
        alefc=(dparameter(1,1,ifactual)*cinterp(dparameter(1,3,ifactual))+...
            dparameter(1,2,ifactual)*cinterp(dparameter(1,4,ifactual)));
        % direita
        
        arelc= (dparameter(2,1,ifactual)*cinterp(dparameter(2,3,ifactual))+...
            dparameter(2,2,ifactual)*cinterp(dparameter(2,4,ifactual)));
        % calculo dos "mu", Eq. 2.8 (resp. eq. 18) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
        if alefc==0 && arelc==0
            mulefc= 0.5;
            murelc=1-mulefc;
        else
            mulefc=abs(arelc)/(abs(alefc)+abs(arelc));
            murelc=1-mulefc;
        end
        
        % calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
        ALLc=norma*mulefc*(dparameter(1,1,ifactual)+dparameter(1,2,ifactual));
        ALRc=norma*murelc*(dparameter(2,1,ifactual)+dparameter(2,2,ifactual));
        
        flowratec(iface+size(bedge,1),1)=(ALLc*Con(lef)-ALRc*Con(rel))+norma*(murelc*arelc-mulefc*alefc);
     
        %Attribute the flow rate to "flowresult"
        %On the left:
        flowresultc(lef) = flowresultc(lef) + flowratec(bedgesize + iface);
        %On the right:
        flowresultc(rel) = flowresultc(rel) - flowratec(bedgesize + iface);
    end
end

end