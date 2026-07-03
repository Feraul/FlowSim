function [flowrate,flowresult]=ferncodes_flowratelfvHP(parameter,...
    weightDMP,pinterp,p,viscosity)
global inedge coord bedge bcflag centelem numcase
%---------------------------------------------------------
% viscosity=viscosity when contaminant transport
% viscosity=mobility when two-phase flow
%---------------------------------------------------------

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(centelem,1),1);

for ifacont=1:bedgesize
    % when 200<numcase<300, is groundwater transport model
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
        % when numcase<200, is two-phase model
    elseif 30<numcase && numcase<200
        % equacao de saturacao "viscosity=mobility"
        visonface=sum(viscosity(ifacont,:));
    else
        % otherwise, the model is groundwater flow model
        visonface=1;
    end
    lef=bedge(ifacont,3);
    
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    
    if bedge(ifacont,5)>200
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        %flowrate(ifacont,1)=normcont*bcflag(r,2);% testes feitos em todos os
        %problemas monofásico
        flowrate(ifacont,1)=-normcont*bcflag(r,2);% problema de buckley leverett Bastian
    else
        facelef1=parameter(1,3,ifacont);
        facelef2=parameter(1,4,ifacont);
        % average hydralic head
        % unconfined aquifer
        flowrate(ifacont,1)= visonface*normcont*((parameter(1,1,ifacont)+parameter(1,2,ifacont))*p(lef)-...
            parameter(1,1,ifacont)*pinterp(facelef1)-parameter(1,2,ifacont)*pinterp(facelef2));
    end
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(ifacont);
    
end
% Montagem da matriz global

for iface=1:inedgesize
    % when 200<numcase<300, is groundwater transport model
    if 200<numcase && numcase<300
        % concentration equation
        if numcase == 246 || numcase == 245 || numcase==247 ||...
                numcase==248 || numcase==249 || numcase==251
            % vicosity on the boundary edge
            visonface = viscosity(bedgesize + iface,:);
            %It is a Two-phase flow
        else
            visonface = 1;
        end  %End of IF
    elseif 30<numcase && numcase<200
        % saturation equation "viscosity=mobility"
        visonface=sum(viscosity(bedgesize + iface,:));
    else
        % single-phase or groundwater flow model
        visonface=1;
    end
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    
    % orientation vector
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma=norm(vd1);
    ifactual=iface+size(bedge,1);
    % Calculo das contribuições do elemento a esquerda
    mulef=weightDMP(ifactual-size(bedge,1),1);
    murel=weightDMP(ifactual-size(bedge,1),2);
    
    % os nós que conforman os pontos de interpolação no elemento a esquerda
    auxfacelef1=parameter(1,3,ifactual);
    auxfacelef2=parameter(1,4,ifactual);
    % os nós que conforman os pontos de interpolação no elemento a direita
    auxfacerel1=parameter(2,3,ifactual);
    auxfacerel2=parameter(2,4,ifactual);
    % calculo dos fluxo parcial a esquerda
    fluxesq=norma*((parameter(1,1,ifactual)+parameter(1,2,ifactual))*p(lef)-...
        parameter(1,1,ifactual)*pinterp(auxfacelef1)-parameter(1,2,ifactual)*pinterp(auxfacelef2));
    % calculo dos fluxo parcial a direita
    fluxdireit=norma*((parameter(2,1,ifactual)+parameter(2,2,ifactual))*p(rel)-...
        parameter(2,1,ifactual)*pinterp(auxfacerel1)-parameter(2,2,ifactual)*pinterp(auxfacerel2));
    % calculo do fluxo unico na face
    flowrate(iface+size(bedge,1),1)=visonface*(murel*fluxesq-mulef*fluxdireit);
    
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);
    
end
end