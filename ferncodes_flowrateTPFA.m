function [flowrate, flowresult,flowratedif]=ferncodes_flowrateTPFA(p,Kde,Kn,Hesq,nflag,viscosity,Kdec,Knc,nflagc,Con)

global coord bedge inedge bcflag bcflagc centelem numcase

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "bedgeamount"
bedgeamount = 1:bedgesize;

%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(centelem,1),1);
flowratedif=0;
for ifacont=1:size(bedge,1);
    lef=bedge(ifacont,3);
    if numcase == 246 || numcase == 245 || numcase==247 || numcase==248 || numcase==249
        % vicosity on the boundary edge
        visonface = viscosity(ifacont,:);
        %It is a Two-phase flow
    else
        visonface = 1;
    end  %End of IF
    nor=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    % calculo das constantes nas faces internas
    A=-Kn(ifacont)/(Hesq(ifacont)*nor);
    if bedge(ifacont,5)<200 % se os nós esteverem na fronteira de DIRICHLET
        c1=nflag(bedge(ifacont,1),2);
        c2=nflag(bedge(ifacont,2),2);
        
        flowrate(ifacont)=visonface*A*(nor^2)*(c1-p(lef));
    else
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        flowrate(ifacont)= nor*bcflag(r,2);
    end
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(ifacont);
    if 200<numcase && numcase<300
        %% ====================================================================
        if bedge(ifacont,7)<200 % se os nós esteverem na fronteira de DIRICHLET
            c1aux=nflagc(bedge(ifacont,1),2);
            Ac=-Knc(ifacont)/(Hesq(ifacont)*nor);
            flowratedif(ifacont)=Ac*(nor^2)*(c1aux-Con(lef));
        else
            x=bcflagc(:,1)==bedge(ifacont,7);
            r=find(x==1);
            flowratedif(ifacont)= nor*bcflagc(r,2);
        end
    end
end

for iface=1:size(inedge,1)
    lef=inedge(iface,3); %indice do elemento a direita da aresta i
    rel=inedge(iface,4); %indice do elemento a esquerda da aresta i
    if numcase == 246 || numcase == 245 || numcase==247 || numcase==248 || numcase==249
        % vicosity on the boundary edge
        visonface = viscosity(bedgesize + iface,:);
        %It is a Two-phase flow
    else
        visonface = 1;
    end  %End of IF
    %-------------------- calculo das vazões e velocidades ---------------%
    
    flowrate(iface+size(bedge,1))=visonface*Kde(iface)*(p(rel)-p(lef));
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);
    if 200<numcase && numcase<300
        %% ====================================================================
        flowratedif(iface+size(bedge,1))=Kdec(iface)*(Con(rel)-Con(lef));
    end
end

end