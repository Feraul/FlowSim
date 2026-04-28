function [flowrate,flowresult,flowratedif,faceaux] = ...
ferncodes_flowrate(p,pinterp,preMPFAD,parmRichardEq,env)

coord=env.geometry.coord;
bedge=env.geometry.bedge;
inedge=env.geometry.inedge;
centelem=env.geometry.centelem;
bcflag=env.config.bcflag;
numcase=env.config.numcase;
viscosity=env.config.visc;
flowrateZ=preMPFAD.flowrateZ;
kmap=parmRichardEq;

%--------------------------------------------------------------------------
Kde=preMPFAD.Kde;
Ded=preMPFAD.Ded;
Kn=preMPFAD.Kn;
Kt=preMPFAD.Kt;
Hesq=preMPFAD.Hesq;
nflag=preMPFAD.nflag;

bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

flowrate    = zeros(bedgesize+inedgesize,1);
flowratedif = zeros(bedgesize+inedgesize,1);

flowresult = zeros(size(centelem,1),1);

%% ==================================================
%% BOUNDARY EDGES  (VETORIZADO)
%% ==================================================

B1 = bedge(:,1);
B2 = bedge(:,2);
lef = bedge(:,3);

coordB1 = coord(B1,:);
coordB2 = coord(B2,:);

edgevec = coordB1 - coordB2;
nor = sqrt(sum(edgevec.^2,2));

O = centelem(lef,:);

dirmask = bedge(:,5) < 200;
neumask = ~dirmask;

c1 = nflag(B1,2);
c2 = nflag(B2,2);

A = Kn ./ (Hesq .* nor);

term1 = sum((O-coordB2).*(coordB1-coordB2),2).*c1;
term2 = sum((O-coordB1).*(coordB2-coordB1),2).*c2;

flowrate_b = -A.*(term1+term2-(nor.^2).*p(lef)) -(c2-c1).*Kt;

%% viscosity
visonface = ones(bedgesize,1);

if 30<numcase && numcase<200
    visonface = sum(viscosity(1:bedgesize,:),2);
elseif 200<numcase && numcase<300
    if any(numcase == [246 245 247 248 249 251])
        visonface = viscosity(1:bedgesize,:);
    end
end

flowrate_b = visonface .* flowrate_b;

%% Neumann boundary
if any(neumask)
    bcvals = zeros(bedgesize,1);
    [~,loc]=ismember(bedge(neumask,5),bcflag(:,1));
    bcvals(neumask)=bcflag(loc,2);
    flowrate_b(neumask) = -nor(neumask).*bcvals(neumask);
end

flowrate(1:bedgesize)=flowrate_b;

%% ==================================================
%% INTERNAL EDGES  (VETORIZADO)
%% ==================================================

node1 = inedge(:,1);
node2 = inedge(:,2);

lef = inedge(:,3);
rel = inedge(:,4);

p1 = pinterp(node1);
p2 = pinterp(node2);

visonface_i = ones(inedgesize,1);

if 30<numcase && numcase<200
    visonface_i = sum(viscosity(bedgesize+1:bedgesize+inedgesize,:),2);
elseif 200<numcase && numcase<300
    if any(numcase == [246 245 247 248 249 251])
        visonface_i = viscosity(bedgesize+1:bedgesize+inedgesize,:);
    end
end

flowrate(bedgesize+1:end,1) = visonface_i .* Kde .* (p(rel)-p(lef)-Ded.*(p2-p1));

if numcase==435 || numcase==431 || numcase==437
    flowrate = flowrate + flowrateZ;
end

%% ==================================================
%% FLOWRESULT (ACCUMARRAY)
%% ==================================================
idx = bedgesize + (1:inedgesize);
auxlef=bedge(:,3);

flowresult = flowresult + accumarray(auxlef,flowrate(1:bedgesize),size(flowresult));


flowresult = flowresult + accumarray(lef,flowrate(idx),size(flowresult))- ...
            accumarray(rel,flowrate(idx),size(flowresult));

%% ==================================================
%% DISPERSIVE FLOW
%% ==================================================

if (200<numcase && numcase<300) || (379<numcase && numcase<400)

    con1 = cinterp(node1);
    con2 = cinterp(node2);

    flowratedif(bedgesize+1:end) = ...
        Kdec .* (Con(rel)-Con(lef)-Dedc.*(con2-con1));

end

faceaux=0;

end