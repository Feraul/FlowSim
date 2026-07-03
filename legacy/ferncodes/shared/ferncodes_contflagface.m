function [nflag]= ferncodes_contflagface
 
% determinar o flag do nó interior e fronteira de Neumann
global bedge bcflag;
nflag = ones(size(bedge,1),2);

for ifacont = 1:size(bedge,1)
    x = logical(bcflag(:,1) == bedge(ifacont,5));
    vertex = bedge(ifacont,1:2);
    %Second column receives the boundary condition value.
    nflag(ifacont,2) = PLUG_bcfunction(vertex,x,0);
    %First column receives the boundary condition flag.
    nflag(ifacont,1) = bcflag(x,1);
end  %End of FOR
end