
function [fonte]=ferncodes_calcfonte_1D
global centelem Nmod varK

[phiGaussNmod10000,wavenumberGauss0Nmod10000,wavenumberGauss1Nmod10000]=parametrosGauss_1D;
phi = phiGaussNmod10000(1:Nmod)  ;
C(:,1) = wavenumberGauss0Nmod10000(1:Nmod);
C(:,2) = wavenumberGauss1Nmod10000(1:Nmod);
KMean = 15;

fonte=zeros(size(centelem,1),1);
for i = 1 : size(centelem,1)
    fonte(i,1) = func(centelem(i,1),Nmod,KMean,varK,C(:,1),C(:,2),phi);
end
end


function f = func(x,Nmod,KMean,varK,wave1,wave2,phi)

S1 = sum( (-2*pi)*wave1.*sin(phi + 2*pi*(wave1*x + wave2))) ;
S2 = sum( cos(phi + 2*pi*(wave1*x + wave2))) ;
f = -(KMean*exp(-varK/2)*exp(sqrt(varK*2/Nmod)*S2)*(sqrt(varK*2/Nmod)*S1*cos(x) - sin(x)) );

end