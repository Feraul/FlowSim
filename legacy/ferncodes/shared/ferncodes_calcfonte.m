function[fonte]=ferncodes_calcfonte

global Nmod centelem
%% 2D exponential case: computation of K fields and source terms f

[phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000]=parametrosGauss;
%[phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000]=parametrosExpo;
phi = phiExpNmod10000(1:Nmod);
C(:,1) = wavenumberExp0Nmod10000(1:Nmod);
C(:,2) = wavenumberExp1Nmod10000(1:Nmod);
KMean = 15;

fonte=zeros(size(centelem,1),1);
for i = 1 : size(centelem,1)
    fonte(i,1) = func(centelem(i,1),centelem(i,2),KMean,C(:,1),C(:,2),phi);
end
end
function F = func(x,y,KMean,C1,C2,phi)
global Nmod varK
S1 = sum( (-2*pi)*C1.*sin(phi + (2*pi)*(C1*x + C2*y))) ;

S2 = sum( cos(phi + (2*pi)*(C1*x + C2*y) ) ) ;

S3 = sum( (-2*pi)*C2.*sin(phi + (2*pi)*(C1*x + C2*y))) ;

F = -(2*KMean*exp(-varK/2)*sqrt(varK*2/Nmod)*S1*exp(sqrt(varK*2/Nmod)*S2)*cos(2*x+y) ...
    - 5*KMean*exp(-varK/2)*exp(sqrt(varK*2/Nmod)*S2)*sin(2*x+y) ...
    + KMean*exp(-varK/2)*sqrt(varK*2/Nmod)*S3*exp(sqrt(varK*2/Nmod)*S2)*cos(2*x+y));

end


