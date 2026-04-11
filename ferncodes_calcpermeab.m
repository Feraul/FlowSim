function[permeabi]=ferncodes_calcpermeab

global Nmod centelem
%% 2D exponential case: computation of K fields and source terms f

[phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000]=parametrosGauss;
%[phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000]=parametrosExpo;
phi = phiExpNmod10000(1:Nmod);
C(:,1) = wavenumberExp0Nmod10000(1:Nmod);
C(:,2) = wavenumberExp1Nmod10000(1:Nmod);
KMean = 15;

permeabi=zeros(size(centelem,1),1);
for i = 1 : size(centelem,1)
    permeabi(i,1) = K(centelem(i,1),centelem(i,2),KMean,C(:,1),C(:,2),phi);
end
end
function sol = K(x,y,KMean,C1,C2,phi)
global Nmod varK
coeff = sqrt(varK*2/Nmod) ;

ak = coeff*sum( cos( (C1*x + C2*y)*(2*pi) + phi) ) ;

sol = KMean * exp(-varK/2) * exp(ak) ;

end

