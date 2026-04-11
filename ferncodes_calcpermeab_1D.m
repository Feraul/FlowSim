

function permeabi=ferncodes_calcpermeab_1D
global Nmod varK centelem
[phiGaussNmod10000,wavenumberGauss0Nmod10000,wavenumberGauss1Nmod10000]=parametrosGauss_1D;
%[phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000]=parametrosExpo;
phi = phiGaussNmod10000(1:Nmod);
C(:,1) = wavenumberGauss0Nmod10000(1:Nmod);
C(:,2) = wavenumberGauss1Nmod10000(1:Nmod);
KMean = 15;

permeabi=zeros(size(centelem,1),1);
for i = 1 : size(centelem,1)
    permeabi(i,1) = K(centelem(i,1),Nmod,KMean,varK,C(:,1),C(:,2),phi);
end
end

function sol = K(x,Nmod,KMean,varK,C1,C2,phi)
coeff = sqrt(varK*2/Nmod);
ak = coeff*sum(cos( (C1*x + C2)*(2*pi) + phi));
sol = KMean * exp(-varK/2) * exp(ak);
end