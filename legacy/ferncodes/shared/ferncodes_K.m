function sol = ferncodes_K(x,y)
global varK Nmod

[phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000]=parametrosGauss;
%[phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000]=parametrosExpo;
phi = phiExpNmod10000(1:Nmod);
C1(:,1) = wavenumberExp0Nmod10000(1:Nmod);
C2(:,1) = wavenumberExp1Nmod10000(1:Nmod);
KMean = 15;

coeff = sqrt(varK*2/Nmod) ;

ak = coeff*sum( cos( (C1*x + C2*y)*(2*pi) + phi) ) ;

sol = KMean * exp(-varK/2) * exp(ak) ;

end