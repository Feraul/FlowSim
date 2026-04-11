%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 03/05/2012
%Modify data:  / /2012
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%This FUNCTION is called by "calctimestep.m", "getmobility.m", 
%"calcnumflux.m" 

%--------------------------------------------------------------------------

function [fw,fo,gama,krw,kro,nw,no,kwmax,komax] = twophasevar2(Sw,...
    numcase,overedgecoord,injecelem)
%Define global parameters:
global centelem satlimit visc;

%Initialize the parameters:
krw = zeros(length(Sw),1);
kro = krw;
fw = krw;
fo = krw;

%"gama" is equal to "fw*kro/visc(2)"
gama = krw;

%Calculate the parameters to each element
    %Special case (Hurtado et al., 2007)
    %Cases 34.3 and 34.4 ==> 1 CV sourse; 46.3 ==> 4 CV source 
    %(boundary free)
    if numcase == 34.3 || numcase == 34.4
        %Define Mobility ratio:
        M = visc(2)/visc(1);
        %Definition of fractional flow (WATER)
        fw = Sw.^2;
        %Definition of fractional flow (OIL)
        fo = 1 - fw;
        %Definition of relative permeability (WATER)
        krw = fw./(M.*(1 - fw) + fw); 
        %Definition of relative permeability (OIL)
        kro = 1 - krw;
    
    %Adapted from Bastian (2002) for the benchmark 31.1, with lambda = 2; 
    elseif numcase == 31.1 || numcase == 31.6 
        %Normalizes the saturation:
        Swn = ((Sw - satlimit(1))./(1 - satlimit(1) - satlimit(2))); 
        
        %Definition of relative permeability (WATER)
        krw = Swn.^4; 
        %Definition of relative permeability (OIL)
        kro = ((1 - Swn).^2).*(1 - (Swn.^2));

        %------------------------------------------------------------------
        %Fractional Flow (equal to all cases)
    
        %Definition of fractional flow (WATER)
        fw = (krw./visc(1))./((krw./visc(1)) + (kro./visc(2)));
        %Definition of fractional flow (OIL)
        fo = (kro./visc(2))./((krw./visc(1)) + (kro./visc(2)));
    
        %------------------------------------------------------------------
        %Define "gama". It is used when gravity effects are account
    
        gama = fw.*kro./visc(2);
        
    %Bucley-Leverett with two materials. (Yu-Shu Wu, 2017) - Livro
    %Evaluate "x" coordinate in face or cnetroid coordinate in order to
    %define the material (rock)
    elseif numcase == 31.9 && length(Sw) ~= length(injecelem)
%         Swi1 = 0.15;
%         Sor1 = 0.1;
%         Swi2 = 0.15;
%         Sor2 = 0.15;
%         %Brooks-Corey expoent
%         nw1 = 2.5;
%         no1 = 1.0;
%         nw2 = 2.0;
%         no2 = 3.0;
%         %Maximum and minimum parameters
%         kwmax1 = 0.9;
%         komax1 = 0.9;
%         kwmax2 = 0.8;
%         komax2 = 0.8;

        Swi2 = 0.15;
        Sor2 = 0.10;
        Swi1 = 0.15;
        Sor1 = 0.10;
        %Brooks-Corey expoent
        nw2 = 2.5;
        no2 = 1.0;
        nw1 = 2.0;
        no1 = 1.5;
        %Maximum and minimum parameters
        kwmax2 = 0.9;
        komax2 = 0.9;
        kwmax1 = 0.8;
        komax1 = 0.8;
           
        %Define an auxiliary position vector of "Sw"
        auxpos = (1:length(Sw))';
        %Points the limits of each region
        %Region 1:
        %Evaluate the faces
        if length(Sw) == size(centelem,1)
            pntreg1 = centelem(:,1) < 60;
        %Evaluate for the centroids
        else
            pntreg1 = overedgecoord(:,1) < 60;
        end  %End of IF
        
        %Get the properties:
        [fw1,fo1,gama1,krw1,kro1] = getproperties(Sw(pntreg1),Swi1,Sor1,...
            nw1,no1,kwmax1,komax1,visc(1),visc(2));
        %Region 2:
        %Get the properties:
        [fw2,fo2,gama2,krw2,kro2] = getproperties(Sw(~pntreg1),Swi2,Sor2,...
            nw2,no2,kwmax2,komax2,visc(1),visc(2));
        %Join both regions
        fw(auxpos(pntreg1)) = fw1; 
        fw(auxpos(~pntreg1)) = fw2;
        fo(auxpos(pntreg1)) = fo1; 
        fo(auxpos(~pntreg1)) = fo2;
        gama(auxpos(pntreg1)) = gama1; 
        gama(auxpos(~pntreg1)) = gama2;
        krw(auxpos(pntreg1)) = krw1; 
        krw(auxpos(~pntreg1)) = krw2;
        kro(auxpos(pntreg1)) = kro1; 
        kro(auxpos(~pntreg1)) = kro2;

    %Boundary in Buckley-Leverett problem
    elseif numcase == 31.9 && length(Sw) == length(injecelem)
%         Swi1 = 0.15;
%         Sor1 = 0.1;
        Swi1 = 0.15;
        Sor1 = 0.10;
        %Brooks-Corey expoent
%         nw1 = 2.5;
%         no1 = 1.0;
        nw1 = 2.0;
        no1 = 1.5;
        %Maximum and minimum parameters
%         kwmax1 = 0.9;
%         komax1 = 0.9;
        kwmax1 = 0.8;
        komax1 = 0.8;
           
        %Get the properties:
        [fw,fo,gama,krw,kro] = getproperties(Sw,Swi1,Sor1,...
            nw1,no1,kwmax1,komax1,visc(1),visc(2));

    %Another examples
    else
        %It Chose the properties according to "numcase"
        %Example 36: Two-Phase Flow case. Lamine and Edwards, 2010.
        booleanlin = (numcase == 36 || numcase == 36.1);
        %Example 48: Adapted from "Khoozan et al. (2011)" (stearscase)
        booleanmean = numcase > 48 && numcase < 49; 
        boolean3rd = (numcase == 35);
        %Example 45: Two-Phase Flow case. Kozdon et al., 2011. 
        %"Multidimensional upstream weighting for multiphase transport 
        %in porous media". Section 5.2 the benchmark 45 to 46.
        boolean4th = ((numcase >= 45 && numcase < 46) || ...
            numcase == 43.1);
        %Example 34.8: Two-Phase Flow case. Adapted from Nikitin and 
        %Vassilevisky, 2012. 
        %It evaluates three different quadrangular meshes.
        %Pinchout problem (Li and Riviere, 2015)
        boolean5th = (numcase == 34.8) || (numcase == 34.9);
        %Example 37.2: Flow in anisotropic media with two also anisotropic 
        %pinchouts (just devised)
        booleanhalf = numcase == 37.2;
        %Any another case
        boolean2nd = (booleanmean + booleanlin + boolean3rd + boolean4th + ...
            boolean5th + booleanhalf) == 0; 
        %Get the exponent:
        nw = 4*booleanmean + booleanlin + 3*boolean3rd + 4*boolean4th + ...
            5*boolean5th + 2*boolean2nd + 2.5*booleanhalf; 
        no = 2*booleanmean + booleanlin + boolean3rd + 2*boolean4th + ...
            boolean5th + 2*boolean2nd + 3*booleanhalf;
        
        %Fit parameter (water and oil)
        kwmax = (1 - (booleanmean + boolean3rd + booleanhalf)) + ...
            (booleanmean + boolean3rd)*0.6 + booleanhalf*1;
        komax = (1 - booleanhalf) + booleanhalf*1;

%             kwmax = 0.6;
        
        %------------------------------------------------------------------
        %Normalized Saturation
    
        Swn = ((Sw - satlimit(1))./(1 - satlimit(1) - satlimit(2))); 
    
        %------------------------------------------------------------------
        %Relative Permeability:
    
        %Definition of relative permeability (WATER)
        krw = kwmax*(Swn).^nw; 
        %Definition of relative permeability (OIL)
        kro = komax*(1 - Swn).^no; 
    
        %------------------------------------------------------------------
        %Fractional Flow (equal to all cases)
    
        %Definition of fractional flow (WATER)
        fw = (krw./visc(1))./((krw./visc(1)) + (kro./visc(2)));
        %Definition of fractional flow (OIL)
        fo = (kro./visc(2))./((krw./visc(1)) + (kro./visc(2)));
    
        %------------------------------------------------------------------
        %Define "gama". It is used when gravity effects are account
    
        gama = fw.*kro./visc(2);
    end  %End of IF (special case)

%--------------------------------------------------------------------------
%Function "getproperties"
%--------------------------------------------------------------------------

function [fw,fo,gama,krw,kro] = getproperties(Sw,Swi,Sor,nw,no,kwmax,komax,...
    miw,mio)            
Swn = ((Sw - Swi)./(1 - Swi - Sor)); 
%Definition of relative permeability (WATER)
krw = kwmax*(Swn).^nw; 
%Definition of relative permeability (OIL)
kro = komax*(1 - Swn).^no; 
    
%Fractional Flow (equal to all cases)
fw = (krw./miw)./((krw./miw) + (kro./mio));
%Definition of fractional flow (OIL)
fo = (kro./mio)./((krw./miw) + (kro./mio));
    
%Define "gama". It is used when gravity effects are account
gama = fw.*kro./mio;
    