%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical code used to simulate fluid flow in porous media. That 
%routine calls several others which defines how the equation will be solved 
%Type of file: FUNCTION
%Criate date: 13/03/2013
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Calculate the water saturation field using an explicit formulation. 

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

function [d2fdS,dgamadS] = calcder2dS(fw,gama,Sw,dertype,overedgecoord)
%Define global parameters
global centelem satlimit visc numcase keygravity;

%Define "d2fdS" and "dgamadS" according to "dertype"
switch dertype
    %Discrete derivative. It needs of two saturation value 
    %(used in timestep definition, for example)
    case 0
        %There saturation difference bigger than zero
        if (Sw(1) - Sw(2)) ~= 0
            d2fdS = (fw(1) - fw(2))/(Sw(1) - Sw(2));
            dgamadS = (gama(1) - gama(2))/(Sw(1) - Sw(2));
        %The saturation difference is zeros
        else
            d2fdS = 0;
            dgamadS = 0;
        end  %End of IF
    %Analitical derivative. It uses only one value into argument.
    case 1
        %Initialize some propoerties (two-phase flow)
        Swi = satlimit(1);
        Sor = satlimit(2);
        miw = visc(1);
        mio = visc(2);
        %Initialize "d2fdS" and "dgamadS"
        d2fdS = zeros(length(Sw),1);
        dgamadS = d2fdS; 
        
        %------------------------------------------------------------------
        %Chose analitical derivative according to "numcase"
        
        %Adapted from Bastian (2002) for the benchmark 31.1 
        %(van Genushten model); 
        if numcase == 31.1 || numcase == 31.6 
            %Calculate "d2fdS"
            d2fdS = ...
                -(4*(Sw.^2).*miw.*mio.*((Sw.^4).*(5 - 6.*Sw + (Sw.^3))*mio ...
                - ((Sw - 1).^4).*(3 + 4*Sw + 4.*(Sw.^2) + (Sw.^3))*miw))./...
                (((Sw.^4)*mio - ((Sw - 1).^3).*(Sw + 1)*miw).^3);
        %Hurtado et al. (2007)
        elseif numcase == 34.3
            d2fdS = 2;
            dgamadS = 0;
        %Bucley-Leverett with two materials. (Yu-Shu Wu, 2017) - Livro
        elseif numcase == 31.9
%             Swi1 = 0.15;
%             Sor1 = 0.1;
%             Swi2 = 0.15;
%             Sor2 = 0.15;
%             %Brooks-Corey expoent
%             nw1 = 2.5;
%             no1 = 1.0;
%             nw2 = 2.0;
%             no2 = 3.0;
%             %Maximum and minimum parameters
%             kwmax1 = 0.9;
%             komax1 = 0.9;
%             kwmax2 = 0.8;
%             komax2 = 0.8;
% 
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
            if  length(Sw) == size(centelem,1)
                pntreg1 = centelem(:,1) < 60;
            else
                pntreg1 = overedgecoord(:,1) < 60;
            end  %End of IF
            
            %Get the second derivative:
            [d2fdS1] = getder(Sw(pntreg1),Swi1,Sor1,nw1,no1,kwmax1,komax1,...
                miw,mio);
            %Region 2:
            %Get the second derivative:
            [d2fdS2] = getder(Sw(~pntreg1),Swi2,Sor2,nw2,no2,kwmax2,komax2,...
                miw,mio);
            %Join both regions
            d2fdS(auxpos(pntreg1)) = d2fdS1; 
            d2fdS(auxpos(~pntreg1)) = d2fdS2;
        
        %Another examples. (Brooks-Corey model). It can vary the 
        %coefficients "nw", "no", ...
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
            %Example 37.2: Flow in anisotropic media with two also 
            %anisotropic pinchouts (just devised)
            booleanhalf = numcase == 37.2;
            %Any another case
            boolean2nd = (booleanmean + booleanlin + boolean3rd + ...
                boolean4th + boolean5th + booleanhalf) == 0; 
            %Get the exponent:
            nw = 4*booleanmean + booleanlin + 3*boolean3rd + ...
                4*boolean4th + 5*boolean5th + 2*boolean2nd + ...
                2.5*booleanhalf; 
            no = 2*booleanmean + booleanlin + boolean3rd + 2*boolean4th + ...
                boolean5th + 2*boolean2nd + 3*booleanhalf;
        
            %Fit parameter (water and oil)
            kwmax = (1 - (booleanmean + boolean3rd + booleanhalf)) + ...
                (booleanmean + boolean3rd)*0.6 + booleanhalf*1;
            komax = (1 - booleanhalf) + booleanhalf*1;

            nw = 4;
            no = 2;
%             kwmax = 0.6;

            %Define some terms:
            term1 = (Sw + Sor - 1);
            term2 = (Swi + Sor - 1);
            term3 = (Swi - Sw);
            term4 = (term1./term2).^no;
            term5 = (term3./term2).^nw;
            term6 = (komax*miw.*term4 + kwmax*mio.*term5);
            denom = ((term1.^2).*((-term3).^2).*(term6.^3)) + 1e-16; 
            %Calculate "dfwdSw"
            d2fdS = (komax*kwmax*miw*mio*term4.*...
                term5.*(-(nw - no).*term1.*term3.*term6 - ...
                term1.*(nw.*term1 + no.*term3).*term6 + nw.*term1.*...
                (nw.*term1 + no.*term3).*term6 + term3.*(nw.*term1 + ...
                no.*term3).*term6 - no*term3.*(nw*term1 + no*term3).*...
                term6 + 2*term1.*term3.*(nw.*term1 + ...
                no.*term3).*(((komax*miw*no*term4)./term1) - ...
                ((kwmax*mio*nw.*term5)./term3))))./denom;

            %Calculate the gravity contribution (when necessary)
            if strcmp(keygravity,'y')
                %Define another terms:
                %Calculate "dgamadSw"
                dgamadS = (komax*kwmax*term4.*term5.*(komax*miw*nw*term1.*...
                    term4 - kwmax*mio*no*term3.*term5))./denom;
            end  %End of IF
        end  %End of IF (is the Bastian example?)
end  %End of SWITCH


%--------------------------------------------------------------------------
%Function "getder"
%--------------------------------------------------------------------------

function [d2fdS] = getder(Sw,Swi,Sor,nw,no,kwmax,komax,miw,mio)            
%Define some terms:
term1 = (Sw + Sor - 1);
term2 = (Swi + Sor - 1);
term3 = (Swi - Sw);
term4 = (term1./term2).^no;
term5 = (term3./term2).^nw;
term6 = (komax*miw.*term4 + kwmax*mio.*term5);
denom = ((term1.^2).*((-term3).^2).*(term6.^3)) + 1e-16; 
%Calculate "dfwdSw"
d2fdS = (komax*kwmax*miw*mio*term4.*term5.*(-(nw - no).*term1.*term3.*...
    term6 - term1.*(nw.*term1 + no.*term3).*term6 + nw.*term1.*...
    (nw.*term1 + no.*term3).*term6 + term3.*(nw.*term1 + no.*term3).*...
    term6 - no*term3.*(nw*term1 + no*term3).*term6 + 2*term1.*term3.*...
    (nw.*term1 + no.*term3).*(((komax*miw*no*term4)./term1) - ...
    ((kwmax*mio*nw.*term5)./term3))))./denom;
