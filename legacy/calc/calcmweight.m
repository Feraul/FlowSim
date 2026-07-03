%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 05/03/2014 (Ash Wednesday)
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: 

%--------------------------------------------------------------------------
%Additional comments: 
%1. It is called by the function "getmultidsatweig" and 
%"getMassWonhalfedge";
%2. "pointrow" is the position of the ADJACENT half-edge in a vector with 
%the "bedge" size + "inedge" size.
%3. Do NOT confuse "pointrow" with "frposition". "frposition" is the
%position of the EVALUATED half-edge in a vector with 2*("bedge" size + 
%"inedge" size). That is, a vector with the amount of half-edges.
%4. "otherposition" is equivalent to "frposition" but for the ADJACENT
%half-edge.

%--------------------------------------------------------------------------

function [mweight] = calcmweight(elemeval,isleft,localflowrate,...
    upstrflowrate,multdlimiter,evalvec,upstrvec)

% fw_adj = fw_adj*(abs(fw_adj) > 1e-12);
% fw_upst = fw_upst*(abs(fw_upst) > 1e-12);

%Flow rate evaluated
freval = abs(localflowrate);

%Define "mweight" by using the sign of flowrate.

%The evaluated flow rate is bigger than zero but the flow rate secundary is 
%lower than zero
if (upstrflowrate < 0 && elemeval ~= isleft) || ...
        (upstrflowrate > 0 && elemeval == isleft) || ...
        (freval == 0) || (upstrflowrate == 0)
    %Attribute "0" to "mweight" (only the cell-centered saturation).
    mweight = 0;
%Any other configuration
else
    %Calculate the weight by a rate (linear combination).
    mweight = abs(upstrflowrate)/(freval + 1e-16);
end  %End of IF

%Choose the limiter (for the convex combination) according to 
%"multdlimiter" value
switch multdlimiter
    %TMU (Tran et al., 2005)
    case 1
        mweight = min(1,mweight);

    %SMU (Hurtado et al., 2007)
    case 2
        mweight = mweight/(1 + mweight);

    otherwise
        angvec = ...
            acosd(dot(evalvec,upstrvec)/(norm(evalvec)*norm(upstrvec)));
        epsilon = abs(angvec/90);
        
%        func1 = 2 - epsilon;
        func1 = 1 - nthroot((epsilon - 1),3);
%        func2 = epsilon;%1 - nthroot((1 - epsilon),3);
%             func2 = 1 - nthroot((1 - epsilon),3);
%             func = 1 - nthroot((epsilon - 1),7);
%             func = 1 - ((epsilon - 1)^3);
%         mweight = mweight/(1 + mweight);

        %Define SMU limiter
        mw_SMU = mweight/(1 + mweight);

        %SMU modified
        if multdlimiter == 3
            mweight = min(1,func1*mw_SMU);
        %Blend of SMU and TMU
        elseif multdlimiter == 4
            %Define TMU limiter
            mw_TMU = min(1,mweight);
            mweight = min(1,max(mw_TMU,func1*mw_SMU));
        end  %End of IF
end  %End of SWITCH
        
