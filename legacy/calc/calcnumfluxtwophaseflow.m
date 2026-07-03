%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with
%high-order resolution
%Type of file: FUNCTION
%Criate date: 13/01/2013
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%.

%--------------------------------------------------------------------------
%Additional comments:
%

%--------------------------------------------------------------------------

function [advecterm,entrineqterm,earlysw,vectorSleft,vectorSright] = calcnumfluxtwophaseflow(Sw,Fg,flowrate,...
    taylorterms,limiterflag,flagknownvert,satonvertices,flagknownedge,...
    satonboundedges,pointbndedg,pointinedg,orderbedgdist,orderinedgdist,...
    constraint,mlplimiter,earlysw,countinter)

                
%Define global parameters:
global coord elem bedge inedge normals dens numcase centelem ;

%Initialize a tolerance. It is a computational zero
tol = 1e-9;
%Initialize "advecterm" and "entrineqterm"
advecterm = zeros(size(elem,1),1);
entrineqterm = advecterm;
%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);

lenginjface = 1.166822;

%--------------------------------------------------------------------------
%Boundary edges (when it exists):

%In cases where the producer well edges are evaluated, it is necessary
%verify if the producer well shares some boundary edge.
if any(pointbndedg)
    %Swept "bedge"
    for i = 1:length(pointbndedg)
        %Initialize some parameters:
        ibedg = pointbndedg(i);
        %Define the order for this edge.
        faceorder = orderbedgdist(i);
        
        %Define the "vertices"
        vertices = bedge(ibedg,1:2);
        %Get the coordinate for the vertices:
        verticescoord = coord(vertices,:);
        %Define left elements
        leftelem = bedge(ibedg,3);
        
        %------------------------------------------------------------------
        %Define velocity due gravity
        
        %There is gravity
        if size(Fg,1) > 1
            %dotvg = Fg(ibedg,1)*(dens(1) - dens(2))/lenginjface;
             dotvg = Fg(ibedg,1)*(dens(1) - dens(2));
            %There is NO gravity
        else
            dotvg = 0;
        end  %End of IF
        
        %------------------------------------------------------------------
        
        %Define the elements that share the edge evaluated
        elemeval = leftelem;
        
        %Verify if there is saturation prescribed on boundary:
        %There is a prescribed saturation
        if flagknownedge(ibedg) == 1
            %Attribute the saturation on boundary
            Sleft = satonboundedges(ibedg);
            if Sleft>1
                Sleft=1;
                
            end
            Sleft = Sleft*(Sleft >= 0);
            %There is no prescribed saturation. It is necessary calculate.
        else
            %Get the statment regarding to mlp limiter:
            boolmlp = (length(mlplimiter) == 1);
            mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
            
            %Get the saturation value recovered
            [Sleftlim,Sleftnonlim] = getsatonedge(elemeval,vertices,verticescoord,...
                taylorterms,Sw,limiterflag,faceorder,constraint,flagknownvert,...
                satonvertices,mlpbyelem,centelem(leftelem,1:2));
            
            if Sleftlim>1
                %Sleftlim=Sw(leftelem);
                Sleftlim=1;
            elseif Sleftlim<0
                %Sleftlim=Sw(leftelem);
                Sleftlim=0;
            end
            if strcmp(limiterflag{9},'on')|| strcmp(limiterflag{11},'on')|| strcmp(limiterflag{12},'on')
                
                if Sleftnonlim>1 || Sleftnonlim<0
                    Sleft=Sleftlim;
                    
                else
                   Sleft= Sleftnonlim;
                end
                % %==========================================================================
                if Sleftnonlim>Sw(leftelem)
                    Sleft=Sleftlim;
                else
                   Sleft= Sleftnonlim; 
                end
            else
                Sleft=Sleftlim;
            end
            %==========================================================================
        end  %End of IF
        
%         Sleft = Sleft*(Sleft >=0);
%         Sleft = Sleft*(abs(Sleft) > tol);
        %Fill "earlysw"
        earlysw(ibedg) = Sleft;
        
        %Calculate the fractional flow in boundary ("fwbound")
        [~,fw,~,gama,] = twophasevar(Sleft);
        
        %Define the normal velocity into face
        dotvn = flowrate(ibedg);
        %Get accuracy for "dotvn"
        dotvn = dotvn*(abs(dotvn) > tol);
                
        %Calculate the numerical flux through interface.
        numflux = dotvn*fw + dotvg*gama;
        %Obtain the contribution of interface over element to LEFT
        advecterm(leftelem) = advecterm(leftelem) + numflux;
        vectorSleft(i,1)=Sleft;
    end  %End of FOR (Swept "bedge")
end  %End of IF (Does evaluate the boundary edges?)

%--------------------------------------------------------------------------
%Internal edges:

%Swept "inedge" evaluating left and right elements by edge. Apply
%approximated Riemann Solver through edge.

for i = 1:length(pointinedg)
    
    mLLF=1;
    
    %Initialize some parameters:
    inedg = pointinedg(i);
    %---------------------------------------
    %Define the normal velocity in each face
    dotvn = flowrate(bedgesize + inedg);
    %Define "vertices"
    vertices = inedge(inedg,1:2);
    v1=coord(vertices,1:2);
    medio=(v1(1,:)+v1(2,:))/2;
    %Get the coordinate for the vertices:
    verticescoord = coord(vertices,:);
    %Define left and right elements
    leftelem = inedge(inedg,3);
    rightelem = inedge(inedg,4);
    
    rijL= medio-centelem(leftelem,1:2);
    rijR= centelem(rightelem,1:2)-medio;
    
    %Left Contribution:
    %Define the order for this edge.
    faceorder = orderinedgdist(i,1);
    %Define the elements that share the edge evaluated
    elemeval = [leftelem rightelem];
    %----------------------------------------------------------------------
    %Calculate the velocity due to GRAVITY effect
    
    %There is gravity
    if size(Fg,1) > 1
        %dotvg = Fg(bedgesize + inedg,1)*(dens(1) - dens(2))/lenginjface;
        dotvg = Fg(bedgesize + inedg,1)*(dens(1) - dens(2));
        %There is NO gravity
    else
        dotvg = 0;
    end  %End of IF
    %----------------------------------------------------------------------
    
    %Get the saturation value recovered on each quadrature point ("on_q")
    %Get the statment regarding to mlp limiter:
    boolmlp = (length(mlplimiter) == 1);
    mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
    %mlpbyelem2=mlplimiter(boolmlp + (1 - boolmlp)*elemeval(2),:);
    %Left Contribution:
     [Sleftlim,Sleftnonlim] = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,...
        limiterflag,faceorder,constraint,flagknownvert,satonvertices,...
        mlpbyelem,centelem(elemeval,1:2));
    
    %Right Contribution:
    %Define the order for this edge.
    faceorder = orderinedgdist(i,2);
    %Define the elements that share the edge evaluated
    elemeval = [rightelem leftelem];
    %Get the statment regarding to mlp limiter:
    mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
    %mlpbyelem2=mlplimiter(boolmlp + (1 - boolmlp)*elemeval(2),:);
    %Get the saturation value recovered on each quadrature point ("on_q")
    [Srightlim,Srightnonlim]= getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,...
        limiterflag,faceorder,constraint,flagknownvert,satonvertices,...
        mlpbyelem,centelem(elemeval,1:2));
   if Sleftlim>1
        %Sleftlim=Sw(leftelem);
        Sleftlim=1;
    elseif Sleftlim<0
        %Sleftlim=Sw(leftelem);
        Sleftlim=0;
    end
    if Srightlim>1
        %Srightlim=Sw(rightelem);
        Srightlim=1;
    elseif Srightlim<0
        %Srightlim=Sw(rightelem);
        Srightlim=0;
    end
    %PAD
    if (strcmp(limiterflag{9},'on')|| strcmp(limiterflag{11},'on')|| strcmp(limiterflag{12},'on') ) && (countinter==0)
        [Sleft,Sright,mLLF]=PhysicalAD(Sw,taylorterms,limiterflag,flagknownvert,satonvertices,...
            constraint,mlplimiter,leftelem,rightelem,vertices,...
            verticescoord,Sleftlim,Sleftnonlim,Srightlim,Srightnonlim,dotvn,tol,rijL,rijR);       
        
    else
        Sleft= Sleftlim;
        Sright=Srightlim;
    end
    %%
%      Sleft = Sleft*(Sleft >=0);
%      Sleft = Sleft*(abs(Sleft) > tol);
%      Sright = Sright*(Sright >=0);
%      Sright = Sright*(abs(Sright) > tol);

    vectorSleft(size(bedge,1)+i,1)=Sleft;
    vectorSright(i,1)=Sright;
    %Discrete:
    %"fw" has three values: fw(Sleft) is fw(1), fw(Sright) is fw(3)
    [~,fw,~,gama,] = twophasevar([Sleft 0.5*(Sleft+Sright) Sright]);
   
    % fw(1) ---> fluxo fracional no elemento fw(Sleft)
    % fw(2) ---> fluxo fracional no elemento fw(Smid)
    % fw(3) ---> fluxo fracional no elemento fw(Sright)
    %Calculate the Rankine-Hugoniot ratio:
    % veja o tese Darlan pag. 165 e sempre é positivo
    [dfwdS_rh,dgamadS_rh] = ...
        calcdfunctiondS([fw(1) fw(3)],[gama(1) gama(3)],[Sleft Sright],0);
    
    %Get accuracy (RH ratio):
    dfwdS_rh = dfwdS_rh*(abs(dfwdS_rh) > tol);
    %Calculate the derivative dfw/dSw:
    %Analytical:
    [dfwdS_analeft,dgamadS_analeft] = calcdfunctiondS(0,0,Sleft,1);
    [dfwdS_analright,dgamadS_analright] = calcdfunctiondS(0,0,Sright,1);
    
    %Get accuracy (second derivative):
    dfwdS_analeft = dfwdS_analeft*(abs(dfwdS_analeft) > tol);
    dfwdS_analright = dfwdS_analright*(abs(dfwdS_analright) > tol);
    
    %Calculate the second derivative d2fw/dSw2:
    [d2fwdS2_discleft,d2gamadS2_discleft] = calcder2dS(0,0,Sleft,1);
    [d2fwdS2_discright,d2gamadS2_discright] = calcder2dS(0,0,Sright,1);
    %Get accuracy (second derivative):
    d2fwdS2_discleft = d2fwdS2_discleft*(abs(d2fwdS2_discleft) > tol);
    d2fwdS2_discright = d2fwdS2_discright*(abs(d2fwdS2_discright) > tol);
    d2gamadS2_discleft = d2gamadS2_discleft*(abs(d2gamadS2_discleft) > tol);
    d2gamadS2_discright = d2gamadS2_discright*(abs(d2gamadS2_discright) > tol);
    
    %Get accuracy:
    dotvn = dotvn*(abs(dotvn) > tol);
    dotvg = dotvg*(abs(dotvg) > tol);
    %---------------------------------------
    
    %Get the sign of the first derivative:
    signder_left = sign(dfwdS_analeft*dotvn + dgamadS_analeft*dotvg);
    signder_right = sign(dfwdS_analright*dotvn + dgamadS_analright*dotvg);
    %Get the sign of second derivative:
    sign2der_left = sign(d2fwdS2_discleft*dotvn + d2gamadS2_discleft*dotvg);
    sign2der_right = sign(d2fwdS2_discright*dotvn + d2gamadS2_discright*dotvg);
    
    %Define the Rankine-Hugoniout velocity
    charvel_rh = dotvn*dfwdS_rh + dotvg*dgamadS_rh;
    
    %----------------------------------------------------------------------
    %Choise according second derivative sign (see Serma, 2009)
    %Get the max value of characteristic velocity
    %Define a range for the saturtion
    
    Sranglr = [Sleft Sright];
   
    %Get the analitical derivative:
    
    [dfwdS,dgamadS] = calcdfunctiondS(0,0,Sranglr,1);
   
    
   
    [numflux, earlysw]=riemannsolvertwophaseflow(signder_left,signder_right,sign2der_left,...
    sign2der_right,Sright,Sleft,bedgesize, inedg,fw,dotvn,dotvg,...
     gama,charvel_rh,dfwdS,dfwdS_rh,dgamadS,mLLF,limiterflag);
    
    %Obtain the contribution of interface over element to LEFT
    advecterm(leftelem) = advecterm(leftelem) + numflux;
    %Obtain the contribution of interface over element to RIGHT
    advecterm(rightelem) = advecterm(rightelem) - numflux;
    
end  %End of FOR ("inedge")

