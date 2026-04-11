%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with
%high-order resolution
%Type of file: FUNCTION
%Programer: Márcio Souza
%Modified: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals:
%.

%--------------------------------------------------------------------------
%Additional comments:
%

%--------------------------------------------------------------------------

function [advecterm,entrineqterm,earlysw,vectorSleft,vectorSright] = calcnumflux(Con,Fg,flowrateadvec,flowratedif,...
    taylorterms,limiterflag,flagknownvert,satonvertices,flagknownedge,...
    satonboundedges,pointbndedg,pointinedg,orderbedgdist,orderinedgdist,...
    constraint,mlplimiter,earlysw,countinter)
%Define global parameters:
global coord elem bedge inedge normals dens numcase centelem order courant;

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
        
        %Define the elements that share the edge evaluated
        elemeval = leftelem;
        
        %Verify if there is concentration prescribed on boundary:
        %There is a prescribed saturation
        if flagknownedge(ibedg) == 1
            %Attribute the saturation on boundary
            Sleft = satonboundedges(ibedg);
            
            %There is no prescribed saturation. It is necessary calculate.
        else
            %Get the statment regarding to mlp limiter:
            boolmlp = (length(mlplimiter) == 1);
            mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
            
            %Get the saturation value recovered
            [Sleftlim,Sleftnonlim] = getsatonedge(elemeval,vertices,verticescoord,...
                taylorterms,Con,limiterflag,faceorder,constraint,flagknownvert,...
                satonvertices,mlpbyelem,centelem(leftelem,1:2));
            
            if strcmp(limiterflag{9},'on')|| strcmp(limiterflag{11},'on')|| strcmp(limiterflag{12},'on')
                
                Sleft= Sleftnonlim;
                % %==========================================================================
                
            else
                Sleft=Sleftlim;
            end
            
            
        end  %End of IF
        
        %Fill "earlysw"
        earlysw(ibedg) = Sleft;
        
        %Define the normal velocity into face
        dotvn = flowrateadvec(ibedg);
        %Get accuracy for "dotvn"
        dotvn = dotvn*(abs(dotvn) > tol);
        %
        dotdif=flowratedif(ibedg);
        %Calculate the numerical flux through interface.
        numflux = dotvn*Sleft + dotdif;
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
    
    %Initialize some parameters:
    inedg = pointinedg(i);
    %---------------------------------------
    %Define the normal velocity in each face
    dotvn = flowrateadvec(bedgesize + inedg);
    dotdif= flowratedif(bedgesize + inedg);
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
    if size(Fg,2) > 1
        dotvg = Fg(bedgesize + inedg,1)*(dens(1) - dens(2))/lenginjface;
        %There is NO gravity
    else
        dotvg = 0;
    end  %End of IF
    %----------------------------------------------------------------------
    
    %Get the saturation value recovered on each quadrature point ("on_q")
    %Get the statment regarding to mlp limiter:
    boolmlp = (length(mlplimiter) == 1);
    mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
    
    %Left Contribution:
    [Sleftlim,Sleftnonlim] = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Con,...
        limiterflag,faceorder,constraint,flagknownvert,satonvertices,...
        mlpbyelem,centelem(elemeval,1:2));
    
    %Right Contribution:
    %Define the order for this edge.
    faceorder = orderinedgdist(i,2);
    %Define the elements that share the edge evaluated
    elemeval = [rightelem leftelem];
    %Get the statment regarding to mlp limiter:
    mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
    
    %Get the saturation value recovered on each quadrature point ("on_q")
    [Srightlim,Srightnonlim]= getsatonedge(elemeval,vertices,verticescoord,taylorterms,Con,...
        limiterflag,faceorder,constraint,flagknownvert,satonvertices,...
        mlpbyelem,centelem(elemeval,1:2));
    
    %PAD
    if (strcmp(limiterflag{9},'on')|| strcmp(limiterflag{11},'on')|| strcmp(limiterflag{12},'on') ) && (countinter==0)
        [Sleft,Sright,mLLF]=PhysicalAD(Con,taylorterms,limiterflag,flagknownvert,satonvertices,...
            constraint,mlplimiter,leftelem,rightelem,vertices,...
            verticescoord,Sleftlim,Sleftnonlim,Srightlim,Srightnonlim,dotvn,tol,rijL,rijR);
    else
        Sleft= Sleftlim;
        Sright=Srightlim;
    end
    %%
    
    vectorSleft(size(bedge,1)+i,1)=Sleft;
    vectorSright(i,1)=Sright;
    %Discrete:
    
    %Get accuracy:
    dotvn = dotvn*(abs(dotvn) > tol);
    dotvg = dotvg*(abs(dotvg) > tol);
    %---------------------------------------
    
    
    %Define the Rankine-Hugoniout velocity
    charvel_rh = dotvn;
    
    method=limiterflag{14};
    
    [numflux, Saproxima]=riemannsolver(Sright,Sleft,method,bedgesize, inedg,dotvn,dotvg,charvel_rh,dotdif);
    earlysw(bedgesize + inedg) = Saproxima;
    %Obtain the contribution of interface over element to LEFT
    advecterm(leftelem) = advecterm(leftelem) + numflux;
    %Obtain the contribution of interface over element to RIGHT
    advecterm(rightelem) = advecterm(rightelem) - numflux;
    
end  %End of FOR ("inedge")

