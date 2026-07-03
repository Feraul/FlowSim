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

function [advecterm] = hyperb_numflux(Sw,flowrate,taylorterms,limiterflag,...
    flagknownvert,satonvertices,satonboundedges,pointbndedg,pointinedg,...
    orderbedgdist,orderinedgdist,constraint,mlplimiter)
%Define global parameters:
global elem bedge inedge numcase centelem coord;
%Initialize a tolerance. It is a computational zero
tol = 1e-9;
countinter=0;
%Initialize "advecterm" and "bodyterm"
advecterm = zeros(size(elem,1),1);
%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);

%--------------------------------------------------------------------------
%Boundary edges (when it exists):

%Chose the strategy (in "bedge") according to benchmark number.
%Linear Advection obtained from Goosh and Van Altena (2002) or
%Wang and Liu (2004)
if numcase == 101 || numcase == 103 || numcase==107 || numcase==108
    %Get the periodic elements.
    periodicpos = getperiodicelem;
    %Calculate the flux:
    %Swept "bedge"
    for i = 1:length(pointbndedg)
        
        %Initialize some parameters:
        ibedg = pointbndedg(i);
        
        %Define the "vertices"
        vertices = bedge(ibedg,1:2);
        %Get the coordinate for the vertices:
        verticescoord = coord(vertices,:);
        %Define left element
        leftelem = bedge(ibedg,3);
        
        if bedge(ibedg,5)<600
            
            %Left contribution:
            %Define the order for this edge.
            faceorder = orderbedgdist(i,1);
            
            elemeval=leftelem;
            
            %Get the saturation value recovered
            [Sleftlim,Sleftnonlim] = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,limiterflag,...
                faceorder,constraint,flagknownvert,satonvertices,mlplimiter,centelem(elemeval,1:2));
            
            %Define the normal velocity into face
            if numcase==103
                if Sleftlim>1
                    %Sleftlim=Sw(leftelem);
                    Sleftlim=1;
                elseif Sleftlim<-1
                    %Sleftlim=Sw(leftelem);
                    Sleftlim=-1;
                end
                
            end
            Sleft=Sleftlim;
            dotvn = flowrate(ibedg);
            
            if numcase==107 
                %UPWIND
                %Calculate Numerical Flux:               
                %Calculate the numerical flux through interface
                numflux = (Sleft^2/2)*dotvn;
                
            else
                %Calculate Numerical Flux:
                numflux = max(dotvn,0)*Sleft;
            end
            %Obtain the contribution of interface over element to LEFT
            advecterm(leftelem) = advecterm(leftelem) + numflux;
            
        else
            
            %Define the vertices for the pseudo right element
            psdvertices = bedge(periodicpos(ibedg),1:2);
            %Define a pseudo right element
            psdrightelem = bedge(periodicpos(ibedg),3);
            
            %Left contribution:
            %Define the order for this edge.
            faceorder = orderbedgdist(i,1);
            %Define the elements that share the edge evaluated
            elemeval = [leftelem psdrightelem];
            
            
            %Get the saturation value recovered
            [Sleftlim,Sleftnonlim] = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,limiterflag,...
                faceorder,constraint,flagknownvert,satonvertices,mlplimiter,centelem(elemeval,1:2));
            
            %Periodic contribution:
            %Define the order for this edge.
            faceorder = orderbedgdist(i,2);
            %Define the elements that share the edge evaluated
            elemeval = [psdrightelem leftelem];
            faceaux=find(bedge(:,3)==psdrightelem);
            if numcase==107
                
                j=1;
                for qqq=1:length(faceaux)
                    mmm=faceaux(qqq);
                    if bedge(mmm,5)>600
                        faceaux1(j)=mmm;
                        j=j+1;
                    end
                end
            else
                
                faceaux1=faceaux;
            end
            %Define the "vertices"
            verticesaux = bedge(faceaux1,1:2);
            %Get the coordinate for the vertices:
            verticescoord = coord(verticesaux,:);
            
            
            %Get the saturation value recovered
            [Srightlim,Srightnonlim] = getsatonedge(elemeval,psdvertices,verticescoord,taylorterms,Sw,...
                limiterflag,faceorder,constraint,flagknownvert,satonvertices,...
                mlplimiter,centelem(elemeval,1:2));
            %Define the normal velocity into face
            if numcase==103
                if Sleftlim>1
                    %Sleftlim=Sw(leftelem);
                    Sleftlim=1;
                elseif Sleftlim<-1
                    %Sleftlim=Sw(leftelem);
                    Sleftlim=-1;
                end
                if Srightlim>1
                    %Srightlim=Sw(rightelem);
                    Srightlim=1;
                elseif Srightlim<-1
                    %Srightlim=Sw(rightelem);
                    Srightlim=-1;
                end
            end
            dotvn = flowrate(ibedg);
            if (strcmp(limiterflag{8},'on')|| strcmp(limiterflag{10},'on')|| strcmp(limiterflag{11},'on') ) && (countinter==0)
                [Sleft,Sright,mLLF]=PhysicalAD(Sw,taylorterms,limiterflag,flagknownvert,satonvertices,...
                    constraint,mlplimiter,leftelem,psdrightelem,vertices,...
                    verticescoord,Sleftlim,Sleftnonlim,Srightlim,Srightnonlim,dotvn,tol,0,0);
            else
                Sleft= Sleftlim;
                Sright=Srightlim;
            end
            if numcase==107 || numcase==108
                
                alfamax = max(abs(0.5*(Sright + Sleft)*dotvn ));
                %alfamax = max(abs(dfwdS*dotvn + dgamadS*dotvg ));
                %Denine the numerical flux
                Fleft = (Sleft^2/2)*dotvn ;
                Fright = (Sright^2/2)*dotvn ;
                %Define Local Lax-Friedrichs Flux
                
                numflux = 0.5*((Fleft + Fright) - alfamax*(Sright - Sleft));
            else
                %Calculate Numerical Flux:
                numflux = max(dotvn,0)*Sleft + min(dotvn,0)*Sright;
            end
            %Obtain the contribution of interface over element to LEFT
            advecterm(leftelem) = advecterm(leftelem) + numflux;
            %Obtain the contribution of interface over element to RIGHT
            advecterm(psdrightelem) = advecterm(psdrightelem) - numflux;
        end
        
        
    end  %End of FOR (Swept "bedge" untill the half)
    
    %Sonar (1994) and Zalezak (1979). These problems are a rotating profile.
else
    %Swept "bedge"
    for i = 1:length(pointbndedg)
        %Initialize some parameters:
        ibedg = pointbndedg(i);
        
        %Get the element on the left
        leftelem = bedge(ibedg,3);
        %There is a prescribed saturation
        %Attribute the saturation on boundary
        Sleft = satonboundedges(ibedg);
        
        %Define the normal velocity into face
        dotvn = flowrate(ibedg);
        
        %Calculate the numerical flux through interface.
        numflux = dotvn*Sleft;
        %Obtain the contribution of interface over element to LEFT
        advecterm(leftelem) = advecterm(leftelem) + numflux;
    end  %End of FOR (Swept "bedge")
end  %End of IF

%--------------------------------------------------------------------------
%Internal edges:

%Swept "inedge" evaluating left and right elements by edge. Apply
%approximated Riemann Solver through edge.
for i = 1:length(pointinedg)
 
    %Initialize some parameters:
    inedg = pointinedg(i);
    
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
    %Get the saturation value recovered on each quadrature point ("on_q")
    %Get the statment regarding to mlp limiter:
    boolmlp = (length(mlplimiter) == 1);
    mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
    %Get the saturation value recovered
    
    [Sleftlim,Sleftnonlim] = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,limiterflag,...
        faceorder,constraint,flagknownvert,satonvertices,mlpbyelem,centelem(elemeval,1:2));
    
    
    %Right Contribution:
    %Define the order for this edge.
    faceorder = orderinedgdist(i,2);
    %Define the elements that share the edge evaluated
    elemeval = [rightelem leftelem];
    %Get the saturation value recovered
    %Define the normal velocity in each face
    dotvn = flowrate(bedgesize + inedg);
    %Get the statment regarding to mlp limiter:
    mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
    [Srightlim,Srightnonlim] = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,limiterflag,...
        faceorder,constraint,flagknownvert,satonvertices,mlpbyelem,centelem(elemeval,1:2));
    if numcase==103
        if Sleftlim>1
            Sleftlim=1;
        elseif Sleftlim<-1
            Sleftlim=-1;
        end
        if Srightlim>1
            Srightlim=1;
        elseif Srightlim<-1
            Srightlim=-1;
        end
    end
    if (strcmp(limiterflag{8},'on')|| strcmp(limiterflag{10},'on')|| strcmp(limiterflag{11},'on') ) && (countinter==0)
        [Sleft,Sright,mLLF]=PhysicalAD(Sw,taylorterms,limiterflag,flagknownvert,satonvertices,...
            constraint,mlplimiter,leftelem,rightelem,vertices,...
            verticescoord,Sleftlim,Sleftnonlim,Srightlim,Srightnonlim,dotvn,tol, rijL,rijR);
        
    else
        Sleft=Sleftlim;
        Sright=Srightlim;
    end
    if numcase==107 || numcase==108
                
        alfamax = max(abs(0.5*(Sright + Sleft)*dotvn ));
        %Denine the numerical flux
        Fleft = (Sleft^2/2)*dotvn ;
        Fright = (Sright^2/2)*dotvn ;
        %Define Local Lax-Friedrichs Flux
        
        numflux = 0.5*((Fleft + Fright) - alfamax*(Sright - Sleft));
        
    else
        
        % UPWIND
        %Calculate Numerical Flux:
        numflux = max(dotvn,0)*Sleft + min(dotvn,0)*Sright;
        
    end
    %Obtain the contribution of interface over element to LEFT
    advecterm(leftelem) = advecterm(leftelem) + numflux;
    %Obtain the contribution of interface over element to RIGHT
    advecterm(rightelem) = advecterm(rightelem) - numflux;
    
    
end  %End of FOR ("inedge")



