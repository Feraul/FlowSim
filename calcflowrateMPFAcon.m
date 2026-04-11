%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 29/06/2013 (mother's visit)
%Modify data:   /  /2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
% Modified by Fernando Contreras 2022
%--------------------------------------------------------------------------
%Goals:
%   

%--------------------------------------------------------------------------
%Additional comments:


%--------------------------------------------------------------------------

function [flowrate,flowresult,flowratedisp] = calcflowrateMPFAcon(pressure,...
      Con, transmvecleftcon,transmvecrightcon,knownvecleftcon,...
      knownvecrightcon,storeinvcon,Bleftcon,Brightcon,Fgcon,mapinvcon,...
      maptransmcon,mapknownveccon,pointedgecon, bodytermcon,...
      transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
            Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge, bodyterm,mobility)
%Define global parameters:
global elem bedge inedge phasekey smethod;

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize "counters"
countmtrx = 0;
countvec = 1;

%Initialise "flowresult". This accumulates the flowrate in all elements 
%associated to producer wells. It's a flow rate resulting for each element. 
flowresult = zeros(size(elem,1),1);
flowresultdisp = zeros(size(elem,1),1);

%--------------------------------------------------------------------------
%Switch the way to calculate "flowrate" according to "phasekey" and
%"satkey" values

%--------------------------------------------------------------------------
%One-phase aplications or Two-phase aplications with non-waveoriented
%schemes (calculate the "flowrate" on whole edges).

if phasekey == 1 || (phasekey == 3 && strcmp(smethod,'stdr')) || ...
        (phasekey == 3 && strcmp(smethod,'goef')) 
    %Initialize "flowrate"
    flowrate = zeros(bedgesize + inedgesize,1);
    flowratedisp = zeros(bedgesize + inedgesize,1);

    %Swept the boundary edges
    for i = 1:bedgesize
        %Define the flow through edge:
        %Get the vertices.
        vertices = bedge(i,1:2);
        %Define the element on the left
        elemleft = bedge(i,3);
    
        %------------------------------------------------------------------
        %Vertex 1:
    
        %Get the elements and the nodes surrounding the vertex 1
        [esurn,] = getsurnode(vertices(1));
    
        %Get "mtxleft" and "vecleft"
        %Get the transmissibility coefficients:
        %Get the transmissibility coefficients:
        [mtxleft,vecleft,] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(1),vertices(2),phasekey,i,1,mobility,bedgesize,...
        inedgesize);
    
        %Calculate the flowrate for edge evaluated
        getflowrate = mtxleft'*pressure(esurn) + vecleft;
        %Attribute to "flowrate" the darcy velocity
        flowrate(i) = flowrate(i) + getflowrate;  
        
        %Attribute to "flowresult" the darcy velocity
        flowresult(elemleft) = flowresult(elemleft) + getflowrate;  
        
        %Update the counters:
        countmtrx = countmtrx + length(esurn);
        countvec = countvec + 1;
        
        %------------------------------------------------------------------
        %Vertex 2:
    
        %Get the elements and the nodes surrounding the vertex 1
        [esurn,] = getsurnode(vertices(2));
    
        %Get the transmissibility coefficients:
        [mtxleft,vecleft,] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(2),vertices(1),phasekey,i,2,mobility,bedgesize,...
        inedgesize);
        
        %Calculate the flowrate for edge evaluated
        getflowrate = mtxleft'*pressure(esurn) + vecleft;
        %Attribute to "flowrate" the darcy velocity
        flowrate(i) = flowrate(i) + getflowrate;  
        
        %Attribute to "flowresult" the darcy velocity
        flowresult(elemleft) = flowresult(elemleft) + getflowrate;  
        
        %------------------------------------------------------------------
        
        %Update the counters:
        countmtrx = countmtrx + length(esurn);
        countvec = countvec + 1;
        %% ==================================================================
        % calculate the flowrate dispersive
        %------------------------------------------------------------------
        %Vertex 1:
    
        %Get the elements and the nodes surrounding the vertex 1
        [esurn,] = getsurnode(vertices(1));
    
        %Get "mtxleft" and "vecleft"
        %Get the transmissibility coefficients:
        %Get the transmissibility coefficients:
        [mtxleftcon,vecleftcon,] = gettransmcoeff(transmvecleftcon,knownvecleftcon,...
        storeinvcon,Bleftcon,1,mapinvcon,maptransmcon,mapknownveccon,pointedgecon,...
        bodytermcon,vertices(1),vertices(2),phasekey,i,1,mobility,bedgesize,...
        inedgesize);
    
        %Calculate the flowrate for edge evaluated
        getflowratecon = mtxleftcon'*Con(esurn) + vecleftcon;
        %Attribute to "flowrate" the darcy velocity
        flowratedisp(i) = flowratedisp(i) + getflowratecon;  
        
        %Attribute to "flowresult" the darcy velocity
        flowresultdisp(elemleft) = flowresultdisp(elemleft) + getflowratecon;  
        
        
        
        %------------------------------------------------------------------
        %Vertex 2:
    
        %Get the elements and the nodes surrounding the vertex 1
        [esurn,] = getsurnode(vertices(2));
    
        %Get the transmissibility coefficients:
        [mtxleftcon,vecleftcon,] = gettransmcoeff(transmvecleftcon,knownvecleftcon,...
        storeinvcon,Bleftcon,1,mapinvcon,maptransmcon,mapknownveccon,pointedgecon,...
        bodytermcon,vertices(2),vertices(1),phasekey,i,2,mobility,bedgesize,...
        inedgesize);
        
        %Calculate the flowrate for edge evaluated
        getflowratecon = mtxleftcon'*Con(esurn) + vecleftcon;
        %Attribute to "flowrate" the darcy velocity
        flowratedisp(i) = flowratedisp(i) + getflowratecon;  
        
        %Attribute to "flowresult" the darcy velocity
        flowresultdisp(elemleft) = flowresultdisp(elemleft) + getflowratecon;  
        
        
        
    end  %End of FOR (Swept boundary edges)
    
    %Swept the internal edges
    for i = 1:inedgesize
        %Define the flow through edge:
        %Get the vertices.
        vertices = inedge(i,1:2);
        %Define the elements on the left and on the right
        elemleft = inedge(i,3);
        elemright = inedge(i,4);
        
        %------------------------------------------------------------------
        %Vertex 1:
    
        %Get the elements and the nodes surrounding the vertex 1
        [esurn,] = getsurnode(vertices(1));
    
        %Get the transmissibility coefficients (left contribution):
        [mtxleft,vecleft,] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(1),vertices(2),phasekey,bedgesize + i,1,...
        mobility,bedgesize,inedgesize);
        
        %Calculate the flowrate for edge evaluated
        getflowrate = mtxleft'*pressure(esurn) + vecleft;
        %Attribute to "flowrate" the darcy velocity
        flowrate(i + bedgesize) = flowrate(i + bedgesize) + getflowrate;  
        
        %Attribute to "flowresult" the darcy velocity (on the left)
        flowresult(elemleft) = flowresult(elemleft) + getflowrate;  
        %Attribute to "flowresult" the darcy velocity (on the right)
        flowresult(elemright) = flowresult(elemright) - getflowrate;  
        
        %Update the counters:
        countmtrx = countmtrx + length(esurn);
        countvec = countvec + 1;
        
        %------------------------------------------------------------------
        %Vertex 2:
    
        %Get the elements and the nodes surrounding the vertex 1
        [esurn,] = getsurnode(vertices(2));
    
        %Get the transmissibility coefficients (left contribution):
        [mtxleft,vecleft,] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(2),vertices(1),phasekey,bedgesize + i,2,...
        mobility,bedgesize,inedgesize);
        
        %Calculate the flowrate for edge evaluated
        getflowrate = mtxleft'*pressure(esurn) + vecleft;
        %Attribute to "flowrate" the darcy velocity
        flowrate(i + bedgesize) = flowrate(i + bedgesize) + getflowrate;  
        
        %Attribute to "flowresult" the darcy velocity (on the left)
        flowresult(elemleft) = flowresult(elemleft) + getflowrate;  
        %Attribute to "flowresult" the darcy velocity (on the right)
        flowresult(elemright) = flowresult(elemright) - getflowrate;  
        
        %Update the counters:
        countmtrx = countmtrx + length(esurn);
        countvec = countvec + 1;
        %% =================================================================
         % calculate the flowrate dispersive
         %------------------------------------------------------------------
        %Vertex 1:
    
        %Get the elements and the nodes surrounding the vertex 1
        [esurn,] = getsurnode(vertices(1));
    
        %Get the transmissibility coefficients (left contribution):
        [mtxleftcon,vecleftcon,] = gettransmcoeff(transmvecleftcon,knownvecleftcon,...
        storeinvcon,Bleftcon,1,mapinvcon,maptransmcon,mapknownveccon,pointedgecon,...
        bodytermcon,vertices(1),vertices(2),phasekey,bedgesize + i,1,...
        mobility,bedgesize,inedgesize);
        
        %Calculate the flowrate for edge evaluated
        getflowratecon = mtxleftcon'*Con(esurn) + vecleftcon;
        %Attribute to "flowrate" the darcy velocity
        flowratedisp(i + bedgesize) = flowratedisp(i + bedgesize) + getflowratecon;  
        
        %Attribute to "flowresult" the darcy velocity (on the left)
        flowresultdisp(elemleft) = flowresultdisp(elemleft) + getflowratecon;  
        %Attribute to "flowresult" the darcy velocity (on the right)
        flowresultdisp(elemright) = flowresultdisp(elemright) - getflowratecon;  
        
        %------------------------------------------------------------------
        %Vertex 2:
    
        %Get the elements and the nodes surrounding the vertex 1
        [esurn,] = getsurnode(vertices(2));
    
        %Get the transmissibility coefficients (left contribution):
        [mtxleftcon,vecleftcon,] = gettransmcoeff(transmvecleftcon,knownvecleftcon,...
        storeinvcon,Bleftcon,1,mapinvcon,maptransmcon,mapknownveccon,pointedgecon,...
        bodytermcon,vertices(2),vertices(1),phasekey,bedgesize + i,2,...
        mobility,bedgesize,inedgesize);
        
        %Calculate the flowrate for edge evaluated
        getflowratecon = mtxleftcon'*Con(esurn) + vecleftcon;
        %Attribute to "flowrate" the darcy velocity
        flowratedisp(i + bedgesize) = flowratedisp(i + bedgesize) + getflowratecon;  
        
        %Attribute to "flowresult" the darcy velocity (on the left)
        flowresultdisp(elemleft) = flowresultdisp(elemleft) + getflowratecon;  
        %Attribute to "flowresult" the darcy velocity (on the right)
        flowresultdisp(elemright) = flowresultdisp(elemright) - getflowratecon;  
        
    end  %End of FOR ("inedge")
    
%--------------------------------------------------------------------------
%Two-phase aplications with Multidimensional scheme or with Spectral Finite 
%Volume (calculate the "flowrate" on half-edges).

elseif (phasekey == 2 && strcmp(smethod,'mwic')) || ...
        (phasekey == 2 && strcmp(smethod,'mwec')) || ...
        (phasekey == 2 && strcmp(smethod,'rtmd'))
    %Initialize "flowrateinhalfedge"
    flowrate = zeros(2*(bedgesize + inedgesize),1);
    
    %Initialize auxiliary counter
    m = 1;
    %Swept the boundary edges
    for i = 1:bedgesize
        %Define the flow through edge:
        %Get the vertices.
        vertices = bedge(i,1:2);
        %Define the element on the left
        elemleft = bedge(i,3);
        
        %------------------------------------------------------------------
        %Vertex 1:
    
        %Get the elements and the nodes surrounding the vertex 1
        [esurn,] = getsurnode(vertices(1));
    
        %Get the transmissibility coefficients:
        [mtxleft,vecleft,] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(1),vertices(2),phasekey,i,1,mobility,bedgesize,...
        inedgesize);
        
        %Calculate the flowrate for edge evaluated
        getflowrate = mtxleft'*pressure(esurn) + vecleft;
        %Attribute to "flowrate" the darcy velocity
        flowrate(m) = getflowrate;  
        
        %Attribute to "flowresult" the darcy velocity
        flowresult(elemleft) = flowresult(elemleft) + getflowrate;  
        
        %Update the counters:
        countmtrx = countmtrx + length(esurn);
        countvec = countvec + 1;
        
        %Increment "m"
        m = m + 1;
    
        %------------------------------------------------------------------
        %Vertex 2:
    
        %Get the elements and the nodes surrounding the vertex 1
        [esurn,] = getsurnode(vertices(2));
    
        %Get the transmissibility coefficients:
        [mtxleft,vecleft,] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(2),vertices(1),phasekey,i,2,mobility,bedgesize,...
        inedgesize);
        
        %Calculate the flowrate for edge evaluated
        getflowrate = mtxleft'*pressure(esurn) + vecleft;
        %Attribute to "flowrate" the darcy velocity
        flowrate(m) = getflowrate;  
        
        %Attribute to "flowresult" the darcy velocity
        flowresult(elemleft) = flowresult(elemleft) + getflowrate;  
        
        %Update the counters:
        countmtrx = countmtrx + length(esurn);
        countvec = countvec + 1;

        %Increment "m"
        m = m + 1;
    end  %End of FOR (Swept boundary edges)
    
    %Swept the internal edges
    for i = 1:inedgesize
        %Define the flow through edge:
        %Get the vertices.
        vertices = inedge(i,1:2);
        %Define the elements on the left and on the right
        elemleft = inedge(i,3);
        elemright = inedge(i,4);
    
        %------------------------------------------------------------------
        %Vertex 1:
    
        %Get the elements and the nodes surrounding the vertex 1
        [esurn,] = getsurnode(vertices(1));
    
        
       %Get the transmissibility coefficients (left contribution):
       [mtxleft,vecleft,] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(1),vertices(2),phasekey,bedgesize + i,1,...
        mobility,bedgesize,inedgesize);
        
        %Calculate the flowrate for edge evaluated
        getflowrate = mtxleft'*pressure(esurn) + vecleft;
        %Attribute to "flowrate" the darcy velocity
        flowrate(m) = getflowrate;  
        
        %Attribute to "flowresult" the darcy velocity (on the left)
        flowresult(elemleft) = flowresult(elemleft) + getflowrate;  
        %Attribute to "flowresult" the darcy velocity (on the right)
        flowresult(elemright) = flowresult(elemright) - getflowrate;  
        
        %Update the counters:
        countmtrx = countmtrx + length(esurn);
        countvec = countvec + 1;

        %Increment "m"
        m = m + 1;

        %------------------------------------------------------------------
        %Vertex 2:
    
        %Get the elements and the nodes surrounding the vertex 1
        [esurn,] = getsurnode(vertices(2));
    
        %Get the transmissibility coefficients (left contribution):
        [mtxleft,vecleft,] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(2),vertices(1),phasekey,bedgesize + i,2,...
        mobility,bedgesize,inedgesize);
        
        %Calculate the flowrate for edge evaluated
        getflowrate = mtxleft'*pressure(esurn) + vecleft;
        %Attribute to "flowrate" the darcy velocity
        flowrate(m) = getflowrate;  
        
        %Attribute to "flowresult" the darcy velocity (on the left)
        flowresult(elemleft) = flowresult(elemleft) + getflowrate;  
        %Attribute to "flowresult" the darcy velocity (on the right)
        flowresult(elemright) = flowresult(elemright) - getflowrate;  
        
        %Update the counters:
        countmtrx = countmtrx + length(esurn);
        countvec = countvec + 1;

        %Increment "m"
        m = m + 1;
    end  %End of FOR ("inedge")
end  %End of IF (calculate "flowrate" in half-edge or whole edge)
