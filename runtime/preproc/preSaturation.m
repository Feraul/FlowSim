%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 04/12/2013
%Modify data: 02/07/2015
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%Modified: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%Additional Comments:

%--------------------------------------------------------------------------

function [wvector,wmap,constraint,massweigmap,othervertexmap,lsw,...
    swsequence,ntriang,areatriang,prodwellbedg,prodwellinedg,mwmaprodelem,...
    vtxmaprodelem,coordmaprodelem,amountofneigvec,rtmd_storepos,...
    rtmd_storeleft,rtmd_storeright,isonbound] = ...
    preSaturation(flagknownedge,injecelem,producelem)
%Define global parameters:
global order recovtype numcase;

%Message to user:
disp(' ');
disp('---------------------------------------------------');
disp('>> Preprocessing Saturation Equation...');

%Initialize the parameters:
wvector = 0;
wmap = 0;
constraint = 0;
lsw = 0;
swsequence = 0;
ntriang = 0;
areatriang = 0;
%Wells edges (Standard schemes)
prodwellbedg = 0;
prodwellinedg = 0;
%Wells edges (MultiD schemes)
mwmaprodelem = 0;
vtxmaprodelem = 0;
coordmaprodelem = 0;
amountofneigvec = 0;
%MultiD parameters:
massweigmap = 0;
othervertexmap = 0;
%Realy Truly MultiDimensional Scheme
rtmd_storepos = 0;
rtmd_storeleft = 0;
rtmd_storeright = 0;
isonbound = 0;

%--------------------------------------------------------------------------
%Chose the preprocessment according to "smethod" and "order"

%Preprocess the edges shared by producer wells. It is common to all
%schemes.
if any(producelem) && numcase<200
    %Get the edges shared by producer wells
    [prodwellbedg,prodwellinedg] = getedgesonprodwell(producelem);
end  %End of IF

%--------------------------------------------------------------------------
%Preprocessing HIGHER ORDER terms:

%Higher Order (MUSCL-Durlofsky or MUSCL-Least Square)    
if order > 1
    %Verify if the strategy of higher order is a Least Squares or that
    %one proposed by Durlofsky (1993)
    %Durlofsky (1993)
    if strcmp(recovtype,'ggg')
        %Get preprocessing terms of second order Durlofisky scheme
        [swsequence,ntriang,areatriang] = ...
            getgggradterms(flagknownedge);
        
    %Least Squares
    else
        %Precalculate the "A" matrix and the weights ("w") for FACE 
        %and FULL neighboring.
        [wvector,wmap,constraint,lsw,amountofneigvec] = ...
            getleastsquarescoeff(flagknownedge,injecelem);
    end  %End of IF
end  %End of IF
    

%Message to user:
disp('>> "preSaturation" was finished with success!');


%Get the coefficients "^". See Goosh et al.(2007)
function[coefficient] = getcoefnuminteg(elemevalcoord,neighcent,vertcoord,...
    neighbour,n,m)
%Initialize parameters:
coefx = zeros(n + 1,1);
secsum = zeros(m + 1,1);
    
%Produce the numerical integration (see Goosh et al., 2007):
%Swept "y"
for l = 0:m
    %Calculate the coefficients of the first summation (Eq. 8 in Goosh 
    %et al., 2007).
    coefy = (calcfact(m)/(calcfact(l)*calcfact(m - l)))*...
        ((neighcent(2) - elemevalcoord(2))^l);
    %Swept "x"
    for k = 0:n
        %Calculate the numerical integral OF difference (x - xj) over 
        %the control volume "j" (neighbour).
        integdiffer = getintegdiff(neighbour,vertcoord,...
            neighcent,n - k,m - l);
        %Calculate the coefficients of the second summation (Eq. 8 in 
        %Goosh et al., 2007).
        coefx(k + 1) = (calcfact(n)/(calcfact(k)*calcfact(n - k)))*...
            ((neighcent(1) - elemevalcoord(1))^k)*integdiffer;
    end  %End of FOR ("x")
    %Calculate the second summation.
    secsum(l + 1) = coefy*sum(coefx);
end  %End of FOR ("y")

%Calculate the first summation
coefficient = sum(secsum);

%--------------------------------------------------------------------------
%Function "attriboundcontrib"
%--------------------------------------------------------------------------

function [A,w] = attriboundcontrib(Ain,win,order,ielem,elemevalcoord,...
    consterms,bedgrownum,flagknownedge)
%Define global variables
global coord elem bedge;

%Define a tolerance (computational zero)
tol = 1e-12;
%Update "A" and "w"
A = Ain;
w = win;

%Verify if the element has boundary face
pointelembound = logical(bedge(:,3) == ielem);
%Verify if there exists any boundary in the element evaluated. 
if any(pointelembound)
    %Get the number of the row ("bedge") which has the element evaluated.
    getbedgrownum = bedgrownum(pointelembound); 
    
    %"getedgesinbound" stores the number of vert. which define each edge.
    getedgesinbound = bedge(pointelembound,1:2);
    %Get "getedgesinbound" size
    edgesinboundsize = size(getedgesinbound,1);

    %Initialize other variables
    Aaux = zeros(1,2);
    neighborcoordaux = Aaux;
    waux = zeros(1,edgesinboundsize);
    
    %Initialize "icount"
    icount = 1;
    %Swept all edges in boundary (for the element evaluated)
    for i = 1:edgesinboundsize
        %Verify if it is a triangle ("elemtype" = 3) or a quaddrangle 
        %("elemtype" = 4)
        elemtype = sum(elem(ielem,1:4) ~= 0); 
        %Get the vertices of the neighbour element.
        vertices = elem(ielem,1:elemtype);
        %Get the vertices coordinate (of the element which we want mirror)
        vertcoord = coord(vertices,:);

        %Get the vertices coordinate:
        firstnodecoord = coord(getedgesinbound(i,1),:);
        secondnodecoord = coord(getedgesinbound(i,2),:);
        
        %Evaluate if the boundary condition is a known saturation on the 
        %edge. It occurs in the Buckley-Leverett application.

        %-------------------------------------
        %It is a known saturation on the edge:
        if flagknownedge(getbedgrownum(i)) == 1
            %End-points coordinate:
            a = firstnodecoord(1:2);
            b = secondnodecoord(1:2);
            
            %For SECOND ORDER (edge with known saturation value)
            if order == 2
                %Get the coordinate of the midpoint (on the edge)
                midedgecoord = (a + b)/2;
                %Get the weight (weighted Least Squares). 
                %This scalar is called how many quadrature points exist.
                wquadpt = calcgeometricweight(elemevalcoord,...
                    midedgecoord);
                %Define Component on "x" and on "y"
                delta = midedgecoord - elemevalcoord(1:2);
                dx = delta(1);
                dy = delta(2);
                %Attribute it to "Aaux"
                Aaux(icount,1:2) = wquadpt*[dx dy];
                
                %Fill "waux"
                waux(icount) = wquadpt;
                %Update "icount" (one quadrature point)
                icount = icount + 1;

            %For THIRD or FOURTH ORDER (edge with known saturation value)
            elseif order == 3 || order == 4 
                %Get the Gauss point coordinates:
                %Point nearest the first vertex:
                c1 = ((a + b)/2) - ((b - a)/(2*sqrt(3)));
                %Point nearest the second vertex:
                c2 = ((a + b)/2) + ((b - a)/(2*sqrt(3)));
                
                %Quadrature point 1:
                
                %Get the weight (weighted Least Squares). 
                %This scalar is called how many quadrature points exist.
                wquadpt1 = calcgeometricweight(elemevalcoord,c1);
                %Define Component on "x" and on "y"
                delta1 = c1 - elemevalcoord(1:2);
                dx1 = delta1(1);
                dy1 = delta1(2);
                
                %Quadrature point 2:
                
                %Get the weight (weighted Least Squares). 
                %This scalar is called how many quadrature points exist.
                wquadpt2 = calcgeometricweight(elemevalcoord,c2);
                %Define Component on "x" and on "y"
                delta2 = c2 - elemevalcoord(1:2);
                dx2 = delta2(1);
                dy2 = delta2(2);
                
                %3rd Order
                if order == 3
                    %Attribute it to "Aaux" (3rd order)
                    Aaux(icount,1:5) = ...
                        wquadpt1*([dx1 dy1 (dx1^2) (dx1*dy1) (dy1^2)] - ...
                        consterms);
                    %Attribute it to "Aaux"
                    Aaux(icount + 1,1:5) = ...
                        wquadpt2*([dx2 dy2 (dx2^2) (dx2*dy2) (dy2^2)] - ...
                        consterms);
                %4th Order
                else
                    %Attribute it to "Aaux" (3rd order)
                    Aaux(icount,1:9) = ...
                        wquadpt1*([dx1 dy1 (dx1^2) (dx1*dy1) (dy1^2) ...
                        (dx1^3) (dx1^2)*dy1 dx1*(dy1^2) (dy1^3)] - ...
                        consterms);
                    %Attribute it to "Aaux"
                    Aaux(icount + 1,1:9) = ...
                        wquadpt2*([dx2 dy2 (dx2^2) (dx2*dy2) (dy2^2) ...
                        (dx2^3) (dx2^2)*dy2 dx2*(dy2^2) (dy2^3)] - ...
                        consterms);
                end  %End of IF (3rd or 4th)
                
                %Fill "waux"
                waux(icount:icount + 1) = [wquadpt1 wquadpt2];
                %Update "icount" (two quadrature points)
                icount = icount + 2;

            %For FIFTH or SIXTH ORDER (edge with known saturation value)
            elseif order == 5 || order == 6
                %Define the Gauss geometric weight, "ggw"
                ggw = 0.774596669241483;
                %Initialize "c"
                c = zeros(3,2);
                %Point nearest the first vertex:
                c(1,:) = ((a + b)/2) - ((b - a)/2)*ggw;
                %Midway point:
                c(2,:) = ((a + b)/2);
                %Point nearest the second vertex:
                c(3,:) = ((a + b)/2) + ((b - a)/2)*ggw;
                
                %Quadrature point 1:
                
                %Get the weight (weighted Least Squares). 
                %This scalar is called how many quadrature points exist.
                wquadpt1 = calcgeometricweight(elemevalcoord,c(1,:));
                %Define Component on "x" and on "y"
                delta1 = c(1,:) - elemevalcoord(1:2);
                dx1 = delta1(1);
                dy1 = delta1(2);
                
                %Quadrature point 2:
                
                %Get the weight (weighted Least Squares). 
                %This scalar is called how many quadrature points exist.
                wquadpt2 = calcgeometricweight(elemevalcoord,c(2,:));
                %Define Component on "x" and on "y"
                delta2 = c(2,:) - elemevalcoord(1:2);
                dx2 = delta2(1);
                dy2 = delta2(2);

                %Quadrature point 3:
                
                %Get the weight (weighted Least Squares). 
                %This scalar is called how many quadrature points exist.
                wquadpt3 = calcgeometricweight(elemevalcoord,c(3,:));
                %Define Component on "x" and on "y"
                delta3 = c(3,:) - elemevalcoord(1:2);
                dx3 = delta3(1);
                dy3 = delta3(2);
                
                %5th Order
                if order == 5
                     %Attribute it to "Aaux"
                     Aaux(icount,1:14) = ...
                        wquadpt1*([dx1 dy1 (dx1^2) (dx1*dy1) (dy1^2) ...
                        (dx1^3) (dx1^2)*dy1 dx1*(dy1^2) (dy1^3) ...
                        (dx1^4) (dx1^3)*dy1 (dx1^2)*(dy1^2) dx1*(dy1^3) ...
                        (dy1^4)] - consterms);
                    %Attribute it to "Aaux"
                    Aaux(icount + 1,1:14) = ...
                        wquadpt2*([dx2 dy2 (dx2^2) (dx2*dy2) (dy2^2) ...
                        (dx2^3) (dx2^2)*dy2 dx2*(dy2^2) (dy2^3) ...
                        (dx2^4) (dx2^3)*dy2 (dx2^2)*(dy2^2) dx2*(dy2^3) ...
                        (dy2^4)] - consterms);
                    %Attribute it to "Aaux"
                    Aaux(icount + 2,1:14) = ...
                        wquadpt3*([dx3 dy3 (dx3^2) (dx3*dy3) (dy3^2) ...
                        (dx3^3) (dx3^2)*dy3 dx3*(dy3^2) (dy3^3) ...
                        (dx3^4) (dx3^3)*dy3 (dx3^2)*(dy3^2) dx3*(dy3^3) ...
                        (dy3^4)] - consterms);
                %6th Order
                else
                     %Attribute it to "Aaux"
                     Aaux(icount,1:20) = ...
                        wquadpt1*([dx1 dy1 (dx1^2) (dx1*dy1) (dy1^2) ...
                        (dx1^3) (dx1^2)*dy1 dx1*(dy1^2) (dy1^3) ...
                        (dx1^4) (dx1^3)*dy1 (dx1^2)*(dy1^2) dx1*(dy1^3) ...
                        (dy1^4) (dx1^5) (dx1^4)*dy1 (dx1^3)*(dy1^2) ...
                        (dx1^2)*(dy1^3) dx1*(dy1^4) (dy1^5)] - consterms);
                    %Attribute it to "Aaux"
                    Aaux(icount + 1,1:20) = ...
                        wquadpt2*([dx2 dy2 (dx2^2) (dx2*dy2) (dy2^2) ...
                        (dx2^3) (dx2^2)*dy2 dx2*(dy2^2) (dy2^3) ...
                        (dx2^4) (dx2^3)*dy2 (dx2^2)*(dy2^2) dx2*(dy2^3) ...
                        (dy2^4) (dx2^5) (dx2^4)*dy2 (dx2^3)*(dy2^2) ...
                        (dx2^2)*(dy2^3) dx2*(dy2^4) (dy2^5)] - consterms);
                    %Attribute it to "Aaux"
                    Aaux(icount + 2,1:20) = ...
                        wquadpt3*([dx3 dy3 (dx3^2) (dx3*dy3) (dy3^2) ...
                        (dx3^3) (dx3^2)*dy3 dx3*(dy3^2) (dy3^3) ...
                        (dx3^4) (dx3^3)*dy3 (dx3^2)*(dy3^2) ...
                        dx3*(dy3^3) (dy3^4) (dx3^5) (dx3^4)*dy3 ...
                        (dx3^3)*(dy3^2) (dx3^2)*(dy3^3) dx3*(dy3^4) ...
                        (dy3^5)] - consterms);
                end  %End of IF (5th or 6th)
                
                %Fill "waux"
                waux(icount:icount + 2) = [wquadpt1 wquadpt2 wquadpt3];
                %Update "icount" (two quadrature points)
                icount = icount + 3;
            end  %End of IF (order)            
            
        %-----------------------------------
        %It is NOT a known saturation value:
        else
            %Get the coordinate of mirrored vertices:
            for ivert = 1:size(vertcoord,1)
                %Verify if the vertex evaluated is over the symmetry axe.
                if norm(vertcoord(ivert,:) - firstnodecoord) > tol && ...
                        norm(vertcoord(ivert,:) - secondnodecoord) > tol
                    %Get a mirrored coordinate.
                    vertcoord(ivert,1:2) = calcmirrorvec(firstnodecoord,...
                        secondnodecoord,vertcoord(ivert,:));
                end  %End of IF
            end  %End of FOR

            %Attribute the coordinate to "neighborcoordaux" (by mirror)
            neighborcoordaux(i,1:2) = mean(vertcoord(:,1:2),1);
        
            %Get the weight (weighted Least Squares). 
            %Vector has a length corresponding to amount og neighbours.
            waux(icount) = calcgeometricweight(elemevalcoord,...
                neighborcoordaux(i,:));

            %For SECOND ORDER
            if order == 2
                %Define "Aaux" according to "order"
                %Calculate the difference between the vectors "mirrorcoord" 
                %and "elemevalcoord".
                delta = neighborcoordaux(i,:) - elemevalcoord(1:2);
                %Component on "x" and on "y"
                Aaux(icount,1:2) = waux(icount)*delta;

                %Update "icount"
                icount = icount + 1;
            
            end
        end  %End of IF (is the saturation known on the edge?)
    end  %End of FOR

    %Concatenate the matrices "A" and "w"
    A = vertcat(Ain,Aaux);
    w = vertcat(win,waux');
end  %End of IF

%--------------------------------------------------------------------------
%Function "calcleastsquarematrix"
%--------------------------------------------------------------------------

function [A,w,consterms] = calcleastsquarematrix(ielem,bedgrownum,...
    flagknownedge,injecelem)
%Define global parameters:
global coord elem centelem order lsneightype;

%Catch the elements surrounding the element evaluated.
[esureface,esurefull] = getsurelem(ielem);
%Verify if it is a triangle ("elemtype" = 3) or a quad. ("elemtype" = 4)
elemtype = sum(elem(ielem,1:4) ~= 0); 
%Get the centroid of the element evaluated ("elemeval")
elemevalcoord = centelem(ielem,:);
%Get the vertices of the neighbour element.
verticeseval = elem(ielem,1:elemtype);
%Get the vertices coordinate
vertcoordeval = coord(verticeseval,:);
%Initialize "consterms"
consterms = 0;

%Choose the strategy according "order" and "lsneightype". See additional 
%comments in the head.
%SECOND ORDER:
%It uses Face Neighbor Least Square (Computational Fluid Dynamics, 
%Bazek, see cap. 5, pp. 162)
if order == 2 && lsneightype == 1
    %Initialize "A".
    A = zeros(length(esureface),2);
            
    %Get the neighboring elements coordinate
    neighcoord = centelem(esureface,1:2);
    %Get the weight (weighted Least Squares). Vector [length("esureface")].
    w = calcgeometricweight(elemevalcoord,neighcoord);

    %Attribute to "A" the difference between position in each element
    %surround "ielem" and the "ielem" position.
    %Component on "x"
    A(1:length(esureface),1) = ...
        w.*(centelem(esureface,1) - elemevalcoord(1));
    %Component on "y"
    A(1:length(esureface),2) = ...
        w.*(centelem(esureface,2) - elemevalcoord(2));
        
%SECOND ORDER:
%It uses Full Neighbor (faces and vertices) Least Square (Computational 
%Fluid Dynamics, Bazek, see cap. 5, pp. 162)
elseif order == 2 && lsneightype == 2
    %Initialize "A".
    A = zeros(length(esurefull),2);
        
    %Get the neighboring elements coordinate
    neighcoord = centelem(esurefull,1:2);
    %Get the weight (weighted Least Squares). Vector [length("esureface")].
    w = calcgeometricweight(elemevalcoord,neighcoord);

    %Attribute to "A" the difference between position in each element
    %surround "ielem" and the "ielem" position.
    %Component on "x"
    A(1:length(esurefull),1) = ...
        w.*(centelem(esurefull,1) - elemevalcoord(1));
    %Component on "y"
    A(1:length(esurefull),2) = ...
        w.*(centelem(esurefull,2) - elemevalcoord(2));
        
%THIRD ORDER (triangle) - Face Neighbour and Face Neighbour of Face Neig.:
%It is a special case (Third-order). It uses Face Neighbor of the 
%Neighbour by Face Least Square (see Goosh and Van Altena, 2002)
elseif order == 3 && elemtype == 3 && lsneightype == 1
    %Evaluate the Neighboring of each element in the "ielem" Face
    %Neighboring
    for i = 1:length(esureface)
        %Catch the elements surrounding the ith neighbour.
        [esurefaceneighb,] = getsurelem(esureface(i));
        %Exclude "ielem" of "esurefaceneighb"
        esurefaceneighb = ...
            esurefaceneighb(logical(esurefaceneighb ~= ielem));
        %Get the union between "esureface" and "esurefaceneighb" for the 
        %neighbour element evaluated
        esureface = union(esureface,esurefaceneighb,'stable');
    end  %End of FOR
    
    %Get the neighboring elements coordinate
    neighcoord = centelem(esureface,1:2);
    %Initialize "A" and "w"
    A = zeros(length(esureface),5);
    w = zeros(length(esureface),1);
    
    %Get the coefficients for the first row (the control volume evaluated).
    %-x:
    coefevalx = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,0);
    %-y:
    coefevaly = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,1);
    %-x2:
    coefevalx2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,2,0);
    %-xy:
    coefevalxy = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,1);
    %-y2:
    coefevaly2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,2);

    %Define the constraint terms
    consterms = [coefevalx coefevaly coefevalx2 coefevalxy coefevaly2];
    
    %Swepts all neighbour elements
    for i = 1:length(esureface)
        %Get the columns of the neighbour element 
        elemcolumns = elem(esureface(i),1:4);
        %Get the vertices of the neighbour element.
        vertices = elemcolumns(logical(elemcolumns ~= 0));
        %Get the vertices coordinate
        vertcoord = coord(vertices,:);
        
        %^x:
        coefneigx = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esureface(i),1,0);
        %^y:
        coefneigy = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esureface(i),0,1);
        %^x2:
        coefneigx2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esureface(i),2,0);
        %^xy:
        coefneigxy = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esureface(i),1,1);
        %^y2:
        coefneigy2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esureface(i),0,2);
    
        %Get the weight (weighted Least Squares). 
        %Vector [length("esureface")].
        w(i) = calcgeometricweight(elemevalcoord,neighcoord(i,:));
    
        %Attribute to "A" the coefficients "^" (numerical integrals).
        %Two first columns (Associated to Gradient terms); 
        %Three last columns (%Hessian Contribution):
        A(i,:) = w(i)*([coefneigx coefneigy coefneigx2 coefneigxy ...
            coefneigy2] - consterms);
    end  %End of FOR
    
%THIRD ORDER (triangle) - Full Neighbour or THIRD ORDER (quadrangle) - Full 
%Neighbour (MANDATORILY):
%It is a special case (Third-order). It uses Full Neighbour (see Mandal and 
%Rao, 2011)
elseif (order == 3 && elemtype == 3 && lsneightype == 2) || ...
        (order == 3 && elemtype == 4)
    %Get the neighboring elements coordinate
    neighcoord = centelem(esurefull,1:2);
    %Initialize "A" and "w"
    A = zeros(length(esurefull),5);
    w = zeros(length(esurefull),1);

    %Get the coefficients for the first row (the control volume evaluated).
    %-x:
    coefevalx = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,0);
    %-y:
    coefevaly = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,1);
    %-x2:
    coefevalx2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,2,0);
    %-xy:
    coefevalxy = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,1);
    %-y2:
    coefevaly2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,2);
    
    %Define the constraint terms
    consterms = [coefevalx coefevaly coefevalx2 coefevalxy coefevaly2];

    %Swepts all neighbour elements
    for i = 1:length(esurefull)
        %Get the columns of the neighbour element 
        elemcolumns = elem(esurefull(i),1:4);
        %Get the vertices of the neighbour element.
        vertices = elemcolumns(logical(elemcolumns ~= 0));
        %Get the vertices coordinate
        vertcoord = coord(vertices,:);

        %Get the coefficients ("d" means "domain", CV inside the domain):
        %^x:
        coefneigx = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esurefull(i),1,0);
        %^y:
        coefneigy = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esurefull(i),0,1);
        %^x2:
        coefneigx2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esurefull(i),2,0);
        %^xy:
        coefneigxy = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esurefull(i),1,1);
        %^y2:
        coefneigy2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esurefull(i),0,2);
    
        %Get the weight (weighted Least Squares). 
        %Vector [length("esureface")].
        w(i) = calcgeometricweight(elemevalcoord,neighcoord(i,:));

        %Attribute to "A" the coefficients "^" (numerical integrals).
        %Two first columns (Associated to Gradient terms); 
        %Three last columns (%Hessian Contribution):
        A(i,:) = w(i)*([coefneigx coefneigy coefneigx2 coefneigxy ...
            coefneigy2] - consterms);
    end  %End of FOR
    
%FOURTH ORDER (triangle) - Full Neighbour MANDATORILY. See Goosh and Van 
%Altena (2002)
elseif order == 4 && elemtype == 3
    %Initialize "esuremod"
    esuremod = esureface;

    %Get the neighboring elements coordinate
    neighcoord = centelem(esurefull,1:2);
    %Initialize "A" and "w"
    A = zeros(length(esurefull),9);
    w = zeros(length(esurefull),1);

    %Get the coefficients for the first row (the control volume evaluated).
    %-x:
    coefevalx = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,0);
    %-y:
    coefevaly = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,1);
    %-x2:
    coefevalx2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,2,0);
    %-xy:
    coefevalxy = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,1);
    %-y2:
    coefevaly2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,2);
    %-x3:
    coefevalx3 = getintegdiff(ielem,vertcoordeval,elemevalcoord,3,0);
    %-x2y:
    coefevalx2y = getintegdiff(ielem,vertcoordeval,elemevalcoord,2,1);
    %-xy2:
    coefevalxy2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,2);
    %-y3:
    coefevaly3 = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,3);
        
    %Define the constraint terms
    consterms = [coefevalx coefevaly coefevalx2 coefevalxy coefevaly2 ...
        coefevalx3 coefevalx2y coefevalxy2 coefevaly3];

    %Swepts all neighbour elements
    for i = 1:length(esurefull)
        %Get the columns of the neighbour element 
        elemcolumns = elem(esurefull(i),1:4);
        %Get the vertices of the neighbour element.
        vertices = elemcolumns(logical(elemcolumns ~= 0));
        %Get the vertices coordinate
        vertcoord = coord(vertices,:);
        
        %Get the coefficients ("d" means "domain", CV inside the domain):
        %^x:
        coefneigx = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esurefull(i),1,0);
        %^y:
        coefneigy = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esurefull(i),0,1);
        %^x2:
        coefneigx2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esurefull(i),2,0);
        %^xy:
        coefneigxy = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esurefull(i),1,1);
        %^y2:
        coefneigy2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esurefull(i),0,2);
        %^x3:
        coefneigx3 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esurefull(i),3,0);
        %^x2y:
        coefneigx2y = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esurefull(i),2,1);
        %^xy2:
        coefneigxy2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esurefull(i),1,2);
        %^y3:
        coefneigy3 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esurefull(i),0,3);
    
        %Get the weight (weighted Least Squares). 
        %Vec. [length("esurefull")].
        w(i) = calcgeometricweight(elemevalcoord,neighcoord(i,:));

        %Attribute to "A" the coefficients "^" (numerical integrals).
        %Two first columns (Associated to Gradient terms); 
        %Three last columns (%Hessian Contribution):
        A(i,:) = w(i)*([coefneigx coefneigy coefneigx2 coefneigxy ...
            coefneigy2 coefneigx3 coefneigx2y coefneigxy2 coefneigy3] - ...
            consterms);
    end  %End of FOR
    
%FOURTH ORDER (quadrangle) - Face Neighbour for the element evaluated and 
%Face Neighbour for each element of Neighboring by Face. See Goosh and Van 
%Altena (2002)
elseif order == 4 && elemtype == 4
    %Get the stencil according amount of neighbours and type of element.
    [esuremod] = getstencil(order,ielem,elemtype,esureface,esurefull,...
        injecelem);

    %Get the neighboring elements coordinate
    neighcoord = centelem(esuremod,1:2);
    %Initialize "A" and "w"
    A = zeros(length(esuremod),9);
    w = zeros(length(esuremod),1);

    %Get the coefficients for the first row (the control volume evaluated).
    %-x:
    coefevalx = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,0);
    %-y:
    coefevaly = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,1);
    %-x2:
    coefevalx2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,2,0);
    %-xy:
    coefevalxy = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,1);
    %-y2:
    coefevaly2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,2);
    %-x3:
    coefevalx3 = getintegdiff(ielem,vertcoordeval,elemevalcoord,3,0);
    %-x2y:
    coefevalx2y = getintegdiff(ielem,vertcoordeval,elemevalcoord,2,1);
    %-xy2:
    coefevalxy2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,2);
    %-y3:
    coefevaly3 = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,3);
        
    %Define the constraint terms
    consterms = [coefevalx coefevaly coefevalx2 coefevalxy coefevaly2 ...
        coefevalx3 coefevalx2y coefevalxy2 coefevaly3];

    %Swepts all neighbour elements
    for i = 1:length(esuremod)
        %Get the columns of the neighbour element 
        elemcolumns = elem(esuremod(i),1:4);
        %Get the vertices of the neighbour element.
        vertices = elemcolumns(logical(elemcolumns ~= 0));
        %Get the vertices coordinate
        vertcoord = coord(vertices,:);
    
        %Get the coefficients ("d" means "domain", CV inside the domain):
        %^x:
        coefneigx = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),1,0);
        %^y:
        coefneigy = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),0,1);
        %^x2:
        coefneigx2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),2,0);
        %^xy:
        coefneigxy = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),1,1);
        %^y2:
        coefneigy2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),0,2);
        %^x3:
        coefneigx3 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),3,0);
        %^x2y:
        coefneigx2y = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),2,1);
        %^xy2:
        coefneigxy2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),1,2);
        %^y3:
        coefneigy3 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),0,3);
    
        %Get the weight (weighted Least Squares). 
        %Vec. [length("esurefacemod")]
        w(i) = calcgeometricweight(elemevalcoord,neighcoord(i,:));

        %Attribute to "A" the coefficients "^" (numerical integrals).
        %Two first columns (Associated to Gradient terms); 
        %Three last columns (%Hessian Contribution):
        A(i,:) = w(i)*([coefneigx coefneigy coefneigx2 coefneigxy ...
            coefneigy2 coefneigx3 coefneigx2y coefneigxy2 coefneigy3] - ...
            consterms);
    end  %End of FOR

%FIFTH (triangle or quadrangle) - See the function "getstencil".
elseif order == 5
    %Get the stencil according amount of neighbours and type of element.
    [esuremod] = getstencil(order,ielem,elemtype,esureface,esurefull);

    %Get the neighboring elements coordinate
    neighcoord = centelem(esuremod,1:2);
    %Initialize "A" and "w"
    A = zeros(length(esuremod),14);
    w = zeros(length(esuremod),1);

    %Get the coefficients for the first row (the control volume evaluated).
    %-x:
    coefevalx = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,0);
    %-y:
    coefevaly = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,1);
    %-x2:
    coefevalx2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,2,0);
    %-xy:
    coefevalxy = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,1);
    %-y2:
    coefevaly2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,2);
    %-x3:
    coefevalx3 = getintegdiff(ielem,vertcoordeval,elemevalcoord,3,0);
    %-x2y:
    coefevalx2y = getintegdiff(ielem,vertcoordeval,elemevalcoord,2,1);
    %-xy2:
    coefevalxy2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,2);
    %-y3:
    coefevaly3 = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,3);
    %-x4:
    coefevalx4 = getintegdiff(ielem,vertcoordeval,elemevalcoord,4,0);
    %-x3y:
    coefevalx3y = getintegdiff(ielem,vertcoordeval,elemevalcoord,3,1);
    %-x2y2:
    coefevalx2y2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,2,2);
    %-xy3:
    coefevalxy3 = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,3);
    %-y4:
    coefevaly4 = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,4);

    %Define the constraint terms
    consterms = [coefevalx coefevaly coefevalx2 coefevalxy coefevaly2 ...
        coefevalx3 coefevalx2y coefevalxy2 coefevaly3 coefevalx4 ...
        coefevalx3y coefevalx2y2 coefevalxy3 coefevaly4];

    %Swepts all neighbour elements
    for i = 1:length(esuremod)
        %Get the columns of the neighbour element 
        elemcolumns = elem(esuremod(i),1:4);
        %Get the vertices of the neighbour element.
        vertices = elemcolumns(logical(elemcolumns ~= 0));
        %Get the vertices coordinate
        vertcoord = coord(vertices,:);
    
        %Get the coefficients ("d" means "domain", CV inside the domain):
        %^x:
        coefneigx = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),1,0);
        %^y:
        coefneigy = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),0,1);
        %^x2:
        coefneigx2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),2,0);
        %^xy:
        coefneigxy = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),1,1);
        %^y2:
        coefneigy2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),0,2);
        %^x3:
        coefneigx3 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),3,0);
        %^x2y:
        coefneigx2y = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),2,1);
        %^xy2:
        coefneigxy2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),1,2);
        %^y3:
        coefneigy3 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),0,3);
        %^x4:
        coefneigx4 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),4,0);
        %^x2y:
        coefneigx3y = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),3,1);
        %^x2y2:
        coefneigx2y2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),2,2);
        %^xy3:
        coefneigxy3 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),1,3);
        %^y4:
        coefneigy4 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),0,4);
    
        %Get the weight (weighted Least Squares). 
        %Vec. [length("esurefacemod")]
        w(i) = calcgeometricweight(elemevalcoord,neighcoord(i,:));

        %Attribute to "A" the coefficients "^" (numerical integrals).
        %Two first columns (Associated to Gradient terms); 
        %Three last columns (%Hessian Contribution):
        A(i,:) = w(i)*([coefneigx coefneigy coefneigx2 coefneigxy ...
            coefneigy2 coefneigx3 coefneigx2y coefneigxy2 coefneigy3 ...
            coefneigx4 coefneigx3y coefneigx2y2 coefneigxy3 coefneigy4] - ...
            consterms);
    end  %End of FOR

%SIXTH (triangle or quadrangle) - See the function "getstencil".
elseif order == 6
    %Get the stencil according amount of neighbours and type of element.
    [esuremod] = getstencil(order,ielem,elemtype,esureface,esurefull);

    %Get the neighboring elements coordinate
    neighcoord = centelem(esuremod,1:2);
    %Initialize "A" and "w"
    A = zeros(length(esuremod),20);
    w = zeros(length(esuremod),1);

    %Get the coefficients for the first row (the control volume evaluated).
    %-x:
    coefevalx = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,0);
    %-y:
    coefevaly = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,1);
    %-x2:
    coefevalx2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,2,0);
    %-xy:
    coefevalxy = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,1);
    %-y2:
    coefevaly2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,2);
    %-x3:
    coefevalx3 = getintegdiff(ielem,vertcoordeval,elemevalcoord,3,0);
    %-x2y:
    coefevalx2y = getintegdiff(ielem,vertcoordeval,elemevalcoord,2,1);
    %-xy2:
    coefevalxy2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,2);
    %-y3:
    coefevaly3 = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,3);
    %-x4:
    coefevalx4 = getintegdiff(ielem,vertcoordeval,elemevalcoord,4,0);
    %-x3y:
    coefevalx3y = getintegdiff(ielem,vertcoordeval,elemevalcoord,3,1);
    %-x2y2:
    coefevalx2y2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,2,2);
    %-xy3:
    coefevalxy3 = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,3);
    %-y4:
    coefevaly4 = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,4);
    %-x5:
    coefevalx5 = getintegdiff(ielem,vertcoordeval,elemevalcoord,5,0);
    %-x4y:
    coefevalx4y = getintegdiff(ielem,vertcoordeval,elemevalcoord,4,1);
    %-x3y2:
    coefevalx3y2 = getintegdiff(ielem,vertcoordeval,elemevalcoord,3,2);
    %-x2y3:
    coefevalx2y3 = getintegdiff(ielem,vertcoordeval,elemevalcoord,2,3);
    %-xy4:
    coefevalxy4 = getintegdiff(ielem,vertcoordeval,elemevalcoord,1,4);
    %-y5:
    coefevaly5 = getintegdiff(ielem,vertcoordeval,elemevalcoord,0,5);

    %Define the constraint terms
    consterms = [coefevalx coefevaly coefevalx2 coefevalxy coefevaly2 ...
        coefevalx3 coefevalx2y coefevalxy2 coefevaly3 coefevalx4 ...
        coefevalx3y coefevalx2y2 coefevalxy3 coefevaly4 coefevalx5 ...
        coefevalx4y coefevalx3y2 coefevalx2y3 coefevalxy4 coefevaly5];

    %Swepts all neighbour elements
    for i = 1:length(esuremod)
        %Get the columns of the neighbour element 
        elemcolumns = elem(esuremod(i),1:4);
        %Get the vertices of the neighbour element.
        vertices = elemcolumns(logical(elemcolumns ~= 0));
        %Get the vertices coordinate
        vertcoord = coord(vertices,:);
    
        %Get the coefficients ("d" means "domain", CV inside the domain):
        %^x:
        coefneigx = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),1,0);
        %^y:
        coefneigy = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),0,1);
        %^x2:
        coefneigx2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),2,0);
        %^xy:
        coefneigxy = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),1,1);
        %^y2:
        coefneigy2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),0,2);
        %^x3:
        coefneigx3 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),3,0);
        %^x2y:
        coefneigx2y = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),2,1);
        %^xy2:
        coefneigxy2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),1,2);
        %^y3:
        coefneigy3 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),0,3);
        %^x4:
        coefneigx4 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),4,0);
        %^x2y:
        coefneigx3y = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),3,1);
        %^x2y2:
        coefneigx2y2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),2,2);
        %^xy3:
        coefneigxy3 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),1,3);
        %^y4:
        coefneigy4 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),0,4);
        %^x5:
        coefneigx5 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),5,0);
        %^x4y:
        coefneigx4y = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),4,1);
        %^x3y2:
        coefneigx3y2 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),3,2);
        %^x2y3:
        coefneigx2y3 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),2,3);
        %^xy4:
        coefneigxy4 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),1,4);
        %^y5:
        coefneigy5 = getcoefnuminteg(elemevalcoord,neighcoord(i,:),...
            vertcoord,esuremod(i),0,5);
    
        %Get the weight (weighted Least Squares). 
        %Vec. [length("esurefacemod")]
        w(i) = calcgeometricweight(elemevalcoord,neighcoord(i,:));

        %Attribute to "A" the coefficients "^" (numerical integrals).
        %Two first columns (Associated to Gradient terms); 
        %Three last columns (%Hessian Contribution):
        A(i,:) = w(i)*([coefneigx coefneigy coefneigx2 coefneigxy ...
            coefneigy2 coefneigx3 coefneigx2y coefneigxy2 coefneigy3 ...
            coefneigx4 coefneigx3y coefneigx2y2 coefneigxy3 coefneigy4 ...
            coefneigx5 coefneigx4y coefneigx3y2 coefneigx2y3 coefneigxy4 ...
            coefneigy5] - consterms);
    end  %End of FOR
end  %End of IF ("order", "lsneightype" and element type)

%-----------------------------------
%Evaluate the boundary contribution:

if elemtype ~= length(esureface)
    %We need of the "ielem" and the number of element neighbour.
    [A,w] = attriboundcontrib(A,w,order,ielem,elemevalcoord,consterms,...
        bedgrownum,flagknownedge);
end  %End of IF

%There is boundary treatment for other elements which not the "elemeval"
if order == 4 && size(A,1) < 16 
    %Swept all face neighbour (or full stencil, above 4th order) of "ielem"
    for ineigh = 1:length(esuremod)
        %Catch the elements surrounding the element evaluated.
        [boundesureface,] = getsurelem(esuremod(ineigh));
        %Verify if it is a triangle ("elemtype" = 3) or a quad. ("elemtype" = 4)
        boundelemtype = sum(elem(esuremod(ineigh),1:4) ~= 0); 
        if boundelemtype ~= length(boundesureface)
            %Verigy the vicinity of each "esureface" element
            [A,w] = attriboundcontrib(A,w,order,esuremod(ineigh),...
                elemevalcoord,consterms,bedgrownum,flagknownedge);
        end %End of IF
    end  %End of FOR
end  %End of IF

%--------------------------------------------------------------------------
%Function "gramschimidt"
%--------------------------------------------------------------------------

function [gs] = gramschimidt(A)
%Get amount of row and column in matrix "A"
[row,column] = size(A);
%Initialize "Q" and "R"
Q = zeros(row,column);
R = zeros(column);

%Gram-Schmidt Algorithm (Fill "Q" and "R")
for j = 1:column
    v = A(:,j);
    for i = 1:j - 1
        R(i,j) = Q(:,i)'*A(:,j);
        v = v - R(i,j)*Q(:,i);
    end  %End of FOR
    R(j,j) = norm(v);
    Q(:,j) = v/R(j,j);
end  %End of FOR
    
%Calculate the Gram-Schimidt coeffiente ("gs")
gs = R\Q';

%--------------------------------------------------------------------------
%Function "getleastsquarescoeff"
%--------------------------------------------------------------------------

function [wvector,wmap,constraint,lsw,amountofneigvec] = ...
    getleastsquarescoeff(flagknownedge,injecelem)
%Define global parameters:
global elem bedge order;

%Initialize "elemsize"
elemsize = size(elem,1);
%Initialize "bedgrownum". It is a vector with 1:size of "bedge" 
bedgrownum = 1:size(bedge,1);
%Initialize "wmap" and counters.
wmap = zeros(elemsize + 1,1);
countw = 0;
%Initialize "constraint"
columns = (order == 2)*2 + (order == 3)*5 + (order == 4)*9 + ...
    (order == 5)*14 + (order == 6)*20;
constraint = zeros(elemsize,columns);
%Initialize "amountofneigvec". It is a vector that stores the amount of
%neighbour used for the recosntruction in each element.
amountofneigvec = zeros(elemsize,1);

%Initialize counter
m = 0;
%Swept all elements
for ielem = 1:elemsize
    %Initialize "wvectoraux"
    wvectoraux = 0;

    %Get the matrix "A"
    [A,w,consterms] = ...
        calcleastsquarematrix(ielem,bedgrownum,flagknownedge,injecelem);
    
    %Fill "amountofneigvec"
    amountofneigvec(ielem) = size(A,1);
    %It stores constraint parameters:
    constraint(ielem,:) = consterms; 
    
    %Calculate the matrices "Q" and "R" (Gram-Schmidit Orthog.)
    gs = gramschimidt(A);
    
    %Convert the matrix "gs" for a vector
    wvectoraux(1:size(A,2)*size(gs,2)) = gs';
    %Fill the vectors "wvector"
    wvector(countw + 1:countw + length(wvectoraux)) = wvectoraux;
    %Stores the Least Squares weight
    lsw(m + 1:m + length(w)) = w;
    
    %Update "countw" and "wmap"
    countw = countw + length(wvectoraux);
    wmap(ielem + 1) = countw;
    %Update "m"
    m = m + length(w);
end  %End of FOR (swept all elements)

%It clear "bedgrownum"
clear bedgrownum;












