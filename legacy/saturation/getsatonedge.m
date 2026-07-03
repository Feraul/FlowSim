%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with
%high-order resolution
%Type of file: FUNCTION
%Criate date: 12/01/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%Modified: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals:
%

%--------------------------------------------------------------------------
%Additional comments:
%"taylorterms" is acctualy the grad of saturation (2nd order applications).
%When we use a third-order scheme, we do not have only a grad. We have also
%a hessian matrix. Thus, for a more generalized understanding we named the
%matrix like that.

%--------------------------------------------------------------------------

function[swonedgelim,swonedgenolim] = getsatonedge(elemeval,vertices,verticescoord,...
    taylorterms,Con,limiterflag,order,constraint,flagknownvert,...
    satonvertices,mlpbyelem,centelemeval)

%Define the strategy recovery according the "smethod" flag
%First-Order Upwind
if order == 1
    %Recovery "Sw" by a simply strategy Upwind.
    
    swonedgelim = Con(elemeval(1));
    swonedgenolim=0;
    %Higher-Order Schemes (Second Order).
elseif order == 2
    %Calculate the second-order saturation value.
    [swonedgelim,swonedgenolim] = get2ndorder(elemeval,vertices,verticescoord,...
        taylorterms(:,1:2),Con,limiterflag,flagknownvert,satonvertices,...
        mlpbyelem,centelemeval);
    %Very Higher-Order Schemes (Third and Fourth Order)
elseif order == 3 || order == 4
    %Define "columns" according to "order"
    columns = 5*(order == 3) + 9*(order == 4);
    %Calculate the third or fourth order saturation value
    swonedge = get3rdand4thorder(order,elemeval,vertices,verticescoord,...
        taylorterms(:,1:columns),Con,constraint(:,1:columns),mlpbyelem,...
        centelemeval);
    
    %Very Higher-Order Schemes (Fifth and Sixth Order)
elseif order == 5 || order == 6
    %Calculate the third or fourth order saturation value
    swonedge = get5thand6thorder(order,elemeval,vertices,verticescoord,...
        taylorterms,Con,constraint,mlpbyelem);
end  %End of IF ("order")

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION:
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "get2ndorder"
%--------------------------------------------------------------------------

function [swonedgelim,swonedgenolim] = get2ndorder(elemeval,vertices,verticescoord,...
    taylorterms,Sw,limiterflag,flagknownvert,satonvertices,mlpbyelem,...
    centelemeval)
%Define global parameters:
global recovtype ;

swonedgeaux=0;
%Initialize parameter "k" (It is "-1" for Upwind approximation)
%and "epsilon" is the machine limit (1e-16)
k = limiterflag{8};
epsilon = 1e-16;
%Get the key for the Edge-based limiter
eblimkey = limiterflag{1};
%Get the type of Edge-based limiter
eblimtype = limiterflag{2};
%Get the key for the Element-based limiter
cvblimkey = limiterflag{3};
%Get the type of Element-based limiter
cvbtype = limiterflag{4};
%Get the tune parameter
delta = limiterflag{5};
%Greek Correction
gckey = limiterflag{7};

%Initialize points coordinate:
pointl = centelemeval(1,:);
%element inside the domain
booleanelem = (length(elemeval) > 1);
pointr = centelemeval(2*booleanelem + (1 - booleanelem),:);

%Verify if exists two verticescoord:
pointi = verticescoord(1,1:2);
booleanvtx = (length(verticescoord) > 1);
pointj = verticescoord(2*booleanvtx + (1 - booleanvtx),1:2);

%Calculate the midedge
midedge = 0.5*(pointi + pointj );

%Initialize "phi" and "grad". Each row is a element gradient.
phi = ones(length(elemeval),1);
grad = taylorterms(elemeval,:);
gradnonlimit=taylorterms(elemeval,:);
%--------------------------------------------------------------------------
%Verify the type of reconstruction (Green-Gauss or Least-Squares)

%------------------------------------
%Green-Gauss (Durlofsky et al., 1992)
if strcmp(recovtype,'ggg')
    %Get the distance between the cell-centered and the midedge:
    %Get the vector distance between the midedge point and the
    %cell-centered point.
    rcm = midedge - centelemeval(1,:);
    
    %Calculate the saturation through the edge (only by projection
    %from cell-center till middle-point)
    swonedge = Sw(elemeval(1)) + dot(grad(1,1:2),rcm);
    
    [esureface,] = getsurelem(elemeval(1));
    %Verify if the Durlofsky Scheme is in use and if there is an extrema
    if swonedge > max(Sw([elemeval(1) esureface])) || ...
            swonedge < min(Sw([elemeval(1) esureface]))
        
        %     if strcmp(recovtype,'ggg') && (swonedge > max(Sw(elemeval)) || ...
        %             swonedge < min(Sw(elemeval)))
        %Put the approximation in FIRST ORDER
        swonedge = Sw(elemeval(1));
    end  %End of IF (Durlofsky)
    
    %----------------------------
    %Least-Squares (Blazek, 2001)
else
    %----------------------------------------------------------------------
    %Verify if there is an Edge-Based Limiter and Which one it is.
    if strcmp(eblimkey,'on')
        %Verify if there is Element-Based limiter
        if strcmp([cvblimkey cvbtype],'onwf')
            %Swept the amount of elements shared by the edge evaluated
            for ishare = 1:length(elemeval)
                %Get Element-based limiter
                phi(ishare) = elemlimiter(elemeval(ishare),Sw,delta,...
                    flagknownvert,satonvertices);
                %Get the gradient (reconstructed by MUSCL strategy) for
                %each element evaluated (in each matrix row)
                grad(ishare,:) = phi(ishare)*grad(ishare,:);
                
            end  %End of FOR
        end  %End of IF (Is there element based limiter?)
        %Calculate the gradient centered ("gradcent")
        pos = 2*booleanelem + (1 - booleanelem);
        % MUSCL and MUSCL-pos-limiter
        if strcmp(limiterflag{6},'on') || strcmp(limiterflag{11},'on')
            %Calculate the gradients centered and upwind.
            %Calculate the vector distance between colocation points.
            dlr = pointr - pointl;            
            gradcent = Sw(elemeval(pos)) - Sw(elemeval(1));
            %Calculate the gradient upwind ("gradupwd"). First "grad"
            %matrix row.
            gradupwd = 2*dot(grad(1,:),dlr) - gradcent;

            %Define "r" (smoothness factor). "epsilon" avoids division by 0
            r = gradcent/(gradupwd + epsilon);
            
            %Define Edge-Based Limiter
            
            qsi = max(0,edgelimiter(r,eblimtype));
            
            %Calculate the gradient limited according "k" value.
            gradslimited = qsi*0.5*((1 - k)*gradupwd + (1 + k)*gradcent);
            swonedgelim = Sw(elemeval(1)) + 0.5*gradslimited;
            %For Modified MUSCL
            if strcmp(limiterflag{11},'on')
                dlr = pointr - pointl;
                gradcent = Sw(elemeval(pos)) - Sw(elemeval(1));
                gradupwdnonlimit= 2*dot(gradnonlimit(1,:),dlr) - gradcent;
                gradsnonlimit = 0.5*((1 - k)*gradupwdnonlimit + (1 + k)*gradcent);
                swonedgenolim=Sw(elemeval(1)) + 0.5*gradsnonlimit;
            else
                swonedgenolim=0;
                
            end
            % For LP MUSCL
        elseif strcmp(limiterflag{13},'on')
            %Calculate the gradients centered and upwind.
            %Calculate the vector distance between colocation points.
            
            dlr = pointr - pointl;
            gradcent = Sw(elemeval(pos)) - Sw(elemeval(1));
            %Calculate the gradient upwind ("gradupwd"). First "grad"
            %matrix row.
            gradupwd = 2*dot(grad(1,:),dlr) - gradcent;
            
            %Define "r" (smoothness factor). "epsilon" avoids division by 0
            r = gradcent/(gradupwd + epsilon);
            
            %Define Edge-Based Limiter
            
            qsi = max(0,edgelimiter(r,eblimtype));
            %Calculate the gradients centered and upwind.
            %Calculate the vector distance between colocation points.
            dlr = midedge - pointl;
            vetor= 2*midedge-(pointr + pointl);
            
            pos = 2*booleanelem + (1 - booleanelem);
            
            if length(elemeval)>1
                Swaux =Sw(elemeval(pos))+ dot(gradnonlimit(2,:),vetor);
                a = 0.5*(Swaux - Sw(elemeval(1)));
                gradslimited = qsi*(dot(gradnonlimit(1,:),dlr) + qsi*k*(a-dot(gradnonlimit(1,:),dlr)));
                %gradslimited = (dot(gradnonlimit(1,:),dlr) + k*(a-dot(gradnonlimit(1,:),dlr)));
                swonedgelim = Sw(elemeval(1)) + gradslimited;
                swonedgenolim=0;
            else
                % edge element
                Swaux =Sw(elemeval(pos));
                a = 0.5*(Swaux - Sw(elemeval(1)));
                gradslimited = qsi*(k*a+(1-k)*dot(gradnonlimit(1,:),dlr)); % método de Burg (eq. 3)
                %gradslimited = k*a+(1-k)*dot(gradnonlimit(1,:),dlr); % método de Burg (eq. 3)
                swonedgelim = Sw(elemeval(1)) + gradslimited;
                swonedgenolim=0;
            end
        end
        %------------------------------------------------------------------
        %Verify if there is a Greek Correction (Delis and Nikolos, 2012).
        
        %There is Greek Correction.
        if strcmp(gckey,'on')
            %Get the crossing point between the vectors "dlr" and "dij"
            %Define the angular coefficients:
            %"rl" segment
            mlr = ...
                (pointr(2) - pointl(2))/(pointr(1) - pointl(1) + epsilon);
            %"ij" segment
            mij = ...
                (pointj(2) - pointi(2))/(pointj(1) - pointi(1) + epsilon);
            
            %"x" component
            crosspoint(1) = (pointl(2) - pointi(2) + pointi(1)*mij - ...
                pointl(1)*mlr)/(mij - mlr + epsilon);
            %"y" component
            crosspoint(2) = pointl(2) + mlr*(crosspoint(1) - pointl(1));
            
            %Get the distance "dlcross"
            dlcross = crosspoint - pointl;
            
            %Finaly, define the "ratio" and calculate the distance between
            %the points "D" and "M"
            ratio = norm(dlcross)/norm(dlr);
            
            %Get the saturation on "D" point
            swonedge = Sw(elemeval(1)) + ratio*gradslimited;
            
            %Catch the middle edge:
            midedge = 0.5*(pointi + pointj);
            %The distance between the point "D" and the point "M"
            rdm = midedge - crosspoint;
            
            %Evaluate if the Greek Correction is necessary.
            if norm(rdm) > 1e-14
                swonedge = bringtomidpoint(vertices,verticescoord,elemeval,...
                    taylorterms,rdm,Sw,swonedge,limiterflag,k);
            end  %End of IF
            
            %There is NO Greek Correction.
        end  %End of IF (Greek Correction)
        
        
        %----------------------------------------------------------------------
        %     elseif strcmp(eblimkey,'on') && strcmp(eblimtype,'ebmlp')
        %         %Get the limiter
        %         phi = calcEBMLP(length(elemeval),elemeval(1),vertices,Sw,...
        %             flagknownvert,grad(1,1:2),satonvertices,centelemeval(1,:),...
        %             verticescoord,'mlp');
        %
        %         grad = phi*grad(1,:);
        %
        %         %Get the vector distance between the midedge point and the
        %         %cell-centered point.
        %         rcm = midedge - centelemeval(1,:);
        %
        %         %Calculate the saturation through the edge (only by projection
        %         %from cell-center till middle-point)
        %         swonedgelim = Sw(elemeval(1)) + dot(grad,rcm);
        
        
        %----------------------------------------------------------------------
        %Only the MLP is turned on (Concentratio recovery on the MIDEDGE -
        %TRADITIONAL application).
    elseif strcmp(cvbtype(1:3),'mlp') || ...
            strcmp([eblimkey cvblimkey],'offoff')
        %Get the gradient (reconstructed by MUSCL strategy) for each
        %element evaluated
        
        pos = 2*booleanelem + (1 - booleanelem);
        if strcmp(limiterflag{6},'on')|| strcmp(limiterflag{11},'on')
            
            
            %Calculate the gradients centered and upwind.
            %Calculate the vector distance between colocation points.
            dlr = pointr - pointl;            
            gradcent = Sw(elemeval(pos)) - Sw(elemeval(1));
            %Calculate the gradient upwind ("gradupwd"). First "grad"
            %matrix row.
            gradupwd = 2*dot(grad(1,:),dlr) - gradcent;    
            %Calculate the gradient limited according "k" value.
            gradslimited = mlpbyelem*0.5*((1 - k)*gradupwd + (1 + k)*gradcent);
            swonedgelim = Sw(elemeval(1)) + 0.5*gradslimited;
           
            if strcmp(limiterflag{11},'on')
                dlr = pointr - pointl;
                gradcent = Sw(elemeval(pos)) - Sw(elemeval(1));
                gradupwdnonlimit= 2*dot(gradnonlimit(1,:),dlr) - gradcent;
                gradsnonlimit = 0.5*((1 - k)*gradupwdnonlimit + (1 + k)*gradcent);
                swonedgenolim=Sw(elemeval(1)) + 0.5*gradsnonlimit;
            else
                swonedgenolim=0;
            end

        % For LP MUSCL    
        elseif strcmp(limiterflag{13},'on')
            % ################################################################
             %grad(1,:)=mlpbyelem*grad(1,:);
            if length(elemeval)>1
                %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                %Calculate the gradients centered and upwind.
                %Calculate the vector distance between colocation points.
                dlr = pointr - pointl;            
                gradcent = Sw(elemeval(pos)) - Sw(elemeval(1));
                %Calculate the gradient upwind ("gradupwd"). First "grad"
                %matrix row.
                gradupwd = 2*dot(grad(1,:),dlr) - gradcent;    
                %Calculate the gradient limited according "k" value.
                %gradslimited1 = 0.25*((1 - k)*gradupwd + (1 + k)*gradcent);
                
                %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                
                %dlr1 = midedge - pointl;
                vetor= 2*midedge-(pointr + pointl);
                
                %Swaux =Sw(elemeval(pos))+ dot(grad(2,:),vetor);
                %a = 0.5*(Swaux - Sw(elemeval(1)));
                aa=0.5*k*dot(grad(2,:),vetor)+0.5*(1-k)*dot(grad(1,:),vetor);
                %gradslimited =mlpbyelem(dot(grad(1,:),dlr1) + k*(a-dot(grad(1,:),dlr1)));
                 gradslimited =mlpbyelem*0.25*((1 - k)*gradupwd + (1 + k)*gradcent)+mlpbyelem*aa;
            else
                
                %Calculate the gradients centered and upwind.
                %Calculate the vector distance between colocation points.
                dlr1 = midedge - pointl;
                Swaux =Sw(elemeval(pos));
                a = 0.5*(Swaux - Sw(elemeval(1)));
                gradslimited =mlpbyelem*(k*a+(1-k)*dot(grad(1,:),dlr1)); % método de Burg (eq. 3)
                %gradslimited =mlpbyelem*(dot(grad(1,:),dlr1)); % método de Burg (eq. 3)
            end
            swonedgelim = Sw(elemeval(1)) + gradslimited;
            swonedgenolim=0;
        end
        % ################################################################
        %No one valid option was chosen. It gives a mensage for the user.
    else
        disp(' ');
        disp('---------------------------------------------------------------');
        disp('---------------------------------------------------------------');
        disp('>>> WARNING!');
        disp('---------------------------------------------------------------');
        disp('---------------------------------------------------------------');
        disp('There exists some problem with the limiter strategy!!!');
        disp('Check the settings into the Start.dat file.');
        disp('---------------------------------------------------------------');
    end  %End of IF (limiter type)
end  %End of IF (reconstruction type)

%--------------------------------------------------------------------------
%Function "get3rdand4thorder"
%--------------------------------------------------------------------------

function [swonedge] = get3rdand4thorder(order,elemeval,vertices,...
    verticescoord,taylorterms,Sw,constraint,mlpbyelem,centelemeval)
%Recovery the saturation through the edge evaluated. It uses two
%Gauss integration points.
%Initialize "phi"
phi = mlpbyelem;

%Initialize the saturation on the element evaluated
satval = Sw(elemeval(1));
%Define the grad of saturation (second and third "taylorterms" columns).
gradsat = taylorterms(elemeval(1),1:2);
%Define the Hessian matrix (3rd to fifth "taylorterms" columns):
hessianterms = taylorterms(elemeval(1),3:5);
%Fill the "hessian" matrix.
hessian = ...
    [hessianterms(1) hessianterms(2); hessianterms(2) hessianterms(3)];
%Define the terms of the constraint equation. See Clain et al. 2011
consteqterm = constraint(elemeval(1),1:5);

%Get the Gauss point coordinates:
%End-points coordinate:
a = verticescoord(1,1:2);
b = verticescoord(2,1:2);
%Point nearest the first vertex:
c1 = ((a + b)./2) - ((b - a)./(2.*sqrt(3)));
%Point nearest the second vertex:
c2 = ((a + b)./2) + ((b - a)./(2.*sqrt(3)));

%Initialize the Gauss weight
gaussw = 0.5;

%Define the Gauss weights for third-order approximation and the
%distance between the quadrature points and the cell-center.
%gaussw = 0.5;
%Get distances:
d1 = c1 - centelemeval(1,:);
d2 = c2 - centelemeval(1,:);

%Reconstruction of saturation on the quadrature points:
%A Third order scheme
if order == 3
    %Calculate the saturation value
    swonedge = gaussw*(satval + (phi(1)*dot(d1,gradsat) + ...
        phi(2)*((d1*hessian)*d1')/2) - ...
        (consteqterm*[phi(1)*gradsat phi(2)*[0.5 1 0.5].*hessianterms]')) ...
        + gaussw*(satval + (phi(1)*dot(d2,gradsat) + ...
        phi(2)*((d2*hessian)*d2')/2) - ...
        (consteqterm*[phi(1)*gradsat phi(2)*[0.5 1 0.5].*hessianterms]'));
    %A Fourth order scheme
else
    %Define the terms of the constraint equation. See Clain et al. 2011
    consteqterm = constraint(elemeval(1),:);
    
    %Define the third derivative terms (sixth to nineth "taylorterms"
    %columns):
    der3terms = taylorterms(elemeval(1),6:9);
    
    %Calculate the saturation value
    swonedge = gaussw*(satval + phi(1)*dot(d1,gradsat) + ...
        phi(2)*(((d1*hessian)*d1')/2 + ...
        phi(3)*((d1(1)^3)*der3terms(1) + 3*(d1(1)^2)*d1(2)*der3terms(2) + ...
        3*(d1(2)^2)*d1(1)*der3terms(3) + (d1(2)^3)*der3terms(4))/6) - ...
        (consteqterm*[phi(1)*gradsat phi(2)*[0.5 1 0.5].*hessianterms ...
        phi(2)*phi(3)*[(1/6) 0.5 0.5 (1/6)].*der3terms]')) + ...
        gaussw*(satval + phi(1)*dot(d2,gradsat) + ...
        phi(2)*(((d2*hessian)*d2')/2 + ...
        phi(3)*((d2(1)^3)*der3terms(1) + 3*(d2(1)^2)*d2(2)*der3terms(2) + ...
        3*(d2(2)^2)*d2(1)*der3terms(3) + (d2(2)^3)*der3terms(4))/6) - ...
        (consteqterm*[phi(1)*gradsat phi(2)*[0.5 1 0.5].*hessianterms ...
        phi(2)*phi(3)*[(1/6) 0.5 0.5 (1/6)].*der3terms]'));
    %     swonedge = [(satval + phi(1)*dot(d1,gradsat) + ...
    %         phi(2)*((d1*hessian)*d1')/2 + ...
    %         phi(3)*((d1(1)^3)*der3terms(1) + 3*(d1(1)^2)*d1(2)*der3terms(2) + ...
    %         3*(d1(2)^2)*d1(1)*der3terms(3) + (d1(2)^3)*der3terms(4))/6 - ...
    %         (consteqterm*[phi(1)*gradsat phi(2)*[0.5 1 0.5].*hessianterms ...
    %         phi(3)*[(1/6) 0.5 0.5 (1/6)].*der3terms]')) ...
    %         (satval + phi(1)*dot(d2,gradsat) + ...
    %         phi(2)*((d2*hessian)*d2')/2 + ...
    %         phi(3)*((d2(1)^3)*der3terms(1) + 3*(d2(1)^2)*d2(2)*der3terms(2) + ...
    %         3*(d2(2)^2)*d2(1)*der3terms(3) + (d2(2)^3)*der3terms(4))/6 - ...
    %         (consteqterm*[phi(1)*gradsat phi(2)*[0.5 1 0.5].*hessianterms ...
    %         phi(3)*[(1/6) 0.5 0.5 (1/6)].*der3terms]'))];
    
    %Calculate the saturation recovered on the edge
    %(mean of both quadrature points)
    %     swonedge = gaussw*sum(swonedge);
end  %End of IF

%--------------------------------------------------------------------------
%Function "get5thand6thorder"
%--------------------------------------------------------------------------

function [swonedge] = get5thand6thorder(order,elemeval,verticescoord,...
    taylorterms,Sw,constraint,mlpbyelem)
%Define global parameters:
global centelem;

%Initialize "phi"
phi = mlpbyelem(elemeval(1),:);
%Recovery the saturation through the edge evaluated. It uses two
%Gauss integration points.
%Define a first order term:

satval = Sw(elemeval(1));
%Define the grad of saturation (second and third "taylorterms" columns).
gradsat = taylorterms(elemeval(1),1:2);
%Define the Hessian matrix (3rd to fifth "taylorterms" columns):
hessianterms = taylorterms(elemeval(1),3:5);
%Fill the "hessian" matrix.
hessian = ...
    [hessianterms(1) hessianterms(2); hessianterms(2) hessianterms(3)];
%Define the third derivative terms (sixth to nineth "taylorterms" columns):
der3terms = taylorterms(elemeval(1),6:9);
%Define the fourth derivative terms (tenth to fourteenth "taylorterms"
%columns):
der4terms = taylorterms(elemeval(1),10:14);

%Define the terms of the constraint equation. See Clain et al. 2011
consteqterm = constraint(elemeval(1),1:14);

%Get the Gauss point coordinates:
%End-points coordinate:
a = verticescoord(1,1:2);
b = verticescoord(2,1:2);
%Define the Gauss geometric weight, "ggw"
ggw = 0.774596669241483;
%Initialize "c"
c = zeros(3,2);
%Point nearest the first vertex:
c(1,:) = ((a + b)./2) - ((b - a)./2).*ggw;
%Midway point:
c(2,:) = ((a + b)./2);
%Point nearest the second vertex:
c(3,:) = ((a + b)./2) + ((b - a)./2).*ggw;

%Define the Gauss weights for third-order approximation and the
%distance between the quadrature points and the cell-center.
gaussw = 0.5*[0.555555555555556 0.888888888888889 0.555555555555556];

%Reconstruction of saturation on the quadrature points:
%A Fifth order scheme
if order == 5
    %Initialize "swbyquadpoint"
    swbyquadpoint = zeros(1,3);
    %Calculate the saturation value
    for i = 1:3
        %Get the distance between the cell-center and the quadrature point.
        d = c(i,:) - centelem(elemeval(1),1:2);
        %Calculate the saturation on the quadrature point
        swbyquadpoint(i) = gaussw(i)*(satval + dot(d,phi(1)*gradsat) + ...
            phi(2)*((((d*hessian)*d')/2) + phi(3)*((((d(1)^3)*der3terms(1) + ...
            3*(d(1)^2)*d(2)*der3terms(2) + 3*(d(2)^2)*d(1)*der3terms(3) + ...
            (d(2)^3)*der3terms(4))/6) + phi(4)*(((d(1)^4)*der4terms(1) + ...
            4*(d(1)^3)*d(2)*der4terms(2) + 4*d(1)*(d(2)^3)*der4terms(4) + ...
            6*(d(1)^2)*(d(2)^2)*der4terms(3) + (d(2)^4)*der4terms(5))/24))) - ...
            (consteqterm*[phi(1)*gradsat phi(2)*[0.5 1 0.5].*hessianterms ...
            phi(2)*phi(3)*[(1/6) 0.5 0.5 (1/6)].*der3terms ...
            phi(2)*phi(3)*phi(4)*[(1/24) (1/6) (1/4) (1/6) (1/24)].*...
            der4terms]'));
    end  %End of FOR
    
    %Get the recovery value
    swonedge = sum(swbyquadpoint);
    
    %A Sixth order scheme
else
    %Define the terms of the constraint equation. See Clain et al. 2011
    consteqterm = constraint(elemeval(1),:);
    
    %Define the fifth derivative terms (fifteenth to twentyth "taylorterms"
    %columns):
    der5terms = taylorterms(elemeval(1),15:20);
    
    %Initialize "swbyquadpoint"
    swbyquadpoint = zeros(1,3);
    %Calculate the saturation value
    for i = 1:3
        %Get the distance between the cell-center and the quadrature point.
        d = c(i,:) - centelem(elemeval(1),1:2);
        %Calculate the saturation on the quadrature point
        swbyquadpoint(i) = gaussw(i)*(satval + dot(d,phi(1)*gradsat) + ...
            phi(2)*((((d*hessian)*d')/2) + phi(3)*((((d(1)^3)*der3terms(1) + ...
            3*(d(1)^2)*d(2)*der3terms(2) + 3*(d(2)^2)*d(1)*der3terms(3) + ...
            (d(2)^3)*der3terms(4))/6) + phi(4)*((((d(1)^4)*der4terms(1) + ...
            4*(d(1)^3)*d(2)*der4terms(2) + 4*d(1)*(d(2)^3)*der4terms(4) + ...
            6*(d(1)^2)*(d(2)^2)*der4terms(3) + (d(2)^4)*der4terms(5))/24) + ...
            phi(5)*(((d(1)^5)*der5terms(1) + 5*(d(1)^4)*d(2)*der5terms(2) + ...
            10*(d(1)^3)*(d(2)^2)*der5terms(3) + ...
            10*(d(1)^2)*(d(2)^3)*der5terms(4) + ...
            5*d(1)*(d(2)^4)*der5terms(5) + (d(2)^5)*der5terms(6))/120)))) - ...
            (consteqterm*[phi(1)*gradsat phi(2)*[0.5 1 0.5].*hessianterms ...
            phi(2)*phi(3)*[(1/6) 0.5 0.5 (1/6)].*der3terms ...
            phi(2)*phi(3)*phi(4)*[(1/24) (1/6) (1/4) (1/6) (1/24)].*...
            der4terms phi(2)*phi(3)*phi(4)*phi(5)*[(1/120) (1/24) (1/12) ...
            (1/12) (1/24) (1/120)].*der5terms]'));
    end  %End of FOR
    
    %Get the recovery value
    swonedge = sum(swbyquadpoint);
end  %End of IF

