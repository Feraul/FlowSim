%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: Element-based Limiter 
%Type of file: FUNCTION
%Criate date: 27/11/2015 (national workers strik)
%Modify data:   /  /2015
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Produce limiter value for 2nd to 6th order in order to guarantee 
%monotonic solution for saturation equation. 

%--------------------------------------------------------------------------
%

%--------------------------------------------------------------------------

function [mlplimiter] = MLPlimiter(auxvecforswept,Sw,taylorterms,cvbtype,...
    flagknownvert,satonvertices,constraint)
%Define global parameters:
global coord elem centelem order esurn1 esurn2;

%Define tolerance. It is like a computational zero.
tol = 1e-12;
%Initialize "coordsize" and "elemsize"
coordsize = size(coord,1);
elemsize = size(elem,1);
%Initialize MLP limiter
mlplimiter = ones(elemsize,order - 1);

if length(auxvecforswept) < elemsize
    getvertex = elem(auxvecforswept,1:4);
    vetxtoeval = unique(reshape(getvertex,4*size(getvertex,1),1));
    vetxtoeval = vetxtoeval(vetxtoeval ~= 0);
else
    vetxtoeval = 1:coordsize;
end  %End of IF

%Obtained from Park and Kim (2012). Look for MLP-u1.
%Swept all vertices
for ivtx = vetxtoeval
    vetx_coord = coord(ivtx,1:2);
    %Get the elements surrounding the evaluated vertex
    esurn = esurn1(esurn2(ivtx) + 1:esurn2(ivtx + 1));
    
    %Initialize "delta_sw"
    delta_sw = zeros(length(esurn),1);
    %Get the projected saturation
    proj_sw = Sw(esurn);
    %Get the gradient for the evaluated element
    localgrad = taylorterms(esurn,1:2);
    %Get the centroid coordinate
    localcentcoord = centelem(esurn,1:2);
    %get auxiliary "vert_coord"
    e = ones(length(esurn),1);
    rv = e;
    vert_coord_aux = [e*vetx_coord(1) e*vetx_coord(2)];
    %Define the vector "r"
    r_esrn = vert_coord_aux - localcentcoord;
    
    %Get max and min value of Sw (projected)
    
    if flagknownvert(ivtx) == 1 % Sw only contour
        extr_sw = dot(localgrad',r_esrn')';
        %extr_sw = satonvertices(ivtx)*e - proj_sw;
        %Get max and min value of Sw (projected)
       
        max_sw = max([satonvertices(ivtx); proj_sw]);
        min_sw = min([satonvertices(ivtx); proj_sw]);
        
    else
        %Get the extrapolated Sw
        extr_sw = dot(localgrad',r_esrn')';
        %Get max and min value of Sw (projected)
        max_sw = max(proj_sw);
        min_sw = min(proj_sw);
    end  %End of IF

    %Verify accuracy
    extr_sw = extr_sw.*(abs(extr_sw) > tol);
    
    isbigger = extr_sw > 0;
    islower = extr_sw < 0;
    nonzero = extr_sw ~= 0;
    delta_sw(isbigger) = max_sw*e(isbigger) - proj_sw(isbigger);
    delta_sw(islower) = min_sw*e(islower) - proj_sw(islower);
    
    rv(nonzero) = delta_sw(nonzero)./extr_sw(nonzero);

    %Calculate "qsi" (It Choces according flag "cvbtype")
    switch cvbtype
        %Original MLP (Park and Kin, 2012)
        case 'mlp'
            qsi = min(1,rv);
        %MLP-slc (Souza, Lyra, Carvalho)
        case 'mlp_vk'
            qsi = (rv.^2 + 2*rv)./(rv.^2 + rv + 2);
        %MLP-Van Albada
        case 'mlp_va'
            qsi = (rv.^2 + rv)./(rv.^2 + 1);
        %MLP-Van Albada 2 (Löhner, 2001)
        case 'mlp_va2'
            qsi = (2*rv)./(rv.^2 + 1);
        %MLP-Van Leer 1 (TVD second order for any value of "rv")
        case 'mlp_vl'
            qsi = 2*rv/(rv + 1);
        %MLP-Van Leer 2 (first-order for "rv" bigger than 1)
        case 'mlp_vl2'
            qsi = 4*rv./((rv + 1).^2);
        %MLP-New Limiter 2 (Mandal, 2008)
        case 'mlp_nl2'
            qsi = (rv.^3 + rv)./(rv.^3 + 1);
        %MLP-New Limiter 2 (Mandal, 2008)
        case 'mlp_nl3'
            qsi = (rv.^2 + 3*rv)./(rv.^2 + rv + 2);
    end  %End of SWITCH
    
    %Attribute "qsi" to the vector
    pointlower = qsi < mlplimiter(esurn);
    mlplimiter(esurn(pointlower)) = qsi(pointlower); 
end  %End of FOR
% %Define global parameters:
% global coord elem elemarea centelem order esurn1 esurn2;
% 
% %Define tolerance. It is like a computational zero.
% tol = 1e-12;
% %Initialize "coordsize" and "elemsize"
% coordsize = size(coord,1);
% elemsize = size(elem,1);
% %Initialize MLP limiter
% mlplimiter = ones(elemsize,order - 1);
% %Initialize some vectors that store the projected value for the order
% %evaluated and a flag to indicate if that vertex was evaluated.
% %It is the maximum value among all the finite volume mean value 
% maxmeanval = zeros(coordsize,1); 
% minmeanval = maxmeanval; 
% 
% %It is the maximum value among all the projected (3rd order) in each 
% %vertex 
% switch order 
%     %Third order
%     case 3
%         %Initialize variable that store the maximum and minimum
%         %projected value
%         maxval3order = maxmeanval;
%         minval3order = maxmeanval;
%     %Fourth order
%     case 4
%         %Initialize variable that store the maximum and minimum
%         %projected value
%         maxval4order = maxmeanval;
%         minval4order = maxmeanval;
% end  %End of IF
% 
% %It is a vector with the flag that indicate if the vertex evaluated 
% %was visited.
% flagwaseval = maxmeanval;
% 
% %Obtained from Park and Kim (2012). Look for MLP-u1.
% %Swept all elements
% for ivec = 1:length(auxvecforswept)
%     %Define the real counter
%     elemeval = auxvecforswept(ivec);
%     %Get the amount of vertices for the element evaluated
%     amountedges = sum(elem(elemeval,1:4) ~= 0);
%     %Initialize "qsi"
%     qsi = zeros(amountedges,1);
%     %Define the vertices
%     vertices = elem(elemeval,1:amountedges);
%     %Initialize "qsi3rd", "qsi4th", "qsi5th" and "qsi6th" 
%     qsi3rd = qsi;
%     qsi4th = qsi;
%     qsi5th = qsi;
%     qsi6th = qsi;
% 
%     %Saturation on the centroid
%     satoncent = Sw(elemeval);
% 
%     %Get some additional parameters if "mlp_vk" is turned on.
%     if strcmp(cvbtype,'mlp_vk')
%         %Define additional parameters (K, n, etc):
%         K = 100;
%         K1 = 0.01;
%         K2 = 0.01;
%         n = 1.5;
%         %Initialize "deltax"
%         deltax = zeros(amountedges,1);
%         %Swept all vertices
%         for i = 1:amountedges
%             %Put the evaluated vertex as the 1st of the vector "vertices". 
%             %It uses "shiftchosen" for that.
%             auxverti = shiftchoosen(vertices,vertices(i),'val');
%             %Get "deltax"
%             deltax(i) = ...
%                 norm(coord(auxverti(2),1:2) - coord(auxverti(1),1:2));
%         end  %End of FOR
%         %Get a mean value for "deltax"
%         deltax = mean(deltax);
%     end  %End of IF
%     
%     %Swept all vertices for the certain element.
%     %The number of edges is equal to number of vertices
%     for i = 1:amountedges
%         %Get the amount of elements surrounding the vertex evaluated
%         esurn = esurn1(esurn2(vertices(i)) + 1:esurn2(vertices(i) + 1));
%         %Reorder "esurn" in order the element evaluated "elemeval" came 
%         %first
%         esurn = shiftchoosen(esurn,elemeval,'val');
%         
%         %Get a little vector with the positions in "esurn"
%         amountesurn = 1:length(esurn);
%         
%         %Verify if the vertex was evaluated
%         %The vertex was NOT evaluated. MUST be calculated
%         if flagwaseval(vertices(i)) == 0
%             %Verify if the node belongs to "bedge" and if the proble has a
%             %known saturation as boundary condition.
%             if flagknownvert(vertices(i)) == 1
%                 %Define saturation on the vertex evaluated
%                 satonverteval = satonvertices(vertices(i));
%                 %Get the bigger saturatio value among those in the 
%                 %surrounding
%                 maxsatval = max([satonverteval; Sw(esurn)]);
%                 %Get the lower saturatio value among those in the 
%                 %surrounding
%                 minsatval = min([satonverteval; Sw(esurn)]);
%             %The vertex does not belong to known boundary 
%             else
%                 %Get the bigger saturatio value among those in the 
%                 %surrounding
%                 maxsatval = max(Sw(esurn));
%                 %Get the lower saturatio value among those in the 
%                 %surrounding
%                 minsatval = min(Sw(esurn));
%             end  %End of IF.
%             
%             %Verify the accuracy of SOME values.
%             %1. Volume average value:
%             maxsatval = maxsatval*(abs(maxsatval) > tol);
%             minsatval = minsatval*(abs(minsatval) > tol);
%             %Attribute the maximum and the minimum value to "maxmeanval" 
%             %and "minmeanval", respectively
%             maxmeanval(vertices(i)) = maxsatval; 
%             minmeanval(vertices(i)) = minsatval;
%             
%         %The vertex was evaluated
%         else
%             %"maxmeanval" and "minmeanval" receive the values already
%             %calculated.
%             maxsatval = maxmeanval(vertices(i)); 
%             minsatval = minmeanval(vertices(i));
%         end  %End of IF
%         
%         %Get the projection value on the vertex
%         vertexcoord = coord(vertices(i),1:2);
%         %Define the distance (centroid -> vertex)
%         distcv = vertexcoord - centelem(elemeval,1:2);
%         %Chose the type of prjection (second-order or third-order)
%         gettaylorsize = size(taylorterms,2);
%             
%         %------------------------------------------------------------------
%         %Evaluate according to "order":
%         
%         %------------
%         %FOURTH ORDER (grads, hessian and 3rd deriv - "taylorterms" has 
%         %nine columns)
%         if gettaylorsize >= 9
%             %Define "projval_4th" according the "flagwaseval" value
%             booleanvecsize = flagwaseval(vertices(i)) == 1;
%             localvecsize = length(esurn)*(1 - booleanvecsize) + ...
%                 booleanvecsize;
%             %Initialize "projval_4th"
%             projval_4th = zeros(localvecsize,1);
% 
%             %Swept all elements surrounding the vertex evaluated.
%             for isurnode = 1:localvecsize
%                 %Define "pointeresurn"
%                 pointeresurn = esurn(isurnode);
%                 
%                 %Get the gradient of saturation for the vertex evaluated:
%                 gradsat = taylorterms(pointeresurn,1:2);
%                 %Get hessian:
%                 hessianterms = taylorterms(pointeresurn,3:5);
%                 %Fill the "hessian" matrix.
%                 hessian = [hessianterms(1) hessianterms(2); ...
%                     hessianterms(2) hessianterms(3)];
%                 %Get third derivative:
%                 der3terms = taylorterms(pointeresurn,6:9);
%             
%                 %Define the terms of the constraint equation. 
%                 %See Clain et al. 2011
%                 consteqterm = constraint(pointeresurn,1:9); 
%         
%                 %Define the distance (centroid -> vertex)
%                 distcv4th = vertexcoord - centelem(pointeresurn,1:2);
% 
%                 %Get "pnprojslope" (see Park and Kim, 2012; ICCFD7 - 4204)
%                 pnprojslope = dot(gradsat,distcv4th);% - ...
%                     %consteqterm(1:2)*gradsat';
%                 %Get "p1filtered" (see Park and Kim, 2012; ICCFD7 - 4204)
%                 p1filtered = ((distcv4th*hessian)*distcv4th'/2) + ...
%                     (((distcv4th(1)^3)*der3terms(1) + ...
%                     3*(distcv4th(1)^2)*distcv4th(2)*der3terms(2) + ...
%                     3*(distcv4th(2)^2)*distcv4th(1)*der3terms(3) + ...
%                     (distcv4th(2)^3)*der3terms(4))/6);% - ...
%                     %(consteqterm(3:9)*[hessianterms der3terms]');
%                 
%                 %Projection of Fourth order on the vertex evaluated.
%                 projval_4th(isurnode) = Sw(pointeresurn) + pnprojslope + ...
%                     p1filtered;
% 
%                 %Attribute "pnprojslope" and "p1filtered" to 
%                 %"pnprojslope_4th" and "p1filtered_4th"
%                 if pointeresurn == elemeval 
%                     pnprojslope_4th = pnprojslope;
%                     p1filtered_4th = p1filtered;
%                     projval_onvtx4th = projval_4th(isurnode);
%                 end  %End of IF
%             end  %End of FOR (swept "esurn")
%                 
%             %Verify if the vertex was evaluated
%             %The vertex was NOT evaluated. MUST be calculated
%             if flagwaseval(vertices(i)) == 0
%                 %Verify if the node belongs to "bedge" and if the proble 
%                 %has a known saturation as boundary condition.
%                 if flagknownvert(vertices(i)) == 1
%                     %Define saturation on the vertex evaluated
%                     satonverteval = satonvertices(vertices(i));
%                     %Get the bigger saturatio value among those in the sur.
%                     maxprojval = max([satonverteval; projval_4th]);
%                     %Get the lower saturatio value among those in the sur.
%                     minprojval = min([satonverteval; projval_4th]);
%                     %Get "projval" from boundary value
%                     projval_onvtx4th = satonvertices(vertices(i));
%                 %The vertex does not belong to known boundary 
%                 else
%                     %Get the bigger saturatio value among those in the sur.
%                     maxprojval = max(projval_4th);
%                     %Get the lower saturatio value among those in the sur.
%                     minprojval = min(projval_4th);
%                 end  %End of IF.
%                 
%                 %Evaluate the parameter "qsi4th"
%                 %1. Projection value:
%                 maxprojval = maxprojval*(abs(maxprojval) > tol);
%                 minprojval = minprojval*(abs(minprojval) > tol);
%             
%                 %Attribute the maximum and the minimum projected value to 
%                 %"maxval3order" and "minval3order", respectively
%                 maxval4order(vertices(i)) = maxprojval; 
%                 minval4order(vertices(i)) = minprojval;
%             
%             %The vertex was evaluated (4th order)
%             else
%                 %"maxprojval" and "minprojval" receive the values already
%                 %calculated.
%                 maxprojval = maxval4order(vertices(i)); 
%                 minprojval = minval4order(vertices(i));
%             end  %End of IF
%                 
%             %Verify the accuracy of all values.
%             %2. The projected value (4th order) value:
%             projval_onvtx4th = projval_onvtx4th*(abs(projval_onvtx4th) > tol);
%             %3. 2nd order projection ("pnprojslope_4th"):
%             pnprojslope_4th = pnprojslope_4th*(abs(pnprojslope_4th) > tol);
%             %4. very higher order projection ("pnprojslope_4th"):
%             p1filtered_4th = p1filtered_4th*(abs(p1filtered_4th) > tol);
% 
%             %Define augmented MLP criteria:
%             qsi4th(i) = ((projval_onvtx4th >= minprojval && ...
%                 minprojval >= minsatval) && ...
%                 (projval_onvtx4th <= maxprojval && ...
%                 maxprojval <= maxsatval));
%             
%             %Verify other condition
%             if qsi4th(i) == 0
%                 
%                 %###################################################
% %                 boolean4th = projval_onvtx4th - satoncent > 0; 
% %                 delta_4th = (maxsatval - satoncent)*boolean4th + ...
% %                     (minsatval - satoncent)*(1 - boolean4th);
% %                 qsi4th(i) = ...
% %                     min(1,delta_4th/(projval_onvtx4th - satoncent)); 
%                 %###################################################
%                 
%                 %Deactivation Condition (for smooth regions)
%                 qsi4th(i) = (abs((projval_onvtx4th - satoncent)/...
%                     (satoncent + 1e-16)) <= max(1e-3,elemarea(elemeval)));
% 
% %                 qsi4th(i) = ...
% %                     ((sign(pnprojslope_4th) ~= sign(p1filtered_4th) && ...
% %                     (projval > minsatval && projval < maxsatval)) || ...
% %                     (abs((projval - satoncent)/(satoncent + 1e-16)) <= ...
% %                     1e-3));
%             end  %End of IF
%             
%             %Attribute the limiter value of 4th order to 3rd order limiter.
%             %It is used only for the conventional MLP
%             %qsi3rd(i) = qsi4th(i);
%         end  %End of IF (fourth order)
%         
%         %-----------
%         %THIRD ORDER (grads and hessian - "taylorterms" has five columns)
%         if gettaylorsize >= 5 %&& qsi3rd(i) ~= 1
%             %Define "projvalvho" according the "flagwaseval" value
%             booleanvecsize = flagwaseval(vertices(i)) == 1;
%             localvecsize = length(esurn)*(1 - booleanvecsize) + ...
%                 booleanvecsize;
%             %Initialize "projvalvho"
%             projvalvho = zeros(localvecsize,1);
%             %Swept all elements surrounding the vertex evaluated.
%             for isurnode = 1:localvecsize
%                 posinesurn = (amountesurn(logical(esurn == elemeval)))*...
%                     booleanvecsize + isurnode*(1 - booleanvecsize);
%                 %Define "pointeresurn"
%                 pointeresurn = esurn(posinesurn);
% 
%                 %Get the gradient of saturation:
%                 gradsat = taylorterms(pointeresurn,1:2);
%                 %Get hessian for each elament of "esurn"
%                 %Define the Hessian matrix (3rd to 5th "taylorterms" 
%                 %columns):
%                 hessianterms = taylorterms(pointeresurn,3:5);
%                 %Fill the "hessian" matrix.
%                 hessian = [hessianterms(1) hessianterms(2); ...
%                     hessianterms(2) hessianterms(3)];
%                 %Define the terms of the constraint equation. 
%                 %See Clain et al. 2011 (Eq. 11).
%                 consteqterm = constraint(pointeresurn,1:5); 
% 
%                 %Define the distance (centroid -> vertex)
%                 distcvho = vertexcoord - centelem(pointeresurn,1:2);
% 
%                 %Get "pnprojslope" (see Park and Kim, 2012; ICCFD7 - 4204)
%                 pnprojslope = dot(gradsat,distcvho);% - ...
%                     %consteqterm(1:2)*gradsat';  
%                 %Get "p1filtered" (see Park and Kim, 2012; ICCFD7 - 4204)
%                 p1filtered = (((distcvho*hessian)*distcvho')/2);% - ...
%                    %(consteqterm(3:5)*hessianterms');
%                 
%                 %Define the position for store the value in "" 
%                 posinvho = booleanvecsize + isurnode*(1 - booleanvecsize);
%                 %Projection of third order on the vertex evaluated.
%                 projvalvho(posinvho) = Sw(pointeresurn) + ...
%                     pnprojslope + p1filtered; 
% 
%                 %Attribute "pnprojslope" and "p1filtered" to 
%                 %"pnprojslope_3rd" and "p1filtered_3rd"
%                 if pointeresurn == elemeval
%                     pnprojslope_3rd = pnprojslope;
%                     p1filtered_3rd = p1filtered;
%                 end  %End of IF
%             end  %End of FOR ("esurn")
% 
%             %--------------------------------------------------------------
%             %Get the maximum and minimum projections on the evaluated
%             %vertex:
%             
%             %Verify if the vertex was evaluated
%             %The vertex was NOT evaluated. MUST be calculated
%             if flagwaseval(vertices(i)) == 0
%                 %Verify if the node belongs to "bedge" and if the proble 
%                 %has a known saturation as boundary condition.
%                 if flagknownvert(vertices(i)) == 1
%                     %Define saturation on the vertex evaluated
%                     satonverteval = satonvertices(vertices(i));
%                     %Get the bigger saturatio value among those in the sur.
%                     maxprojval = max([satonverteval; projvalvho]);
%                     %Get the lower saturatio value among those in the sur.
%                     minprojval = min([satonverteval; projvalvho]);
%                     %Get "projval" from boundary value
%                     projval_onvtx3rd = satonvertices(vertices(i));
%                 %The vertex does not belong to known boundary 
%                 else
%                     %Get the bigger saturatio value among those in the sur.
%                     maxprojval = max(projvalvho);
%                     %Get the lower saturatio value among those in the sur.
%                     minprojval = min(projvalvho);
%                     %Get the projection for the "elemeval"
%                     projval_onvtx3rd = ...
%                         projvalvho(logical(esurn == elemeval));
%                 end  %End of IF.
%             
%                 %Evaluate the parameter "qsi3rd"
%                 %1. Projection value:
%                 maxprojval = maxprojval*(abs(maxprojval) > tol);
%                 minprojval = minprojval*(abs(minprojval) > tol);
%                 
%                 %Attribute the maximum and the minimum projected value to 
%                 %"maxval3order" and "minval3order", respectively
%                 maxval3order(vertices(i)) = maxprojval; 
%                 minval3order(vertices(i)) = minprojval;
%             
%             %The vertex was evaluated
%             else
%                 %"maxprojval" and "minprojval" receive the values already
%                 %calculated.
%                 maxprojval = maxval3order(vertices(i)); 
%                 minprojval = minval3order(vertices(i));
%             end  %End of IF
%             
%             %2. Projected value on the vertex evaluated (from the element 
%             %evaluated)
%             projval_onvtx3rd = ...
%                 projval_onvtx3rd*(abs(projval_onvtx3rd) > tol);
%             %3. 2nd order projection ("pnprojslope_4th"):
%             pnprojslope_3rd = pnprojslope_3rd*(abs(pnprojslope_3rd) > tol);
%             %4. very higher order projection ("pnprojslope_4th"):
%             p1filtered_3rd = p1filtered_3rd*(abs(p1filtered_3rd) > tol);
%             
%             %Define augmented MLP criteria:
%             qsi3rd(i) = ((projval_onvtx3rd >= minprojval && ...
%                 minprojval >= minsatval) && ...
%                 (projval_onvtx3rd <= maxprojval && ...
%                 maxprojval <= maxsatval));
% 
%             %Verify other condition
%             if qsi3rd(i) == 0  
% 
%                 %###################################################
% %                 boolean3rd = projval_onvtx3rd - satoncent > 0; 
% %                 delta_3rd = (maxsatval - satoncent)*boolean3rd + ...
% %                     (minsatval - satoncent)*(1 - boolean3rd);
% %                 qsi3rd(i) = ...
% %                     min(1,delta_3rd/(projval_onvtx3rd - satoncent)); 
%                 %###################################################
%                 
%                 %Deactivation Condition (for smooth regions)
%                 qsi3rd(i) = (abs((projval_onvtx3rd - satoncent)/...
%                     (satoncent + 1e-16)) <= max(1e-3,elemarea(elemeval)));
% %                 qsi3rd(i) = ...
% %                     ((sign(pnprojslope_3rd) ~= sign(p1filtered_3rd) && ...
% %                     (projval > minsatval && projval < maxsatval)) || ...
% %                     (abs((projval - satoncent)/(satoncent + 1e-16)) <= ...
% %                     1e-3));
%             end  %End of IF
%         end  %End of IF (third order)
%     
%         %------------
%         %SECOND ORDER (only grads - "grads" has two columns)
% 
% %         if gettaylorsize >= 5 && (qsi3rd(i) ~= 1 || qsi4th(i) ~= 1)
% %             %Call "get2ndorderecovery.m"
% %             grad = get2ndorderecovery(elemeval,Sw,flagknownedge,...
% %                 satonboundedges,bedgrownum);  
% %         %The gradiente is taken for a preprocessing way.
% %         else
%            %Get the gradient from "taylorterms"
%            grad = taylorterms(elemeval,1:2);
% %         end  %End of IF 
%             
%         %Verify if the saturation on the vertex evaluated is known
%         if flagknownvert(vertices(i)) == 1
%             projval_onvtx2nd = satonvertices(vertices(i)) - satoncent;
%         %The vertex does not belong to known boundary 
%         else
%             %Projection on the vertex
%             projval_onvtx2nd = dot(grad,distcv);
%         end  %End of IF
% 
%         %Verify if there is a null projection and the very higher order is
%         %turned off.
%         if (abs(projval_onvtx2nd) > tol && gettaylorsize <= 5) || ...
%                 (gettaylorsize >= 5 && qsi3rd(i) ~= 1 && ...
%                 abs(projval_onvtx2nd) > tol) || (gettaylorsize >= 9 && ...
%                 qsi4th(i) ~= 1 && abs(projval_onvtx2nd) > tol)
%             %The gradient is calculated right now in order to get a more
%             %accurated second order approximation.
% 
%             %Get the ratio "rv":
%             %Define a boolean condition
%             booleanrv = projval_onvtx2nd > 0;
%             %Defien "deltasatval"
%             deltasatval = (maxsatval - satoncent)*booleanrv + ...
%                 (1 - booleanrv)*(minsatval - satoncent);
%             
%             %Get the ratio "rv":
%             rv = deltasatval/projval_onvtx2nd;
% 
%             %Calculate "qsi" (It Choces according flag "cvbtype")
%             switch cvbtype
%                 %Original MLP (Park and Kin, 2012)
%                 case 'mlp'
%                     qsi(i) = min(1,rv);
%                 %MLP u2, Venkatakrishnan (see Park and Kin, 2012)
%                 case 'mlp_vk'
%                     %Define "deltaplus" according "projval" sign.
%                     booleanrv = projval_onvtx2nd > 0;
%                     %Defien "deltaplus"
%                     deltaplus = (maxsatval - satoncent)*booleanrv + ...
%                         (1 - booleanrv)*(minsatval - satoncent);
%                     
%                     %Define "deltaless"
%                     deltaless = projval_onvtx2nd;
%                     %Delta sw (difference between the max and min 
%                     %saturation values):
%                     deltasonvertex = maxsatval - minsatval;
%                     %Calculate "teta"
%                     teta = deltasonvertex/(K2*(deltax^n));
%                     %Calculate epsilon^2:
%                     epsilonsqre = (K1*(deltasonvertex^2))/(1 + teta);  %K*(deltax^3);%    
%                     
%                     %Use Venkathakrishnen approximation
%                     qsi(i) = (1/deltaless)*((((deltaplus^2) + epsilonsqre)*...
%                         deltaless + 2*(deltaless^2)*deltaplus)/...
%                         ((deltaplus^2) + 2*(deltaless^2) + ...
%                         deltaless*deltaplus + epsilonsqre));
%                 %MLP-slc (Souza, Lyra, Carvalho)
%                 case 'mlp_slc'
% %                     qsi(i) = min(1,(rv^2 + 2*rv)/(rv^2 + (rv^1) + 2));
%                     %qsi(i) = (rv^3 + 2*rv)/(rv^3 + (rv^2) + 3);
%                     qsi(i) = (rv^3 + 3*rv)/(rv^3 + (rv^2) + 4);
%                 %MLP-Van Albada
%                 case 'mlp_va'
%                     qsi(i) = (rv^2 + rv)/(rv^2 + 1);
%                 %MLP-Van Albada 2 (Löhner, 2001)
%                 case 'mlp_va2'
%                     qsi(i) = (2*rv)/(rv^2 + 1);
%                 %MLP-Van Leer 1 (TVD second order for any value of "rv")
%                 case 'mlp_vl'
%                     qsi(i) = 2*rv/(rv + 1);
%                 %MLP-Van Leer 2 (first-order for "rv" bigger than 1)
%                 case 'mlp_vl2'
%                     qsi(i) = 4*rv/((rv + 1)^2);
%                 %MLP-New Limiter 2 (Mandal, 2008)
%                 case 'mlp_nl2'
%                     qsi(i) = (rv^3 + rv)/(rv^3 + 1);
%                 %MLP-New Limiter 2 (Mandal, 2008)
%                 case 'mlp_nl3'
%                     qsi(i) = (rv^2 + 3*rv)/(rv^2 + rv + 2);
%             end  %End of SWITCH
%         
%         %There is a null projection
%         else
%             qsi(i) = 1;
%         end  %End of IF
%         
%         %Update "flagwaseval"
%         flagwaseval(vertices(i)) = 1;
%     end  %End of FOR (swept the vertices)
% 
%     %Define the limiter ("phi")
%     %Choose according to size of "taylorterms"
%     %Second order:
%     if gettaylorsize == 2
%         phi = min(qsi);
%     %Third order limited
%     elseif gettaylorsize == 5
%         phi = [min(qsi) min(qsi3rd)];
%     %Fourth order limited
%     elseif gettaylorsize == 9
%         phi = [min(qsi) min(qsi3rd) min(qsi4th)];        
%     %Fifth order limited
%     elseif gettaylorsize == 14
%         phi = [min(qsi) min(qsi3rd) min(qsi4th) min(qsi5th)];        
%     %Sixth order UNlimited
%     elseif gettaylorsize == 20 
%         phi = [min(qsi) min(qsi3rd) min(qsi4th) min(qsi5th) min(qsi6th)]; 
%     end  %End of IF (MLP limiter)
% 
%     %Fill "mlplimiter"
%     mlplimiter(elemeval,:) = phi; 
% end  %End of FOR (for all domain)

%========================================================================

% line = [4 12 20 28 36 44 52 60 68 76]; 
%     %84 92 100 108 116 124 132 140 148 156];
%     if timelevel == 132 && ismember(elemeval,line)
%         elemeval
%         phi
%     end



%--------------------------------------------------------------------------
%Garbage

%         %-----------
%         %SIXTH ORDER (grads, hessian, 3rd, 4th, 5th deriv - "taylorterms" 
%         %has twenty columns)
%         if gettaylorsize == 20
%             %Get the gradient of saturation for the vertex evaluated:
%             gradsat = taylorterms(elemeval,1:2);
%             %Get hessian:
%             hessianterms = taylorterms(elemeval,3:5);
%             %Fill the "hessian" matrix.
%             hessian = [hessianterms(1) hessianterms(2); ...
%                 hessianterms(2) hessianterms(3)];
%             %Get third derivative:
%             der3terms = taylorterms(elemeval,6:9);
%             %Get the fourth derivative:
%             der4terms = taylorterms(elemeval,10:14);
%             %Get the fifth derivative:
%             der5terms = taylorterms(elemeval,15:20);
%             
%             %Define the terms of the constraint equation. 
%             %See Clain et al. 2011
%             consteqterm = constraint(elemeval,:); 
%         
%             %Get "pnprojslope" (see Park and Kim, 2012; ICCFD7 - 4204)
%             pnprojslope = dot(gradsat,distcv) - ...
%                 consteqterm(1:2)*gradsat';
%             %Get "p1filtered" (see Park and Kim, 2012; ICCFD7 - 4204)
%             p1filtered = ((distcv*hessian)*distcv'/2) + ...
%                 (((distcv(1)^3)*der3terms(1) + ...
%                 3*(distcv(1)^2)*distcv(2)*der3terms(2) + ...
%                 3*(distcv(2)^2)*distcv(1)*der3terms(3) + ...
%                 (distcv(2)^3)*der3terms(4))/6 + ...
%                 (((distcv(1)^4)*der4terms(1) + ...
%                 4*(distcv(1)^3)*distcv(2)*der4terms(2) + ...
%                 4*distcv(1)*(distcv(2)^3)*der4terms(4) + ...
%                 6*(distcv(1)^2)*(distcv(2)^2)*der4terms(3) + ...
%                 (distcv(2)^4)*der4terms(5))/24) + ...
%                 (((distcv(1)^5)*der5terms(1) + ...
%                 5*(distcv(1)^4)*distcv(2)*der5terms(2) + ...
%                 10*(distcv(1)^3)*(distcv(2)^2)*der5terms(3) + ...
%                 10*(distcv(1)^2)*(distcv(2)^3)*der5terms(4) + ...
%                 5*distcv(1)*(distcv(2)^4)*der5terms(5) + ...
%                 (distcv(2)^5)*der5terms(6))/120)) - ...
%                 (consteqterm(3:20)*[[0.5 1 0.5].*hessianterms ...
%                 [(1/6) 0.5 0.5 (1/6)].*der3terms ...
%                 [(1/24) (1/6) (1/4) (1/6) (1/24)].*der4terms ...
%                 [(1/120) (1/24) (1/12) (1/12) (1/24) (1/120)].*der5terms]');
% 
%             %Verify if the node belongs to "bedge" and if the proble has a
%             %known saturation as boundary condition.
%             if flagknownvert(vertices(i)) == 1
%                 %Get "projval" from boundary value
%                 projval = satonvertices(vertices(i));
%             %The node does not belong to boundary. "projval" is calculated.
%             else
%                 %Get the projected value untill the veretex evaluated.
%                 projval = satoncent + pnprojslope + p1filtered;
%             end  %End of IF
% 
%             %Verify the accuracy of all values.
%             %1. The projected value (4th order) value:
%             projval = projval*(abs(projval) > tol);
%             %2. 2nd order projection ("pnprojslope"):
%             pnprojslope = pnprojslope*(abs(pnprojslope) > tol);
%             %3. very higher order projection ("pnprojslope"):
%             p1filtered = p1filtered*(abs(p1filtered) > tol);
% 
%             %Define augmented MLP criteria:
%             qsi6th(i) = 1 - ((pnprojslope > 0 && p1filtered < 0 && ...
%                 projval > minsatval) && (pnprojslope < 0 && ...
%                 p1filtered > 0 && projval < maxsatval) && ...
%                 abs((projval - satoncent)/(satoncent + 1e-16)) <= 1e-3);
%         end  %End of IF (sixth order)
% 
%         %-----------
%         %FIFTH ORDER (grads, hessian, 3rd and 4th deriv - "taylorterms" has 
%         %fourteen columns)
%         if gettaylorsize >= 14
%             %Initialize "projval_4th"
%             projval_5th = zeros(length(esurn),1);
%             %Swept all elements surrounding the vertex evaluated.
%             for isurnode = 1:length(esurn)
%                 %Define "pointeresurn"
%                 pointeresurn = esurn(isurnode);
%                 
%                 %Get the gradient of saturation for the vertex evaluated:
%                 gradsat = taylorterms(pointeresurn,1:2);
%                 %Get hessian:
%                 hessianterms = taylorterms(pointeresurn,3:5);
%                 %Fill the "hessian" matrix.
%                 hessian = [hessianterms(1) hessianterms(2); ...
%                     hessianterms(2) hessianterms(3)];
%                 %Get third derivative:
%                 der3terms = taylorterms(pointeresurn,6:9);
%                 %Get the fourth derivative:
%                 der4terms = taylorterms(pointeresurn,10:14);
%             
%                 %Define the terms of the constraint equation. 
%                 %See Clain et al. 2011
%                 consteqterm = constraint(pointeresurn,1:14); 
%         
%                 %Define the distance (centroid -> vertex)
%                 distcv5th = vertexcoord - centelem(pointeresurn,1:2);
% 
%                 %Get "pnprojslope" (see Park and Kim, 2012; ICCFD7 - 4204)
%                 pnprojslope = dot(gradsat,distcv5th) - ...
%                     consteqterm(1:2)*gradsat';
%                 %Get "p1filtered" (see Park and Kim, 2012; ICCFD7 - 4204)
%                 p1filtered = ((distcv5th*hessian)*distcv5th'/2) + ...
%                 (((distcv5th(1)^3)*der3terms(1) + ...
%                 3*(distcv5th(1)^2)*distcv5th(2)*der3terms(2) + ...
%                 3*(distcv5th(2)^2)*distcv5th(1)*der3terms(3) + ...
%                 (distcv5th(2)^3)*der3terms(4))/6 + ...
%                 (((distcv5th(1)^4)*der4terms(1) + ...
%                 4*(distcv5th(1)^3)*distcv5th(2)*der4terms(2) + ...
%                 4*distcv5th(1)*(distcv5th(2)^3)*der4terms(4) + ...
%                 6*(distcv5th(1)^2)*(distcv5th(2)^2)*der4terms(3) + ...
%                 (distcv5th(2)^4)*der4terms(5))/24)) - ...
%                 (consteqterm(3:14)*[[0.5 1 0.5].*hessianterms ...
%                 [(1/6) 0.5 0.5 (1/6)].*der3terms ...
%                 [(1/24) (1/6) (1/4) (1/6) (1/24)].*der4terms]');
%                 
%                 %Projection of Fourth order on the vertex evaluated.
%                 projval_5th(isurnode) = Sw(pointeresurn) + pnprojslope + ...
%                     p1filtered;
% 
%                 %Attribute "pnprojslope" and "p1filtered" to 
%                 %"pnprojslope_5th" and "p1filtered_5th"
%                 if pointeresurn == elemeval
%                     pnprojslope_5th = pnprojslope;
%                     p1filtered_5th = p1filtered;
%                     projval = projval_5th(isurnode);
%                 end  %End of IF
%             end  %End of FOR
%                 
%             %Verify if the node belongs to "bedge" and if the proble has a
%             %known saturation as boundary condition.
%             if flagknownvert(vertices(i)) == 1
%                 %Define saturation on the vertex evaluated
%                 satonverteval = satonvertices(vertices(i));
%                 %Get the bigger saturatio value among those in the sur.
%                 maxprojval = max([satonverteval; projval_5th]);
%                 %Get the lower saturatio value among those in the sur.
%                 minprojval = min([satonverteval; projval_5th]);
%                 %Get "projval" from boundary value
%                 projval = satonvertices(vertices(i));
%             %The vertex does not belong to known boundary 
%             else
%                 %Get the bigger saturatio value among those in the sur.
%                 maxprojval = max(projval_5th);
%                 %Get the lower saturatio value among those in the sur.
%                 minprojval = min(projval_5th);
%             end  %End of IF.
%             
%             %Evaluate the parameter "qsi4th"
%             %1. Projection value:
%             maxprojval = maxprojval*(abs(maxprojval) > tol);
%             minprojval = minprojval*(abs(minprojval) > tol);
%             %Verify the accuracy of all values.
%             %2. The projected value (4th order) value:
%             projval = projval*(abs(projval) > tol);
%             %3. 2nd order projection ("pnprojslope_4th"):
%             pnprojslope_5th = pnprojslope_5th*(abs(pnprojslope_5th) > tol);
%             %4. very higher order projection ("pnprojslope_4th"):
%             p1filtered_5th = p1filtered_5th*(abs(p1filtered_5th) > tol);
% 
%             %Define augmented MLP criteria:
%             qsi5th(i) = ...
%                 ((projval >= minprojval && minprojval >= minsatval) && ...
%                 (projval <= maxprojval && maxprojval <= maxsatval));
%             %Verify other condition
%             if qsi5th(i) == 0
%                 qsi5th(i) = ...
%                     (sign(pnprojslope_5th) ~= sign(p1filtered_5th));
%             end  %End of IF
%         end  %End of IF (fifth order)
