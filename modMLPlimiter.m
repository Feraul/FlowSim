function[mlplimiter] = modMLPlimiter(Sw,taylorterms,flagknownvert,...
    satonvertices)
%Define global parameters
global coord bedge inedge centelem elem esurn1 esurn2

%Define a tolerance
tol = 1e-12;
%Initialize "mlplimiter"
mlplimiter = ones(size(elem,1),1);
%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
coordsize = size(coord,1);
%Max and min values for each vertex
% maxvtxval = zeros(coordsize,1); 
% minvtxval = maxmeanval; 

%Swept "bedge":
for ibedg = 1:bedgesize
    %Initialize
    maxbyvtx = zeros(2,1);
    minbyvtx = maxbyvtx;
    %Get the vertices
    vertices = bedge(ibedg,1:2);
    elemeval = bedge(ibedg,3);
    
    %Define the grad
    grad = taylorterms(elemeval,1:2);
    %Get the distance between the centroid and the midpoint
    distcv = mean(coord(vertices,1:2)) - centelem(elemeval,1:2);

    for ivtx = 1:2
        esurn = ...
            esurn1(esurn2(vertices(ivtx)) + 1:esurn2(vertices(ivtx) + 1));
        %Verify if the node has a known saturation as boundary condi.
        if flagknownvert(vertices(ivtx)) == 1
            %Define saturation on the vertex evaluated
            satonverteval = satonvertices(vertices(ivtx));
            %Get the bigger saturatio value among those in the surrounding
            maxsatval = max([satonverteval; Sw(esurn)]);
            %Get the lower saturatio value among those in the surrounding
            minsatval = min([satonverteval; Sw(esurn)]);
            %Get the saturation extrapolated on the vertex
            sw_onmid = satonvertices(vertices(ivtx)) - Sw(elemeval);

        %The vertex does not have a known boundary 
        else
            %Get the bigger saturatio value among those in the surrounding
            maxsatval = max(Sw(esurn));
            %Get the lower saturatio value among those in the surrounding
            minsatval = min(Sw(esurn));
            %Get the saturation extrapolated on the vertex
            sw_onmid = dot(grad,distcv);
        end  %End of IF.
        
        maxbyvtx(ivtx) = maxsatval;
        minbyvtx(ivtx) = minsatval;
    end  %End of FOR
    
    %Verify if there is a null projection and the very higher order is
    %turned off.
    if abs(sw_onmid) > tol
        %Get the ratio "rv":
        %Define a boolean condition
        booleanrv = sw_onmid > 0;
        %Defien "deltasatval"
        deltasatval = (min(maxbyvtx) - Sw(elemeval))*booleanrv + ...
            (1 - booleanrv)*(max(minbyvtx) - Sw(elemeval));
            
        %Get the ratio "rv":
        rv = deltasatval/sw_onmid;

        %Calculate "qsi" (It Choces according flag "cvbtype")
        qsi = (2*rv)/(rv^2 + 1);
%         qsi = min(1,rv);
            
    %There is a null projection
    else
        qsi = 1;
    end  %End of IF
    
    %Verify the final limiter (on the element)
    if qsi < mlplimiter(elemeval)
        mlplimiter(elemeval) = qsi;
    end  %End of IF
end  %End of FOR ("bedge")

%Swept "inedge"
for iinedg = 1:inedgesize
    %Initialize
    maxbyvtx = zeros(2,1);
    minbyvtx = maxbyvtx;
    %Get the vertices
    vertices = inedge(iinedg,1:2);
    leftelem = inedge(iinedg,3);
    rightelem = inedge(iinedg,4);
    
    %Define the grad
    grad_left = taylorterms(leftelem,1:2);
    grad_right = taylorterms(rightelem,1:2);
    %Get the distance between the centroid and the midpoint
    distcv_left = mean(coord(vertices,1:2)) - centelem(leftelem,1:2);
    distcv_right = mean(coord(vertices,1:2)) - centelem(rightelem,1:2);

    for ivtx = 1:2
        esurn = ...
            esurn1(esurn2(vertices(ivtx)) + 1:esurn2(vertices(ivtx) + 1));
        %Verify if the node has a known saturation as boundary condi.
        if flagknownvert(vertices(ivtx)) == 1
            %Define saturation on the vertex evaluated
            satonverteval = satonvertices(vertices(ivtx));
            %Get the bigger saturatio value among those in the surrounding
            maxsatval = max([satonverteval; Sw(esurn)]);
            %Get the lower saturatio value among those in the surrounding
            minsatval = min([satonverteval; Sw(esurn)]);
            %Get the saturation extrapolated on the vertex
            sw_onmid_left = satonvertices(vertices(ivtx)) - Sw(leftelem);
            %Get the saturation extrapolated on the vertex
            sw_onmid_right = satonvertices(vertices(ivtx)) - Sw(rightelem);

        %The vertex does not have a known boundary 
        else
            %Get the bigger saturatio value among those in the surrounding
            maxsatval = max(Sw(esurn));
            %Get the lower saturatio value among those in the surrounding
            minsatval = min(Sw(esurn));
            %Get the saturation extrapolated on the vertex
            sw_onmid_left = dot(grad_left,distcv_left);
            %Get the saturation extrapolated on the vertex
            sw_onmid_right = dot(grad_right,distcv_right);
        end  %End of IF.
        
        maxbyvtx(ivtx) = maxsatval;
        minbyvtx(ivtx) = minsatval;
    end  %End of FOR
    
    %-------------------
    %ELEMENT on the LEFT

    %Verify if there is a null projection and the very higher order is
    %turned off.
    %For the control volume on the left
    if abs(sw_onmid_left) > tol
        %Get the ratio "rv":
        %Define a boolean condition
        booleanrv = sw_onmid_left > 0;
        %Defien "deltasatval"
        deltasatval = (min(maxbyvtx) - Sw(leftelem))*booleanrv + ...
            (1 - booleanrv)*(max(minbyvtx) - Sw(leftelem));
            
        %Get the ratio "rv":
        rv = deltasatval/sw_onmid_left;

        %Calculate "qsi" (It Choces according flag "cvbtype")
        qsi = (2*rv)/(rv^2 + 1);
%         qsi = min(1,rv);
            
    %There is a null projection
    else
        qsi = 1;
    end  %End of IF

    %Verify the final limiter (on the element)
    if qsi < mlplimiter(leftelem)
        mlplimiter(leftelem) = qsi;
    end  %End of IF

    %--------------------
    %ELEMENT on the RIGHT

    %Verify if there is a null projection and the very higher order is
    %turned off.
    %For the control volume on the left
    if abs(sw_onmid_right) > tol
        %Get the ratio "rv":
        %Define a boolean condition
        booleanrv = sw_onmid_right > 0;
        %Defien "deltasatval"
        deltasatval = (min(maxbyvtx) - Sw(rightelem))*booleanrv + ...
            (1 - booleanrv)*(max(minbyvtx) - Sw(rightelem));
            
        %Get the ratio "rv":
        rv = deltasatval/sw_onmid_right;

        %Calculate "qsi" (It Choces according flag "cvbtype")
        qsi = (2*rv)/(rv^2 + 1);
%         qsi = min(1,rv);
            
    %There is a null projection
    else
        qsi = 1;
    end  %End of IF

    %Verify the final limiter (on the element)
    if qsi < mlplimiter(rightelem)
        mlplimiter(rightelem) = qsi;
    end  %End of IF
end  %End of FOR ("inedge")

