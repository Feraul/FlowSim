function[phi] = calcEBMLP(numelem,elemeval,vertices,Sw,flagknownvert,grad,...
    satonvertices,centcoord,vtxcoord,ebtype)
%Define global parameters
global esurn1 esurn2

%Define a tolerance
tol = 1e-12;
%Initialize "qsy"
qsi = zeros(1,2);

    distcv = mean(vtxcoord(:,1:2)) - centcoord(1:2);

%Evaluate the vertices:
for ivtx = 1:2
    %Get the distance between the centroid and the vertex evaluated
%     distcv = vtxcoord(ivtx,1:2) - centcoord(1:2);
    esurn = ...
        esurn1(esurn2(vertices(ivtx)) + 1:esurn2(vertices(ivtx) + 1));

    %Verify if the node has a known saturation as boundary condi.
    if numelem < 2 && flagknownvert(vertices(ivtx)) == 1
        %Define saturation on the vertex evaluated
        satonverteval = satonvertices(vertices(ivtx));
        %Get the bigger saturatio value among those in the surrounding
        maxsatval = max([satonverteval; Sw(esurn)]);
        %Get the lower saturatio value among those in the surrounding
        minsatval = min([satonverteval; Sw(esurn)]);
        %Get the saturation extrapolated on the vertex
        sw_onvtx = satonvertices(vertices(ivtx)) - Sw(elemeval);

    %The vertex does not have a known boundary 
    else
        %Get the bigger saturatio value among those in the surrounding
        maxsatval = max(Sw(esurn));
        %Get the lower saturatio value among those in the surrounding
        minsatval = min(Sw(esurn));
        %Get the saturation extrapolated on the vertex
        sw_onvtx = dot(grad,distcv);
    end  %End of IF.
    
    %Verify if there is a null projection and the very higher order is
    %turned off.
    if abs(sw_onvtx) > tol
        %Get the ratio "rv":
        %Define a boolean condition
        booleanrv = sw_onvtx > 0;
        %Defien "deltasatval"
        deltasatval = (maxsatval - Sw(elemeval))*booleanrv + ...
            (1 - booleanrv)*(minsatval - Sw(elemeval));
            
        %Get the ratio "rv":
        rv = deltasatval/sw_onvtx;

        %Calculate "qsi" (It Choces according flag "cvbtype")
        switch ebtype
            %Original MLP (Park and Kin, 2012)
            case 'mlp'
                qsi(ivtx) = min(1,rv);
            %MLP u2, Venkatakrishnan (see Park and Kin, 2012)
            case 'mlp_vk'
                %Define "deltaplus" according "projval" sign.
                booleanrv = sw_onvtx > 0;
                %Defien "deltaplus"
                deltaplus = (maxsatval - Sw(elemeval))*booleanrv + ...
                    (1 - booleanrv)*(minsatval - Sw(elemeval));
                    
                %Define "deltaless"
                deltaless = sw_onvtx;
                %Delta sw (difference between the max and min 
                %saturation values):
                deltasonvertex = maxsatval - minsatval;
                %Calculate "teta"
%                 teta = deltasonvertex/(K2*(deltax^n));
                %Calculate epsilon^2:
                epsilonsqre = 1e-10;%(K1*(deltasonvertex^2))/(1 + teta);  %K*(deltax^3);%    
                    
                %Use Venkathakrishnen approximation
                qsi(ivtx) = (1/deltaless)*((((deltaplus^2) + epsilonsqre)*...
                    deltaless + 2*(deltaless^2)*deltaplus)/...
                    ((deltaplus^2) + 2*(deltaless^2) + ...
                    deltaless*deltaplus + epsilonsqre));
            %MLP-slc (Souza, Lyra, Carvalho)
            case 'mlp_slc'
                qsi(ivtx) = (rv^3 + 3*rv)/(rv^3 + (rv^2) + 4);
            %MLP-Van Albada
            case 'mlp_va'
                qsi(ivtx) = (rv^2 + rv)/(rv^2 + 1);
            %MLP-Van Albada 2 (Löhner, 2001)
            case 'mlp_va2'
                qsi(ivtx) = (2*rv)/(rv^2 + 1);
            %MLP-Van Leer 1 (TVD second order for any value of "rv")
            case 'mlp_vl'
                qsi(vtx) = 2*rv/(rv + 1);
            %MLP-Van Leer 2 (first-order for "rv" bigger than 1)
            case 'mlp_vl2'
                qsi(ivtx) = 4*rv/((rv + 1)^2);
            %MLP-New Limiter 2 (Mandal, 2008)
            case 'mlp_nl2'
                qsi(ivtx) = (rv^3 + rv)/(rv^3 + 1);
            %MLP-New Limiter 2 (Mandal, 2008)
            case 'mlp_nl3'
                qsi(ivtx) = (rv^2 + 3*rv)/(rv^2 + rv + 2);
        end  %End of SWITCH
            
    %There is a null projection
    else
        qsi(ivtx) = 1;
    end  %End of IF
end  %End of FOR

%Get the limiter for the edge
phi = min(qsi);
    