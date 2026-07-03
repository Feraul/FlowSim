%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 12/09/2013
%Modify data: 20/10/2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%It obtain both mobility and saturation over midege. It is necessary to 
%update the total mobility in pressure equation.    

%--------------------------------------------------------------------------
%Additional comments:
%The calculum of mobility depends of the option marked.

%1. Get the SATURATION through edge using VOLUME AVERAGE for each vertex 
    %and arithmetic average for saturation on midedge (So calculate the 
    %MOBILITY on the midedge).
%2. Get the SATURATION using VOLUME AVERAGE for each vertex, calculate the
    %MOBILITY on each vertex and arithmetic average of MOBILITY on the 
    %midedge).
%3. Get the SATURATION through midedge using FACE SHARED elements VOLUMES 
    %AVERAGE (So calculate the MOBILITY on the midedge with the saturation 
    %obtained).
%4. Get the SATURATION through edge using LINEAR INTERPOLATION for each 
    %vertex and arithmetic average for SATURATION on midedge. After that, 
    %calculate the mobilities with the saturation obtained.
%5. Get the SATURATION using LINEAR INTERPOLATION for each vertex, calcul.
    %the MOBILITY on each vertex and arithmetic average of MOBILITY on the 
    %midedge mobility on the midedge).
%6. Get the SATURATION using MUSCL-like projection, get a average 
    %SATURATION by INVERSE DISTANCE WEIGHTED and calculate the MOBILITY 
    %with the saturation calculated.

%7. Get the SATURATION using the saturation in a early iteraction. In this
    %case, for the first "timelevel" the saturation on the each edge is 
    %obtained using the strategy number "1". 

%8. Get the saturation for each Control Volume and calculate the MOBILITY 
    %through edge using MOBILITY VOLUME AVERAGE for each vertex and 
    %arithmetic average (by using the mobility into vertices) for MOBILITY 
    %on midedge.
    
%9. Get the saturation for each Control Volume and calculate the MOBILITY 
    %through edge using using FACE SHARED elements VOLUMES AVERAGE with the 
    %MOBILITY.
%10. Get the MOBILITY in each cell-centered, get the MOBILITY using LINEAR 
    %INTERPOLATION for each vertex and arithmetic average of MOBILITY on 
    %the midedge).

%--------------------------------------------------------------------------

function [mobility] = getmobility2(satinbound,injecelem,Sw,earlysw,...
    smethod,timelevel,benchkey,Sw_left,Sw_right,c,overedgecoord)
%Define global parameters:
global coord bedge inedge visc elemarea;

%Initialize "interpmobility". It is the type of strategy used for calculate 
%the MOBILITY.
interpmobility = 12;

%Initialize "bedgesize" and "inedge"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize again "mobility"
multipfac = 2*((interpmobility == 7 || interpmobility == 11) && ...
    (strcmp(smethod,'mwec') || strcmp(smethod,'mwic'))) + ...
    ((interpmobility ~= 7 && interpmobility ~= 11) || ...
    ((interpmobility == 7 || interpmobility == 11) && ...
    (strcmp(smethod,'mwec') == 0 && strcmp(smethod,'mwic') == 0)));

%Initialize "mobility"
mobility = zeros(multipfac*(bedgesize + inedgesize),2);

%Choose the Mobility Calculate Strategy according to "interpmobility"
switch interpmobility
    %Alternative 1:
    case 1
        %Get the saturation on the vertices and midedges
        [~,satinedges,~,flagknownedge] = getsatandflag(satinbound,...
            injecelem,Sw,1);

        %Get the relative permeability for water and oil (It uses the sat. 
        %on the MIDEDGES)
        [~,~,~,krw,kro,] = twophasevar(satinedges,benchkey);
        
        %Fill "mobility"
        for i = 1:length(satinedges)
            %Verify if there exists saturation prescribed on the boundary
            %There exists
            if i <= bedgesize && flagknownedge(i) == 1
                %Attribute WATER mobility
                mobility(i,1) = 1;
                %Attribute OIL mobility
                mobility(i,2) = 0;
            %There is NO saturation prescribed on the boundary
            else
                %Attribute WATER mobility
                mobility(i,1) = krw(i)./visc(1);
                %Attribute OIL mobility
                mobility(i,2) = kro(i)./visc(2);
            end  %End of IF
        end  %End of FOR (Swept all edges)
    
    %Alternative 2:
    case 2
        %Initialize "mobinvertices"
        mobinvertices = zeros(size(coord,1),2);
        
        %Get the saturation on the vertices (Volume Average)
        [satinvertices,] = getsatandflag(satinbound,injecelem,Sw,1);

        %Get the relative permeability for water and oil (It uses the satur. 
        %on the VERTICES)
        [~,~,~,krw,kro,] = twophasevar(satinvertices,benchkey);
        %Fill "mobinvertices"
        i = 1:length(satinvertices);
        %Attribute WATER mobility
        mobinvertices(i,1) = krw(i)./visc(1);
        %Attribute OIL mobility
        mobinvertices(i,2) = kro(i)./visc(2);
        
        %Calculate the MOBILITY on the midedges
        %Swept "bedge"
        for i = 1:size(bedge,1)
            %Water mobility
            mobility(i,1) = mean(mobinvertices(bedge(i,1:2),1));
            %Oil mobility
            mobility(i,2) = mean(mobinvertices(bedge(i,1:2),2));
        end  %End of FOR ("bedge")
            
        %Swept "inedge"
        for i = 1:size(inedge,1)
            bedgesize = size(bedge,1);
            %Water mobility
            mobility(bedgesize + i,1) = ...
                mean(mobinvertices(inedge(i,1:2),1));
            %Oil mobility
            mobility(bedgesize + i,2) = ...
                mean(mobinvertices(inedge(i,1:2),2));
        end  %End of FOR ("bedge")

    %Alternative 3:
    case 3
        %Get the saturation on the vertices and midedges ("interptype" = 2)
        [~,satinedges,] = getsatandflag(satinbound,injecelem,Sw,2);
    
        %Get the relative permeability for water and oil (It uses the sat. 
        %on the MIDEDGES)
        [~,~,~,krw,kro,] = twophasevar(satinedges,benchkey);
        
        %Fill "mobility"
        i = 1:length(satinedges);
        %Attribute WATER mobility
        mobility(i,1) = krw(i)./visc(1);
        %Attribute OIL mobility
        mobility(i,2) = kro(i)./visc(2);

    %Alternative 4:
    case 4
        %Get the relative permeability for water and oil (It uses the sat. 
        %on the control volume centroid)
        [~,~,~,krw,kro,] = twophasevar(Sw,benchkey);
        
        %Define lambda_w and lambda_o for the control volumes:
        lambda_w = krw/visc(1);
        lambda_o = kro/visc(2);
        
        %Get the mobility on the vertices and midedges ("interptype" = 2)
        [~,mob_w_onface,] = getsatandflag(satinbound,injecelem,lambda_w,2);
        [~,mob_o_onface,] = getsatandflag(satinbound,injecelem,lambda_o,2);
        
        %Attribute WATER mobility
        mobility(:,1) = mob_w_onface;
        %Attribute OIL mobility
        mobility(:,2) = mob_o_onface;
        
        %Clear auxiliary variables
        clear mob_w_onface mob_o_onface lambda_w lambda_o krw kro;
        
    %Alternative 5:
    case 5
        %Initialize "mobinvertices"
        mobinvertices = zeros(size(coord,1),2);
        
        %Get the saturation on the vertices (Linear Interpolation)
        [satinvertices,] = getsatandflag(satinbound,injecelem,Sw,3);

        %Get the relative permeability for water and oil (It uses the satur. 
        %on the VERTICES)
        [~,~,~,krw,kro,] = twophasevar(satinvertices,benchkey);
        %Fill "mobinvertices"
        i = 1:length(satinvertices);
        %Attribute WATER mobility
        mobinvertices(i,1) = krw(i)./visc(1);
        %Attribute OIL mobility
        mobinvertices(i,2) = kro(i)./visc(2);
        
        %Calculate the MOBILITY on the midedges
        %Swept "bedge"
        for i = 1:size(bedge,1)
            %Water mobility
            mobility(i,1) = mean(mobinvertices(bedge(i,1:2),1));
            %Oil mobility
            mobility(i,2) = mean(mobinvertices(bedge(i,1:2),2));
        end  %End of FOR ("bedge")
            
        %Swept "inedge"
        for i = 1:size(inedge,1)
            bedgesize = size(bedge,1);
            %Water mobility
            mobility(bedgesize + i,1) = ...
                mean(mobinvertices(inedge(i,1:2),1));
            %Oil mobility
            mobility(bedgesize + i,2) = ...
                mean(mobinvertices(inedge(i,1:2),2));
        end  %End of FOR ("bedge")

    %Alternative 7 (Get MOBILITY from saturation in early time iteraction):
    case 7
        %Use the first strategy (1) for the initial timelevel:
        %Get the saturation on the vertices and midedges
        [~,satinedges_aux,~,flagknownedge] = ...
            getsatandflag(satinbound,injecelem,Sw,1);
        %Verify the numerical scheme:
        %When we have a multidimensional scheme with half-edges loop.
        if strcmp(smethod,'mwec') || strcmp(smethod,'mwic')
            %Verify if the iteraction ("timelevel") is the first
            if timelevel == 1
                %Initialize "satinedges" and "m"
                satinedges = zeros(size(mobility,1),1);
                m = 0;
                %Distribute "satinedges" which the size is equal to "bedge" +
                %"inedge" to "earlysw" size which is the double.
                for i = 1:length(satinedges_aux)
                    satinedges(m + 1:m + 2) = satinedges_aux(i);
                    %Update "m"
                    m = m + 2;
                end  %End of FOR
            %Any other timelevel.
            else
                satinedges = earlysw;
            end  %End of IF

            %Get the relative permeability for water and oil (It uses the sat. 
            %on the MIDEDGES)
            [~,~,~,krw,kro,] = twophasevar(satinedges,benchkey);

            %Initialize "j" and "c"
            j = 1;
            c = 1;
            %Fill "mobility"
            for i = 1:size(mobility,1)
                %Verify if there exists saturation prescribed on the boundary
                %There exists
                if i <= 2*bedgesize && flagknownedge(j) == 1
                    %Attribute WATER mobility
                    mobility(i,1) = 1;
                    %Attribute OIL mobility
                    mobility(i,2) = 0;
                %There is NO saturation prescribed on the boundary
                else
                    %Attribute WATER mobility
                    mobility(i,1) = krw(i)/visc(1);
                    %Attribute OIL mobility
                    mobility(i,2) = kro(i)/visc(2);
                end  %End of IF
            
                %Update "j"
                if c == 2
                    j = j + 1;
                    c = 0;
                end  %end of IF
                %Update "c"
                c = c + 1;
            end  %End of FOR (Swept all edges)
        
        %Any other scheme. It uses the whole edge, not half edges.
        else
            %Verify if the iteraction ("timelevel") is the first
            if timelevel == 1
                %Initialize "satinedges" and "m"
                satinedges = satinedges_aux;
            %Any other timelevel.
            else
                satinedges = earlysw;
            end  %End of IF

            %Get the relative permeability for water and oil (It uses the sat. 
            %on the MIDEDGES)
            [~,~,~,krw,kro,] = twophasevar(satinedges,benchkey);
        
            %Fill "mobility"
            for i = 1:size(mobility,1)
                %Verify if there exists saturation prescribed on the boundary
                %There exists
                if i <= bedgesize && flagknownedge(i) == 1
                    %Attribute WATER mobility
                    mobility(i,1) = 1;
                    %Attribute OIL mobility
                    mobility(i,2) = 0;
                %There is NO saturation prescribed on the boundary
                else
                    %Attribute WATER mobility
                    mobility(i,1) = krw(i)/visc(1);
                    %Attribute OIL mobility
                    mobility(i,2) = kro(i)/visc(2);
                end  %End of IF
            end  %End of FOR (Swept all edges)
        end  %End of IF
        
    %Alternative 8: Get the MOBILITY for each control volume surrounding 
    %the vertices "i" and "j". After that, calculate the mean value on the 
    %mid-point on the evaluated face.
    case 8
        %Initialize "mobilityoncentroid". In first column stores the water
        %mobility, in second column stores the oil mobility.
        mobilityoncentroid = zeros(length(Sw),2);
        %Initialize "mobinbound"
        mobinbound = zeros(length(satinbound),2);
        
        %Get the relative permeability for water and oil (It uses the sat. 
        %on the Cell-center)
        [~,~,~,krw,kro,] = twophasevar(Sw,benchkey);
        
        %Fill "mobilityoncentroid"
        i = 1:length(Sw);
        %Attribute WATER mobility
        mobilityoncentroid(i,1) = krw(i)./visc(1);
        %Attribute OIL mobility
        mobilityoncentroid(i,2) = kro(i)./visc(2);

        %Verify if there is saturation attributed to boundary
        if any(satinbound) && length(satinbound) > 1
            %Get the relative permeability for water and oil on the 
            %boundary (It uses the sat. on the Cell-center)
%             [null1,null2,null3,krw,kro,] = twophasevar(satinbound,benchkey);
%             %Clear the null variable
%             clear null1 null2 null3;
            
            %Fill "mobinbound"
            i = 1:length(satinbound);
            %Attribute WATER mobility
            mobinbound(i,1) = 1;%krw(i)./visc(1);
            %Attribute OIL mobility
            mobinbound(i,2) = 0;%kro(i)./visc(2);
        end  %End of IF

        %Get the Water MOBILITY on the midedges
        [~,watermob,] = getsatandflag(mobinbound(:,1),injecelem,...
            mobilityoncentroid(:,1),1);
    
        %Get the Oil MOBILITY on the midedges
        [~,oilmob,] = getsatandflag(mobinbound(:,2),injecelem,...
            mobilityoncentroid(:,2),1);

        %Fill "mobility"
        mobility(:,1:2) = [watermob oilmob];
        
    %Alternative 9: Get the SATURATION for each control volume surrounding 
    %the vertices "i" and "j". After that, calculate the mean value on the 
    %mid-point on the evaluated face.
    case 9
        %Get the Water saturation on the mid-edges
        [~,swonedge,] = getsatandflag(satinbound,injecelem,Sw,1);
        
        %Get the relative permeability for water and oil (It uses the sat. 
        %on the Cell-center)
        [~,~,~,krw,kro,] = twophasevar(swonedge,benchkey);
        
        %Fill "mobility"
        mobility(:,1:2) = [(krw/visc(1)) (kro/visc(2))];

    %Alternative 10 (Get MOBILITY instead saturation) - Linear Interpol.:
    case 10
        %Initialize "mobilityoncentroid". In first column stores the water
        %mobility, in second column stores the oil mobility.
        mobilityoncentroid = zeros(length(Sw),2);
        %Initialize "mobinbound"
        mobinbound = zeros(length(satinbound),2);
        
        %Get the relative permeability for water and oil (It uses the sat. 
        %on the Cell-center)
        [~,~,~,krw,kro,] = twophasevar(Sw,benchkey);
        
        %Fill "mobilityoncentroid"
        i = 1:length(Sw);
        %Attribute WATER mobility
        mobilityoncentroid(i,1) = krw(i)./visc(1);
        %Attribute OIL mobility
        mobilityoncentroid(i,2) = kro(i)./visc(2);

        %Verify if there is saturation attributed to boundary
        if any(satinbound) && length(satinbound) > 1
            %Get the relative permeability for water and oil on the 
            %boundary (It uses the sat. on the Cell-center)
            [~,~,~,krw,kro,] = twophasevar(satinbound,...
                benchkey);
            
            %Fill "mobinbound"
            i = 1:length(satinbound);
            %Attribute WATER mobility
            mobinbound(i,1) = krw(i)./visc(1);
            %Attribute OIL mobility
            mobinbound(i,2) = kro(i)./visc(2);
        end  %End of IF

        %Get the Water MOBILITY on the vertices ("interptype" = 2)
        [watermobinvert,] = getsatandflag(mobinbound(:,1),injecelem,...
            mobilityoncentroid(:,1),3);
    
        %Get the Water MOBILITY on the vertices ("interptype" = 2)
        [oilmobinvert,] = getsatandflag(mobinbound(:,2),injecelem,...
            mobilityoncentroid(:,2),3);

        %Calculate the MOBILITY on the midedges
        %Swept "bedge"
        for i = 1:size(bedge,1)
            %Water mobility
            mobility(i,1) = mean(watermobinvert(bedge(i,1:2)));
            %Oil mobility
            mobility(i,2) = mean(oilmobinvert(bedge(i,1:2)));
        end  %End of FOR ("bedge")
            
        %Swept "inedge"
        for i = 1:size(inedge,1)
            bedgesize = size(bedge,1);
            %Water mobility
            mobility(bedgesize + i,1) = ...
                mean(watermobinvert(inedge(i,1:2)));
            %Oil mobility
            mobility(bedgesize + i,2) = ...
                mean(oilmobinvert(inedge(i,1:2)));
        end  %End of FOR ("inedge")

    %Alternative 11 (Get MOBILITY instead saturation):
    case 11
        %Initialize "mobilityoncentroid". In first column stores the water
        %mobility, in second column stores the oil mobility.
        mobilityoncentroid = zeros(length(Sw),2);
        %Initialize "mobinbound"
        mobinbound = zeros(length(satinbound),2);
        
        %Get the relative permeability for water and oil (It uses the sat. 
        %on the Cell-center)
        [~,~,~,krw,kro,] = twophasevar(Sw,benchkey);
        
        %Fill "mobilityoncentroid"
        i = 1:length(Sw);
        %Attribute WATER mobility
        mobilityoncentroid(i,1) = krw(i)./visc(1);
        %Attribute OIL mobility
        mobilityoncentroid(i,2) = kro(i)./visc(2);

        %Verify if there is saturation attributed to boundary
        if any(satinbound) && length(satinbound) > 1
            %Get the relative permeability for water and oil on the 
            %boundary (It uses the sat. on the Cell-center)
%             [null1,null2,null3,krw,kro,] = twophasevar(satinbound,benchkey);
%             %Clear the null variable
%             clear null1 null2 null3;
            
            %Fill "mobinbound"
            i = 1:length(satinbound);
            %Attribute WATER mobility
            mobinbound(i,1) = 1;%krw(i)./visc(1);
            %Attribute OIL mobility
            mobinbound(i,2) = 0;%kro(i)./visc(2);
        end  %End of IF

        %Get the Water MOBILITY on the midedges
        [watermob_vtx,watermob_face,] = getsatandflag(mobinbound(:,1),injecelem,...
            mobilityoncentroid(:,1),1);
    
        %Get the Oil MOBILITY on the midedges
        [oilmob_vtx,oilmob_face,] = getsatandflag(mobinbound(:,2),injecelem,...
            mobilityoncentroid(:,2),1);

        c = 0;
        %Fill "mobility"
        %Swept "bedge"
        for j = 1:bedgesize
            %Get vertices
            vertices = bedge(j,1:2);
            %Water
            watermob1 = (watermob_vtx(vertices(1)) + watermob_face(j))/2; 
            watermob2 = (watermob_vtx(vertices(2)) + watermob_face(j))/2; 
            %Oil
            oilmob1 = (oilmob_vtx(vertices(1)) + oilmob_face(j))/2; 
            oilmob2 = (oilmob_vtx(vertices(2)) + oilmob_face(j))/2; 
            
            mobility(c + 1:c + 2,:) = ...
                [watermob1 oilmob1; watermob2 oilmob2];
            
            c = c + 2;
        end  %End of FOR

        %Swept "inedge"
        for j = 1:inedgesize
            %Get vertices
            vertices = inedge(j,1:2);
            %Water
            watermob1 = (watermob_vtx(vertices(1)) + ...
                watermob_face(j + bedgesize))/2; 
            watermob2 = (watermob_vtx(vertices(2)) + ...
                watermob_face(j + bedgesize))/2; 
            %Oil
            oilmob1 = (oilmob_vtx(vertices(1)) + ...
                oilmob_face(j + bedgesize))/2; 
            oilmob2 = (oilmob_vtx(vertices(2)) + ...
                oilmob_face(j + bedgesize))/2; 
            
            mobility(c + 1:c + 2,:) = ...
                [watermob1 oilmob1; watermob2 oilmob2];
            
            c = c + 2;
        end  %End of FOR

    %Alternative 12: Volume average with the elements that shares the face.
    %It uses extrapolated value up to the quadrature point instead the 
    %finite volume value over the element. 
    case 12
        %Verify the numerical scheme:
        %When we have a multidimensional scheme with half-edges loop.
%         if strcmp(smethod,'mwec') || strcmp(smethod,'mwic')
            %Initialize the vectors "satonvertices" and "satonedges"
            satonedges = zeros(bedgesize + inedgesize,1);
            %Verify if the iteraction ("timelevel") is the first
            if timelevel == 1 || (c == 0 && timelevel > 1)
                %Get the saturation on the vertices and midedges ("interptype" = 2)
                [~,satonedges_aux,] = getsatandflag(satinbound,injecelem,Sw,2);
                satonedges = satonedges_aux;
                clear satonedges_aux;
            %Any other timelevel.
            else
                %A standard formulation is in use. The quadrature point is
                %on the mid-point of the edge.
                if length(Sw_left) == (bedgesize + inedgesize) 
                    %Swept "bedge"
                    %Calculate the saturation into each unknown edge 
                    %("bedge")    
                    satonedges(1:bedgesize) = Sw_left(1:bedgesize);    

                    %Swept "inedge"
                    %Calculate the saturation into each unknown edge 
                    %("inedge")    
                    satonedges(bedgesize + 1:bedgesize + inedgesize) = ...
                        mean([Sw_left(bedgesize + 1:...
                        bedgesize + inedgesize) Sw_right],2); 
                
                %A MultiD formulation is in use. The quadrature point is
                %on the mid-point of the HALF-edge.
                else
                    %Swept "bedge"
                    %Calculate the saturation into each unknown edge 
                    %("bedge")    
                    satonedges(1:bedgesize) = Sw_left(1:2:2*bedgesize);    
                    satonedges(1:bedgesize) = (satonedges(1:bedgesize) + ...
                        Sw_left(2:2:2*bedgesize))/2;    

                    %Swept "inedge"
                    %Calculate the saturation into each unknown edge 
                    %("inedge")    
                    meansat = mean([Sw_left(2*bedgesize + 1:...
                        2*(bedgesize + inedgesize)) Sw_right],2); 
                    %Attribute to "satonedges" the mean of states for the
                    %first half-surface on the surface.
                    satonedges(bedgesize + 1:bedgesize + inedgesize) = ...
                        meansat(1:2:length(meansat));
                    %Attribute to "satonedges" the mean of both states for 
                    %the whole surface.
                    satonedges(bedgesize + 1:bedgesize + inedgesize) = ...
                        (satonedges(bedgesize + 1:bedgesize + inedgesize) + ...
                        meansat(2:2:length(meansat)))/2;
                end  %End of IF
            end  %End of IF
            
            %Get the relative permeability for water and oil 
            %(It uses the sat on the MIDEDGES)
            [~,~,~,krw,kro,] = twophasevar2(satonedges,benchkey,...
                overedgecoord(:,1),injecelem);

            %Fill "mobility"
            %Attribute WATER mobility
            mobility(:,1) = krw/visc(1);
            %Attribute OIL mobility
            mobility(:,2) = kro/visc(2);
%         else
%         end  %End of IF        


    %Alternative 12: Volume average with the elements that shares the face.
    %It uses extrapolated value up to the quadrature point instead the 
    %finite volume value over the element. 
    case 13
        %Verify the numerical scheme:
        %When we have a multidimensional scheme with half-edges loop.
%         if strcmp(smethod,'mwec') || strcmp(smethod,'mwic')
            %Initialize the vectors "satonvertices" and "satonedges"
            satonedges = zeros(bedgesize + inedgesize,1);
            %Verify if the iteraction ("timelevel") is the first
            if timelevel == 1
                %Get the saturation on the vertices and midedges ("interptype" = 2)
                [~,satonedges_aux,] = getsatandflag(satinbound,injecelem,Sw,2);
                satonedges = satonedges_aux;
                clear satonedges_aux;
                
                %Get the relative permeability for water and oil (It uses the sat. 
                %on the MIDEDGES)
                [~,~,~,krw,kro,] = twophasevar(satonedges,benchkey);
        
                %Fill "mobility"
                %Attribute WATER mobility
                mobility(:,1) = krw/visc(1);
                %Attribute OIL mobility
                mobility(:,2) = kro/visc(2);
                %Any other timelevel.
            else
                %A standard formulation is in use. The quadrature point is
                %on the mid-point of the edge.
                if length(Sw_left) == (bedgesize + inedgesize) 
                    %Swept "bedge"
                    %Calculate the saturation into each unknown edge 
                    %("bedge")    
                    satonedges(1:bedgesize) = Sw_left(1:bedgesize);    

                    %Swept "inedge"
                    %Calculate the saturation into each unknown edge 
                    %("inedge")    
                    satonedges(bedgesize + 1:bedgesize + inedgesize) = ...
                        mean([Sw_left(bedgesize + 1:...
                        bedgesize + inedgesize) Sw_right],2); 
                
                %A MultiD formulation is in use. The quadrature point is
                %on the mid-point of the HALF-edge.
                else
                    %Swept "bedge"
                    %Calculate the saturation into each unknown edge 
                    %("bedge")    
                    %Get mobility on the left:
                    [~,~,~,krw_left,kro_left,] = ...
                        twophasevar(Sw_left(1:2*bedgesize),benchkey);
                    %Water
                    mobility(1:bedgesize,1) = ...
                        krw_left(1:2:2*bedgesize)/visc(1); 
                    mobility(1:bedgesize,1) = (mobility(1:bedgesize,1) + ...
                        krw_left(2:2:2*bedgesize)/visc(1))/2; 
                    %Oil
                    mobility(1:bedgesize,2) = ...
                        kro_left(1:2:2*bedgesize)/visc(2); 
                    mobility(1:bedgesize,2) = (mobility(1:bedgesize,2) + ...
                        kro_left(2:2:2*bedgesize)/visc(2))/2; 

                    %Swept "inedge"
                    %Calculate the saturation into each unknown edge 
                    %("inedge")    
                    %Get mobility on the left:
                    [~,~,~,krw_left,kro_left,] = ...
                        twophasevar(Sw_left(2*bedgesize + 1:...
                        2*(bedgesize + inedgesize)),benchkey);
                    %Get mobility on the left:
                    [~,~,~,krw_right,kro_right,] = ...
                        twophasevar(Sw_right,benchkey);

                    meanlamb_w = mean([(krw_left/visc(1)) (krw_right/visc(1))],2); 
                    meanlamb_o = mean([(kro_left/visc(2)) (kro_right/visc(2))],2); 
                    %Attribute to "satonedges" the mean of states for the
                    %first half-surface on the surface.
                    %Water
                    mobility(bedgesize + 1:bedgesize + inedgesize,1) = ...
                        meanlamb_w(1:2:length(meanlamb_w));
                    mobility(bedgesize + 1:bedgesize + inedgesize,1) = ...
                        (mobility(bedgesize + 1:bedgesize + inedgesize,1) + ...
                        meanlamb_w(2:2:length(meanlamb_w)))/2;
                    %Oil
                    mobility(bedgesize + 1:bedgesize + inedgesize,2) = ...
                        meanlamb_o(1:2:length(meanlamb_w));
                    mobility(bedgesize + 1:bedgesize + inedgesize,2) = ...
                        (mobility(bedgesize + 1:bedgesize + inedgesize,2) + ...
                        meanlamb_o(2:2:length(meanlamb_w)))/2;
                end  %End of IF
            end  %End of IF
end  %End of SWITCH




