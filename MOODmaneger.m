%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 07/06/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%Define the MOOD strategy (see Clain et al., 2011).   

%--------------------------------------------------------------------------
%Additional comments: 
%It is called by the function "calcnewsatfield.m"

%--------------------------------------------------------------------------

function [pointelemeval,pntonlyelemeval,pointbndedg,pointinedg,...
    orderelemdist,orderbedgdist,orderinedgdist,entropyptr,dmp] = ...
    MOODmaneger(moodkey,Sw,newSw,injecelem,producelem,orderelemdist,...
    satinbound,entrineqterm,dt,entropyptr,flowrate)
%Define global parameters:
global elem bedge inedge bcflag satlimit elemarea pormap;

%Define a tolerance (it is like a zero for the computer)
tol = 1e-12;

%Initialize "elemsize"
elemsize = size(elem,1);
%Initialize "orderdtec_elemsur"
orderdtec_elemsur = zeros(elemsize,1);
orderdtect = orderdtec_elemsur;
%Initialize "elemvec", "pointelemeval", "bedgvec" and "inedgvec"
elemvec = 1:elemsize;
pointelemeval = elemvec;
pntonlyelemeval = elemvec;
bedgvec = 1:size(bedge,1);
inedgvec = 1:size(inedge,1);
%Attribute "1" to "numcol"
numcol = 1;

%There is PERIODIC BOUNDARY CONDITION.
if any(bcflag(:,1) > 600)
    periodicpos = getperiodicelem;
%There is no PERIODIC BOUNDARY CONDITION.
else
    periodicpos = 0;
end  %End of IF

%Take off the wells elements, when they exist (cell-centered well)
%It is a five-spot problem
if any(satinbound) == 0
    %Define the elements in "elem" which will be verified (injector well 
    %elements are out)
    pointelemeval_aux = setdiff(pointelemeval,injecelem,'stable');
%It is a Buckley-Leverett problem (injector surface)
else
    pointelemeval_aux = pointelemeval;    
end  %End of IF

%Swept all elements
for i = 1:length(pointelemeval_aux)
    %Define "ielem"
    ielem = pointelemeval_aux(i);
    %Get the elements which surround the element evaluated.
    [esureface,esurefull] = getsurelem(ielem);

    %Get periodic contribution for the "esurefull". Only Hyperbolic Eq.
    if any(periodicpos)
        %Attribute "2" to "numcol"
        numcol = 2;
        %Get a pointer for the "bedge" rows which have elements with
        %periodic boundary.
        pointperiodbedg = logical(bedge(:,3) == ielem);
        %Define additional elements:
        addelem = bedge(periodicpos(pointperiodbedg),3);
        %It changes "esurefull"
        esurefull = union(esurefull,addelem','stable');
    end  %End of IF
    
    %Define "Sw" on the element evaluated and on the vecinity
    %There is prescribed saturation on the boundary (Bickley-Leverett case)
    if any(satinbound) && ismember(ielem,injecelem)
        Swvicinity = [1 - satlimit(2); Sw([ielem esurefull])];
        Swvicinityproximo=[1 - satlimit(2); Sw([ielem esureface])];
        
    %There is NO prescribed saturation on the boundary.
    else
        Swvicinity = Sw([ielem esurefull]);
        Swvicinityproximo=Sw([ielem esureface]);
        
    end  %End of IF
    ielem
     [F]=ferncodes_element_face(ielem);
    m=1;h=1;
    for ii=F
        
        if ii>size(bedge,1)
            leftelem=inedge(ii-size(bedge,1),3);
            rightelem=inedge(ii-size(bedge,1),4);
            dotvn=flowrate(ii,1);
            if (ielem==leftelem)&& (dotvn>0)
                Sjus(m,1)=newSw(rightelem);
                m=m+1;
            elseif (ielem==leftelem)&& (dotvn<0)
                Smon(h,1)=newSw(rightelem);
                h=h+1;
            end
            
            if (ielem==rightelem)&& (dotvn>0)
                Smon(h,1)=newSw(leftelem);
                h=h+1;
            elseif (ielem==rightelem)&& (dotvn<0)
                Sjus(m,1)=newSw(leftelem);
                m=m+1;
            end
        elseif ii<=size(bedge,1) && bedge(ii,5)==202
            
            Smon(h,1)=1-satlimit(2);
            h=h+1;
        elseif ii<=size(bedge,1) && (bedge(ii,5)==101 || bedge(ii,5)==201 ) % revisar com mais cuidado
            Smon(h,1)=newSw(ielem);
            % Sjus(m,1)=newSw(ielem);
            %m=m+1;
            h=h+1;
        end
    end
    
        
        
        
    %It stabilish conditions for DMP violation. It take care to avoid trsh
    %values (values in the difference above 1e-12)
    %For the new value:
    Sweval = logical(abs(newSw(ielem)) > tol)*newSw(ielem);
    %For the old values:
    Swvicinity = logical(abs(Swvicinity) > tol).*Swvicinity;
    %======================================================================
    
    %It stabilish conditions for DMP violation. It take care to avoid trsh
    %values (values in the difference above 1e-12)
    %For the new value:
    Sweval = logical(abs(newSw(ielem)) > tol)*newSw(ielem);
    %For the old values:
    Swvicinityproximo = logical(abs(Swvicinityproximo) > tol).*Swvicinityproximo;
    %======================================================================

    if (max(Swvicinity)< Sweval || Sweval < min(Swvicinity)) %|| Sweval<Sw(ielem)
        %Update "orderdtec_elemsur" for the element evaluated
        %("0" and "1"). It marks the element own and its surrounding.
        orderdtec_elemsur([ielem esureface]) = 1;
        %It mark only the evaluated element.
        orderdtect(ielem) = 1;
        %Decremente the order in "orderelemdist" (... -> 4 -> 3 -> 2 -> 1).
        orderelemdist(ielem) = max(1,orderelemdist(ielem) - 1);
    %else
    elseif (max(Swvicinityproximo)< Sweval || Sweval < min(Swvicinityproximo))
        %Update "orderdtec_elemsur" for the element evaluated
        %("0" and "1"). It marks the element own and its surrounding.
        orderdtec_elemsur([ielem esureface]) = 1;
        %It mark only the evaluated element.
        orderdtect(ielem) = 1;
        %Decremente the order in "orderelemdist" (... -> 4 -> 3 -> 2 -> 1).
        orderelemdist(ielem) = max(1,orderelemdist(ielem) - 1);
    elseif max(Smon)<Sweval || Sweval<min(Sjus)
        %Update "orderdtec_elemsur" for the element evaluated
        %("0" and "1"). It marks the element own and its surrounding.
        orderdtec_elemsur([ielem esureface]) = 1;
        %It mark only the evaluated element.
        orderdtect(ielem) = 1;
        %Decremente the order in "orderelemdist" (... -> 4 -> 3 -> 2 -> 1).
        orderelemdist(ielem) = max(1,orderelemdist(ielem) - 1);
    end  %End of IF (DMP)
    %ielem
    clear Sjus Smon
end  %End of FOR

%     %Attribute the order to "ordermap" according "moodkey" value.
%     ordermap(ielem) = orderelemdist(ielem)*(moodkey == 1) + ...
%         min(orderelemdist([ielem esureface]))*(moodkey == 2); 

%Define again "pointelemeval" and "orderelemdist". It points to element 
%which does not respect the DMP condition (receives the number of the 
%elements).
pointelemeval = pointelemeval(logical(orderdtec_elemsur));

%"pntonlyelemeval" receives only the evaluated elements 
%(that ones problematics)
pntonlyelemeval = pntonlyelemeval(logical(orderdtect));

%Update "entropyptr":
%Update the non-marked positions
% entropyptr(logical(entropyptr == 1)) = 0;
%Update the marked positions
% entropyptr(logical(entropyptr > 1)) = 1;

%Define "pointbndedg". It points to "bedge" row whose element does 
%not respect the DMP condition (receives the number of the "bedge" row).
pointbndedg = bedgvec(logical(ismember(bedge(:,3),pointelemeval)));
%Initialize "orderbedgdist"
orderbedgdist = zeros(length(pointbndedg),numcol);

% pointelemeval
% length(pointelemeval)
%orderelemdist
%pause

%Swept the BOUNDARY edges
for i = 1:length(pointbndedg)
    %Define the real "bedge" row counter
    ibedg = pointbndedg(i);
    %Verify if there is a periodic boundary condition
    %There is periodic b.c.
    if any(periodicpos)
        %Define "elemleft" and "elemright" (periodic bound.)
        elemleft = bedge(ibedg,3);
        elemright = bedge(periodicpos(ibedg),3);
         if moodkey==0
             %Define the order for the elements on the left and on the right
            ordonleft = orderelemdist(elemleft);
            ordonright = orderelemdist(elemright);
            orderbedgdist(i,:) = [ordonleft ordonright];
        %Choose according "moodkey"
         elseif moodkey == 1
            %Define the order for the elements on the left and on the right
            ordonleft = orderelemdist(elemleft);
            ordonright = orderelemdist(elemright); 
            %Attribute it to "orderbedgdist"
            orderbedgdist(i,:) = min(ordonleft,ordonright);
        %MOOD type 2 (see Clain et al, 2011)
        elseif moodykey == 2
            %Get the elements which surround the element evaluated.
            [esurefaceleft,] = getsurelem(elemleft);
            [esurefaceright,] = getsurelem(elemright);
            %Define the order for the elements on the left and on the right
            ordonleft = ...
                min(orderelemdist([elemleft elemright esurefaceleft]));
            ordonright = ...
                min(orderelemdist([elemright elemleft esurefaceright])); 
            %Attribute it to "orderbedgdist"
            orderbedgdist(i,:) = [ordonleft ordonright];
        end  %End of IF
    %There is NO periodic b.c.
    else
        %Define "elemleft"
        elemleft = bedge(ibedg,3);
        %Get the elements which surround the element evaluated.
        [esurefaceleft,] = getsurelem(elemleft);
        %Attribute the order to the edge according to "moodkey" value
        orderbedgdist(i) = orderelemdist(elemleft)*(moodkey == 1) + orderelemdist(elemleft)*(moodkey == 0)+...
            min(orderelemdist([elemleft esurefaceleft]))*(moodkey == 2);
    end  %End of IF (evaluate periodic boundary condition)
end  %End of IF ("bedge")

%Define again "pointinedg". It points to "inedge" row whose element does 
%not respect the DMP condition (receives the number of the "inedge" row).
pointinedg = inedgvec(any(ismember(inedge(:,3:4),pointelemeval),2));

%Initialize "orderinedgdist"
orderinedgdist = zeros(length(pointinedg),2);

%Swept the INTERNAL edges
for i = 1:length(pointinedg)
    %Define the real "inedge" row counter
    inedg = pointinedg(i);
    %Define "elemleft" and "elemright"
    elemleft = inedge(inedg,3);
    elemright = inedge(inedg,4);
    if moodkey==0
        
            %Define the order for the elements on the left and on the right
            ordonleft = orderelemdist(elemleft);
            ordonright = orderelemdist(elemright);
            %Attribute it to "orderinedgdist"
            orderinedgdist(i,:) = [ordonleft ordonright];
        
    %Choose according "moodkey"
    elseif moodkey == 1
        %Define the order for the elements on the left and on the right
        ordonleft = orderelemdist(elemleft);
        ordonright = orderelemdist(elemright); 
        %Attribute it to "orderinedgdist"
        orderinedgdist(i,:) = min(ordonleft,ordonright);

    %MOOD type 2 (see Clain et al, 2011)
    elseif moodkey == 2
        %Get the elements which surround the element evaluated.
        [esurefaceleft,] = getsurelem(elemleft);
        [esurefaceright,] = getsurelem(elemright);
        %Define the order for the elements on the left and on the right
        ordonleft = min(orderelemdist([elemleft esurefaceleft]));
        ordonright = min(orderelemdist([elemright esurefaceright])); 
        %Attribute it to "orderinedgdist"
        orderinedgdist(i,:) = [ordonleft ordonright];
    end  %End of IF
end  %End of IF ("inedge")

%Verify if is necessary one more loop (if "dmp" == 0, it is not necessary).
dmp = any(orderdtec_elemsur);

