%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 25/03/2016 (Passion Fryday)
%Modify data:   /  /2016
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: 

%--------------------------------------------------------------------------
%Additional comments: 
%This function is called by the function "calcnewsatfield.m"

%--------------------------------------------------------------------------

function [multidsworfw,nonmultidsworfw] = ...
    getMultiDsat_padmec(sw_md,sw_nmd,mweightvec,halfedgepos,...
    coordmaprodelem,wellprodkey,multidsign,amthe_well)
%Define global parameters:
global coord inedge bedge nsurn1 nsurn2;

%Get the size of "coord", "bedge" and "inedge"
%Get the size of "bedge" and "coord"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%The "coordsize" changes according to "wellprodvalue" (1 --> evaluate only 
%the wells, 0 --> evaluate all elements)
coordsize = (1 - wellprodkey)*size(coord,1) + ...
    wellprodkey*length(coordmaprodelem);

%Define the size of "multidsworfw" and "nonmultidsworfw" according to well
%treatment ("wellprodkey").
storemultidsize = 2*(bedgesize + inedgesize)*(1 - wellprodkey) + ...
    wellprodkey;
%Initialize "multidsworfw". It stores the MULTID saturation in each 
%half-edge on boundary and inside the demoin.
multidsworfw = zeros(storemultidsize,1);
%Initialize "nonnmultidsworfw". It stores the saturation in each half-edge
%calculated from side non MultiD.
nonmultidsworfw = multidsworfw;

%Verify if ther is a well treatment
if wellprodkey == 1
    %Initialize "hesequence" and "hesorted"
    hesequence = 1:length(halfedgepos);
    hesorted = sort(halfedgepos);    
end  %End of IF

%Initialize the auxiliary counter "m"
m = 1;

%Swept all vertices (for each one there is an interaction reg.)
for j = 1:coordsize
    %Define "inode" according "wellprodkey"
    %It is not restrict to producer well
    if wellprodkey == 0
        inode = j;
    %It is restrict to producer well
    else
        inode = coordmaprodelem(j);
    end  %End of IF

    %Verify if the treatment of producer well is turned ON
    %It is Turned OFF (conventional)
    if wellprodkey == 0
        %It catches the amount of elements and nodes surrounding each vertex.
        nsurn = nsurn1(nsurn2(inode) + 1:nsurn2(inode + 1));
        localsw_md = sw_md(nsurn2(inode) + 1:nsurn2(inode + 1));
        localmweight = mweightvec(nsurn2(inode) + 1:nsurn2(inode + 1));
        localmdflag = multidsign(nsurn2(inode) + 1:nsurn2(inode + 1));
        
        %It catches the length of "nsurn"
        lengnsurn = length(nsurn);
    %It is Turned ON
    else
        %It catches the length of a type of "nsurn" with the half-edges
        %associated to the producer well
        lengnsurn = amthe_well(j);
    end  %End of IF

    %Initialize "localmdsatinhalfedge", "anothersatinhalfedge" and 
    %"halfedgepos"
    localmdsatinhalfedge = zeros(lengnsurn,1);
    anothersatinhalfedge = localmdsatinhalfedge;
    localhalfedgepos = localmdsatinhalfedge;

    %It swepts the halfedges surrounding the vertex evaluated.
    for i = 1:lengnsurn
        %Get the MultiD variable ("Sw" or "fw")
        %Evaluate if the way is on the counter clockwise or clockwise

        %Get the length of the each local vector
        lengvec = length(localsw_md);
        
        %Define "booleancw"
        booleancw = (localmdflag(i) == 1);

        %Get an auxiliary vector for the local variables 
        %Saturation sequence:
        aux_sw_md = ...
            booleancw*[localsw_md(i:lengvec); localsw_md(1:i - 1)] + ...
            (1 - booleancw)*[localsw_md(i:-1:1); ...
            localsw_md(lengvec:-1:i + 1)];
        %Weight sequence:
        aux_mweight = ...
            booleancw*[localmweight(i:lengvec); localmweight(1:i - 1)] + ...
            (1 - booleancw)*[localmweight(i:-1:1); ...
            localmweight(lengvec:-1:i + 1)]; 
            
        %---------------------------------
        %Calculate the Convex Combination:
        
        %MultiD approximation:
        %Define multip. factor
        multfact = (1/(1 - prod(aux_mweight)));
        
        %Calculate the last term
        mdcount = (1 - aux_mweight(lengvec))*aux_sw_md(lengvec);
        %Calculate another terms
        for imd = lengvec - 1:-1:1
            mdcount = mdcount*aux_mweight(imd) + ...
                (1 - aux_mweight(imd))*aux_sw_md(imd);
        end  %End of FOR
        
        %Get the final expression:
        localmdsatinhalfedge(i) = multfact*mdcount;
        %Non MultiD variable ("Sw" or "fw")
        anothersatinhalfedge(i) = sw_nmd(m);  %ok 
        
        %Get the half-edge position. It helps to put the variable in its
        %correct place (in "multidsworfw").
        %There is NOT Well Treatment
        if wellprodkey == 0
            %Fill "localhalfedgepos" with global half-edges
            localhalfedgepos(i) = halfedgepos(m);
        %There is a Well Treatment
        else
            %Fill "localhalfedgepos" with local half-edges
            localhalfedgepos(i) = ...
                hesequence(logical(halfedgepos(m) == hesorted));
        end  %End of IF

        %Update the auxiliary counters "c" and "m"
        m = m + 1;
    end  %End of FOR (swept the half-edges)
    
    %Alocate "localmdsatinhalfedge" in the vector "multidsworfw".
    multidsworfw(localhalfedgepos) = localmdsatinhalfedge;
    %Alocate "anothersatinhalfedge" in the vector "nonmultidsworfw".
    nonmultidsworfw(localhalfedgepos) = anothersatinhalfedge;
end  %End of FOR (swept the vertices)




% b = [a(pos:length(a)) a(1:pos - 1)]
% c = [a(pos:-1:1) a(length(a):-1:pos + 1)]