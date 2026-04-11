%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 06/03/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: This function calculate the saturation on the adjacent half-edge
%and the surrounding half-edge (more distant half-edge). In adition, this
%function get the weight (for the multidimensionality) and the position of
%each half-edge considered.

%--------------------------------------------------------------------------
%Additional comments: 
%This function is called by the function "calcnewsatfield.m"

%--------------------------------------------------------------------------

function [satonadjsuredge,satonanotherstate,mweightvec,halfedgepos,...
    amthe_well] = getmultidsatweig(Sw,taylorterms,limiterflag,flowrate,...
    mobility,flagknownedge,satonedges,massweigmap,othervertexmap,...
    multdlimiter,constraint,flagknownvert,satonvertices,mlplimiter,...
    coordmaprodelem,wellprodkey)
%Define global parameters:
global coord centelem inedge bedge order rowposit nsurn1 nsurn2;

%Define "tol"
tol = 1e-12;
%Get the size of "bedge" and "coord"
bedgesize = size(bedge,1);
%The "coordsize" changes according to "wellprodvalue" (1 --> evaluate only 
%the wells, 0 --> evaluate all elements)
coordsize = (1 - wellprodkey)*size(coord,1) + ...
    wellprodkey*length(coordmaprodelem);
%Initialize "satonadjsuredge" and "mweightvec". They store, respectively,
%the saturation in adjedge and suredge (in sequence) and the weights.
satonadjsuredge = zeros(2*length(nsurn1),1);
%"satonanotherstate" stores the state on the other side of the half-edge
%evaluated. We verify the side where the state is MultiD (store it in 
%"satonadjsuredge") and the another is stored into "satonanotherstate".
satonanotherstate = zeros(length(nsurn1),1);   %zeros(length(nsurn1),1);
mweightvec = satonanotherstate; 
halfedgepos = satonanotherstate;

%Initialize the auxiliary counters "m", "l"
m = 0;
l = 1;
amthe_well = 0;

%Initialize a counter
vtxcount = 0;
%Swept all vertices (for each one there is an interaction reg.)
for j = 1:coordsize
    %Initialize "u". It is a counter of the amount of half-edges taken in
    %account from "massweigmap". It fills "amthe_well".
    u = 0;
    %Define "inode" according "wellprodkey"
    %There is NOT WELL treatment
    if wellprodkey == 0
        inode = j;
    %There EXISTS WELL treatment
    else
        inode = coordmaprodelem(j);
    end  %End of IF

    %It catches the amount of elements and nodes surrounding each vertex.
    nsurn = nsurn1(nsurn2(inode) + 1:nsurn2(inode + 1));
    localrowposit = rowposit(nsurn2(inode) + 1:nsurn2(inode + 1));
    
    %It catches the length of "nsurn"
    lengnsurn = length(nsurn);
    
    %Get the "inode" coordinate
    inodecoord = coord(inode,1:2);

    %It swepts the halfedges surrounding the vertex evaluated.
    for i = 1:lengnsurn
        %Initialize "getvalkey". It allows the increment of "c", "m" and 
        %"l"
        getvalkey = 0;
        %Get "nodesur"
        nodesur = nsurn(i);
    
        %Get the "inode" coordinate
        nodesurcoord = coord(nodesur,1:2);
        %Get a vector correnponding to edge evaluated
        evalvec = nodesurcoord - inodecoord;
        
        selectednode = [inode nodesur];
        
        %The half-edge belongs to "bedge"
        if localrowposit(i) <= bedgesize
            %Define the vertices
            vertices = bedge(localrowposit(i),1:2);
            
            %Verify if the command below must be done (for "bedge" half-ed)
            if (wellprodkey == 1 && all(ismember(vertices,selectednode))) ...
                    || wellprodkey == 0
                %Get the coordinate of the vertices
                verticescoord = coord(vertices,:);
                %Define the control volume on the left.
                leftelem = bedge(localrowposit(i),3);
                %Define the adjacent element ("adjelem")
                adjelem = leftelem;

                %--------------------------------
                %Flow rate in UPSTREAM half-edge

                %Calculate the weight:                
                mweight = 0;
                surelem = adjelem;

                %Verify if there exists a boundary cond. (prescribed 
                %saturation).
                %Obs.: It is used in Buckley-Leverett applications.
                booleankey = (flagknownedge(localrowposit(i)) == 1);            
                %Attribute to "knownsat" the boundary saturation value.
                knownsat = [satonedges(localrowposit(i)) 0];

                %auxiliary "verticescoord"
%                 auxvertcoord = [coord(inode,:); mean(verticescoord,1)];
            
                %Get the saturation in the:  
                mlpbyelem = ...
                    mlplimiter((order == 1) + (order > 1)*adjelem,:);
                %1. adjacent half-edge:
                swonadjedge = getsatonedge(adjelem,vertices,verticescoord,...
                    taylorterms,Sw,limiterflag,order,constraint,...
                    flagknownvert,satonvertices,mlpbyelem,...
                    centelem(adjelem,1:2));
%                 swonadjedge = getsatonedge(adjelem,vertices,auxvertcoord,...
%                     taylorterms,Sw,limiterflag,order,constraint,flagknownvert,...
%                     satonvertices,mlpbyelem,centelem(adjelem,1:2));
                %2. surounding half-edge:
                %Define the vertices of the distant half-edge
%                 vertdistant = [inode othervertex(1)];
%                 %Get the coordinate of the distant half-edge
%                 vertdistantcoord = coord(vertdistant,:);

                %auxiliary "vertdistantcoord"
%                 auxvertdstcoord = [coord(inode,:); mean(vertdistantcoord,1)];
            
                %Calculate the saturation on the half-edge.
%                 swonsuredge = getsatonedge(surelem,vertdistant,...
%                     vertdistantcoord,taylorterms,Sw,limiterflag,order,...
%                     constraint,flagknownvert,satonvertices,...
%                     mlplimiter(surelem(1),:),centelem(surelem,1:2)); 
%                 swonsuredge = getsatonedge(surelem,vertdistant,...
%                     auxvertdstcoord,taylorterms,Sw,limiterflag,order,...
%                     constraint,flagknownvert,satonvertices,...
%                     mlplimiter(surelem(1)),centelem(surelem,1:2)); 
                %Attribute to "unknowsat" the saturation value calculated  
                %by the function "getsatonedge".
%                 unknownsat = [swonadjedge swonsuredge];  
                unknownsat = [swonadjedge 0];  
            
                %Store the position saturation in "adjedge" and "suredge" 
                %(in sequence)
                satonadjsuredge(m + 1:m + 2) = knownsat*booleankey + ...
                    unknownsat*(1 - booleankey);
                %Store the weight calculated.
                mweightvec(l) = mweight*(1 - booleankey);
                %Store the position of the half-edge. It is used for store 
                %the saturation in each half-edge in a unic vector.
                halfedgepos(l) = localrowposit(i);
            
                %Turn "getvalkey" on
                getvalkey = 1;
            end  %End of IF (maybe well treatment)
        
        %The half-edge belongs to "inedge"
        else
            %Define
            inlocrowpos = localrowposit(i) - bedgesize;
            %Define the vertices
            vertices = inedge(inlocrowpos,1:2);
    
            %Verify if the command below must be done (for "inedge" half-e)
            if (wellprodkey == 1 && all(ismember(vertices,selectednode))) ...
                    || wellprodkey == 0
                %Get the coordinate of the vertices
                verticescoord = coord(vertices,:);
                %Define the control volume on the left and on the right.
                leftelem = inedge(inlocrowpos,3);
                rightelem = inedge(inlocrowpos,4);

                booleankey = (inode == inedge(inlocrowpos,1));
                %Get the "flowrate" position
                frposition = 2*localrowposit(i) - booleankey;
                %Get the local flowrate
                localflowrate = flowrate(frposition);

                %It rounds the flow rate value
                booleanfr = abs(localflowrate) >= tol;
                %Define again "localflowrate"
                localflowrate = booleanfr*localflowrate;

                %Define boolean conditions
                booleanfr = (localflowrate >= 0);
                %Define "adjelem" and "rownumber"
                adjelem = booleanfr*leftelem + (1 - booleanfr)*rightelem;
                %Define "elemeval"
                elemeval = booleanfr*[leftelem rightelem] + ...
                    (1 - booleanfr)*[rightelem leftelem];
                
                %--------------------------------
                %Get data from UPSTREAM half-edge
                
                booleanupst = (inode < nodesur && localflowrate >= 0) || ...
                    (inode > nodesur && localflowrate < 0);
                next = (i + 1)*(i < lengnsurn) + (i == lengnsurn);
                before = (i - 1)*(i > 1) + lengnsurn*(i == 1);
                
                %Define the real other vertex:
                othervtx = ...
                    nsurn(next*booleanupst + before*(1 - booleanupst));
                %Get the number of adjacent half-edge row 
                %(for "bedge" or "inedge")
                rownumber = localrowposit(next*booleanupst + ...
                    before*(1 - booleanupst));
                
                %Get the "inode" coordinate
                othervtxcoord = coord(othervtx,1:2);
                %Get a vector correnponding to edge evaluated
                upstrvec = othervtxcoord - inodecoord;

                if rownumber <= bedgesize
                    booleankey = (inode == bedge(rownumber,1));
                    %Define "isleft"
                    isleft = bedge(rownumber,3);
                    surelem = adjelem;
                else
                    booleankey = (inode == inedge(rownumber - bedgesize,1));
                    %It is used insteady "setdiff"
                    possiblelem = inedge(rownumber - bedgesize,3:4);
                    surelem = [possiblelem(logical(possiblelem ~= adjelem)) ...
                        adjelem];
                    %Define "isleft"
                    isleft = possiblelem(1);
                end  %End of IF
            
                %Get the "flowrate" position
                upstfrposition = 2*rownumber - booleankey;
                %Get the local flowrate
                upstrflowrate = flowrate(upstfrposition);
                %It rounds the flow rate value
                booleanfr = abs(upstrflowrate) >= tol;
                %Define again "localflowrate"
                upstrflowrate = booleanfr*upstrflowrate;

                %Calculate the weight:                
                mweight = calcmweight(adjelem,isleft,localflowrate,...
                    upstrflowrate,multdlimiter,evalvec,upstrvec);

                %auxiliary "verticescoord"
                auxvertcoord = [coord(inode,:); mean(verticescoord,1)];

            
                %Get the saturation in the:
                mlpbyelem = ...
                    mlplimiter((order == 1) + (order > 1)*elemeval(1),:);
                %1. adjacent half-edge:
%                 swonadjedge = getsatonedge(elemeval,vertices,verticescoord,...
%                     taylorterms,Sw,limiterflag,order,constraint,...
%                     flagknownvert,satonvertices,mlpbyelem,...
%                     centelem(elemeval,1:2)); 
                swonadjedge = getsatonedge(elemeval,vertices,auxvertcoord,...
                    taylorterms,Sw,limiterflag,order,constraint,flagknownvert,...
                    satonvertices,mlpbyelem,centelem(elemeval,1:2)); 
                %2. surounding half-edge:
                %Define the vertices of the distant half-edge
                vertdistant = [inode othervtx];
                %Get the coordinate of the distant half-edge
                vertdistantcoord = coord(vertdistant,:);
            
                %auxiliary "vertdistantcoord"
                auxvertdstcoord = [coord(inode,:); mean(vertdistantcoord,1)];
            
                mlpbyelem = ...
                    mlplimiter((order == 1) + (order > 1)*surelem(1),:);
                %Calculate the saturation on the distant half-edge.
%                 swonsuredge = getsatonedge(surelem,vertdistant,...
%                     vertdistantcoord,taylorterms,Sw,limiterflag,order,...
%                     constraint,flagknownvert,satonvertices,...
%                     mlpbyelem,centelem(surelem,1:2)); 
                swonsuredge = getsatonedge(surelem,vertdistant,...
                    auxvertdstcoord,taylorterms,Sw,limiterflag,order,...
                    constraint,flagknownvert,satonvertices,mlpbyelem,...
                    centelem(surelem,1:2)); 
                mlpbyelem = ...
                    mlplimiter((order == 1) + (order > 1)*elemeval(2),:);
                %3. Saturation on the another side of the half-edge 
                %(another state):
%                 swonanotherside = getsatonedge(fliplr(elemeval),vertices,...
%                     verticescoord,taylorterms,Sw,limiterflag,order,...
%                     constraint,flagknownvert,satonvertices,...
%                     mlpbyelem,centelem(fliplr(elemeval),1:2)); 
                swonanotherside = getsatonedge(fliplr(elemeval),vertices,...
                    auxvertcoord,taylorterms,Sw,limiterflag,order,...
                    constraint,flagknownvert,satonvertices,...
                    mlpbyelem,centelem(fliplr(elemeval),1:2)); 

                %Attribute to "satonadjsur" the saturation value calculated  
                %by the function "getsatonedge".
                satonadjsur = [swonadjedge swonsuredge];  

                %Store the position saturation in adjedge and suredge 
                %(in sequence)
                satonadjsuredge(m + 1:m + 2) = satonadjsur;
                %Store the saturation calculated by the another state.
                satonanotherstate(l) = swonanotherside;
                %Store the weight calculated.
                mweightvec(l) = mweight;
                %Store the position of the half-edge. It is used for store 
                %the saturation in each half-edge in a unic vector.
                halfedgepos(l) = frposition;

                %Turn "getvalkey" on
                getvalkey = 1;
            end  %End of IF (maybe well treatment)
        end  %End of IF (half-edges: internal or on the boundary)

        %Define the boolean for increment
        booleaninc = (getvalkey == 1);
        %Update the auxiliary counters "m", "l" and "vtxcount"
        m = m + 2*booleaninc;
        l = l + 1*booleaninc;
        vtxcount = vtxcount + 2*booleaninc;
        u = u + 1*(getvalkey == 1);
    end  %End of FOR (amount of half-edges surrounding each node)

    %Fill "amthe_well" (amount of half-edges surrounding the vertex 
    %associated to producer well(s))
    if wellprodkey == 1
        amthe_well(j) = u;
    end  %End of IF
end  %End of FOR (swept the vertices)
