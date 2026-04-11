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

function [satonadjsuredge,satonanotherstate,halfedgepos,amthe_well] = ...
    getmultidsatweig_implct(Sw,taylorterms,limiterflag,flowrate,mobility,...
    flagknownedge,satonedges,massweigmap,othervertexmap,multdlimiter,...
    constraint,flagknownvert,satonvertices,mlplimiter,coordmaprodelem,...
    wellprodkey)
%Define global parameters:
global coord centelem inedge bedge order rowposit esurn1 esurn2 nsurn1 ...
    nsurn2 timelevel;

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
satonadjsuredge = zeros(length(nsurn1),1);
%"satonanotherstate" stores the state on the other side of the half-edge
%evaluated. We verify the side where the state is MultiD (store it in 
%"satonadjsuredge") and the another is stored into "satonanotherstate".
satonanotherstate = satonadjsuredge;   %zeros(length(nsurn1),1);
halfedgepos = satonanotherstate;

%Initialize the auxiliary counters "l" and "c" and "vtxcount"
c = 0;
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
    %"esurn"
    esurn = esurn1(esurn2(inode) + 1:esurn2(inode + 1));
    %"nsurn":
    nsurn = nsurn1(nsurn2(inode) + 1:nsurn2(inode + 1));
    %"localrowposit"
    localrowposit = rowposit(nsurn2(inode) + 1:nsurn2(inode + 1));
    
    %It catches the length of "esurn" and "nsurn"
    lengesurn = length(esurn);
    lengnsurn = length(nsurn);
    
    %Get the "inode" coordinate
    inodecoord = coord(inode,1:2);

    %Initialize the matrices "A" and "B"
    A = zeros(lengnsurn);  
    B = zeros(lengnsurn,lengesurn);
    %Initialize "knownsat"
    knownsat = zeros(lengnsurn,1);
    %Initialize "anotheunknw"
    anotheunknw = knownsat; 
    %Initialize "pointknown"    
    pointknown = knownsat;
    
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
%                 verticescoord = coord(vertices,:);
                %Define the control volume on the left.
                leftelem = bedge(localrowposit(i),3);
                %Define the adjacent element ("adjelem")
%                 adjelem = leftelem;

                %--------------------------------
                %Flow rate in UPSTREAM half-edge

                %Calculate the weight:                
                mweight = 0;

                %Fill "A" and "B"
                %Fill "A".
                A(i,i) = 1;   
                %Fill "B"
                %Points to posit. in "esurn" of element "leftelem"
                pointposes = logical(esurn == leftelem); 
                %Attribute the term "(1 - mweight)" to corrected position
                B(i,pointposes) = (1 - mweight);

                %Verify if there exists a boundary cond. (prescribed 
                %saturation).
                %Obs.: It is used in Buckley-Leverett applications.
                if flagknownedge(localrowposit(i)) == 1;
                    knownsat(i) = satonedges(localrowposit(i)); 
                    %Fill "pointknown"
                    pointknown(i) = 1;
                end  %End of IF

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
            
%             mobtot_adjaux = mobility(localrowposit(i),:);
%             mobw_adj = mobtot_adjaux(1);
%             mobtot_adj = sum(mobtot_adjaux);
            fw_adj = 1;%mobw_adj/(mobtot_adj + 1e-16);
            
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
                %Get the number of upstream half-edge row 
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
                else
                    booleankey = (inode == inedge(rownumber - bedgesize,1));
                    %It is used insteady "setdiff"
                    possiblelem = inedge(rownumber - bedgesize,3:4);
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

%                 mobtot_upstaux = mobility(rownumber,:);
%                 mobw_upst = mobtot_upstaux(1);
%                 mobtot_upst = sum(mobtot_upstaux);
                fw_upst = 1;%mobw_upst/(mobtot_upst + 1e-16);

                
                %Calculate the weight:                
                mweight = calcmweight(adjelem,isleft,localflowrate,...
                    upstrflowrate,multdlimiter,evalvec,upstrvec);

                %auxiliary "verticescoord"
%                 auxvertcoord = [coord(inode,:); mean(verticescoord,1)];


                %Fill the matrix ref. to sat. on half-edges.
                %Term associated to Sk (see Kozdon et al., 2011)
                %Fill "A":   
                A(i,i) = 1;   
                %Points to posit. in "nsurn" of vector elem "otherhalfedge"
                pointposns = logical(nsurn == othervtx);
                %Attribute the term "-mweight" to corrected position
                A(i,pointposns) = -mweight;

                %Fill "B"
                %Points to posit. in "esurn" of element "leftelem"
                pointposes = logical(esurn == adjelem); 
                %Attribute the term "(1 - mweight)" to corrected position
                B(i,pointposes) = (1 - mweight);
                
                %Get the saturation in the:
                mlpbyelem = ...
                    mlplimiter((order == 1) + (order > 1)*elemeval(1),:);
                %1. adjacent half-edge:
%                 swonadjedge = getsatonedge(elemeval,vertices,verticescoord,...
%                     taylorterms,Sw,limiterflag,order,constraint,...
%                     flagknownvert,satonvertices,mlpbyelem,...
%                     centelem(elemeval,1:2)); 
%                 swonadjedge = getsatonedge(elemeval,vertices,auxvertcoord,...
%                     taylorterms,Sw,limiterflag,order,constraint,flagknownvert,...
%                     satonvertices,mlpbyelem,centelem(elemeval,1:2)); 

                mlpbyelem = ...
                    mlplimiter((order == 1) + (order > 1)*elemeval(2),:);
                %2. Saturation on the another side of the half-edge 
                %(another state):
                swonanotherside = getsatonedge(fliplr(elemeval),vertices,...
                    verticescoord,taylorterms,Sw,limiterflag,order,...
                    constraint,flagknownvert,satonvertices,...
                    mlpbyelem,centelem(fliplr(elemeval),1:2)); 
%                 swonanotherside = getsatonedge(fliplr(elemeval),vertices,...
%                     auxvertcoord,taylorterms,Sw,limiterflag,order,...
%                     constraint,flagknownvert,satonvertices,...
%                     mlpbyelem,centelem(fliplr(elemeval),1:2); 

                %Attribute to "satonadjsur" the saturation value calculated  
                %by the function "getsatonedge".
                anotheunknw(i) = swonanotherside; 

                %Store the position saturation in adjedge and suredge 
                %(in sequence)
%                 satonadjsuredge(m + 1:m + 2) = satonadjsur;
%                 %Store the saturation calculated by the another state.
%                 satonanotherstate(l) = swonanotherside;

                %Store the position of the half-edge. It is used for store 
                %the saturation in each half-edge in a unic vector.
                halfedgepos(l) = frposition;

                %Turn "getvalkey" on
                getvalkey = 1;
            end  %End of IF (maybe well treatment)
        end  %End of IF (half-edges: internal or on the boundary)

        %Define the boolean for increment
        booleaninc = (getvalkey == 1);
        %Update the auxiliary counters "c", "m", "l" and "vtxcount"
        l = l + 1*booleaninc;
        vtxcount = vtxcount + 2*booleaninc;
        u = u + 1*(getvalkey == 1);
    end  %End of FOR (amount of half-edges surrounding each node)

    %Treat the know saturation values (if they exist)
    pointknown = logical(pointknown == 1);
    %Null the columns of known saturation (which corresponds to "i"th row)
    A(pointknown,:) = 0;
    %Put "1" in the known diagonal terms.
    A(logical(diag(pointknown))) = 1;
    %Null the columns of known saturation (which corresponds to "i"th row)
    B(pointknown,:) = 0;
    
%     %Get only the position that there is no known value
%     unknownsat = unknownsat(logical(pointknown == 1));

    %Solve the local algebraic system
    vecsol = A\B;
    %Calculate "swinhalfedge"
    satinhalfedge = vecsol*Sw(esurn) + knownsat;

    %Stores this little vector
    satonadjsuredge(c + 1:c + lengnsurn) = satinhalfedge;
    %Stores the saturation with the another state (non-MultiD)
    satonanotherstate(c + 1:c + lengnsurn) = anotheunknw; 
    
    %Update the auxiliary counters
    c = c + lengnsurn;
    
    %Fill "amthe_well" (amount of half-edges surrounding the vertex 
    %associated to producer well(s))
    if wellprodkey == 1
        amthe_well(j) = u;
    end  %End of IF
    
%     if timelevel == 26 && (j == 452 || j == 201 || j == 454 || j == 457)
%         j
%         nsurn
%         bedgesize
%         localrowposit
%         w'
%         a'
%         e'
%         pause
%     end
        
end  %End of FOR (swept the vertices)
