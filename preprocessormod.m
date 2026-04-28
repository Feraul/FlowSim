
%--------------------------------------------------------------------------

function env = preprocessormod(isimu)

%--------------------------------------------------------------------------
%It reads the data file (from "Start.dat")
[config] = getdatafile;
%==================================================
% VERIFICAÇÃO DE MESH
%==================================================
mainpathfolders = build_directories(config,isimu);

verifymshfile = fullfile(config.keypathmesh,config.mesh);
%--------------------------------------------------------------------------
%"coord" matrix - The coordinates (x,y,z) of each point in a carts.
[coord,nnode] = getcoord(verifymshfile);
%Gives the information of "coord" generated
disp('"coord" was generated!');
%--------------------------------------------------------------------------
%"elem" matrix - The points that constitute each element
[elem,nbe,nelem,nodelim,flaglim] = ...
    getelem(verifymshfile,nnode);
%Gives the information of "elem" generated
disp('"elem" was generated!');
%--------------------------------------------------------------------------
%"centelem" matrix - It get the centroid coordinate for each control volume

%Fill the matrix "centelemcoord"
centelem = getcentelem(coord,elem);

%Gives the information of "elem" generated
disp('"centelem" was generated!');
%--------------------------------------------------------------------------
%Fill the vector "elemarea"
elemarea = calcelemarea(coord,elem);
%Gives the information of "elem" generated
disp('"elemarea" was generated!');
%--------------------------------------------------------------------------
%"bedge" and "inedge" - It get inform. about boundary and internal edges
[bedge,inedge,klb,auxpar] = getinfoedge(coord, elem, nnode, nbe, nelem, ...
    nodelim, flaglim, verifymshfile, config.bcflag);

%[inedge,bedge] = reordbedge(inedge,bedge,coord,centelem);
%Gives the information of "bedge" and "inedge" generated
disp('"bedge" was generated!');
disp('"inedge" was generated!');
%Get the normals
[normals,unitnormals,medEdge] = calcnormals(coord,centelem,bedge,inedge);
%Gives the information of "normals" generated
disp('"normals" was generated!');

%[node2edge, node2elem] = build_connectivities(bedge, inedge, elem, size(coord,1));

%--------------------------------------------------------------------------
% ============================================================
% 1) CONTAR QUANTOS ELEMENTOS TOCAM CADA NÓ (esurnqnt)
% ============================================================

% ============================================================
% 2) CONSTRUIR ESURN2 (ponteiros CSR começando em zero)
% ============================================================
% esurn2 = zeros(nnode + 1,1);   
% 
%     %"ielem" is a counter of elements (swept the rows of "elem")
%     for ielem = 1:nelem
%         %Get "elemcontents"
%         elemcontents = elem(ielem,1:4);
%         %Get the amount of edges by element.
%         amountedges = sum(elemcontents ~= 0);
%         %Get a vector with the vertices for the element evaluated. 
%         vertices = elemcontents(1:amountedges);
%         %Get the number of vertices for this element.
%         numvert = length(vertices);
%         %This loop swept all nodes who constitute each element (columns of 
%         %"elem" matrix). 
%         for inode = 1:numvert    
%             %"nodeval" is the receives the number of each node. 
%             %The elem(1,1) has nodeval = 1, elem(1,2) has nodeval = 5 
%             %(see basic exemple in Marcio's notebook) 
%             nodeval = vertices(inode);
%             %Add one to position which corresponds to number of node. We
%             %have to consider that "esurn2" start with value 0 for a better
%             %use as will be saw in the future (see Lohner, 2001 - Chap 2)
%             esurn2(nodeval + 1) = esurn2(nodeval + 1) + 1;
%         end  %End of FOR (in each element)
%     end  %End of FOR (in all "elem" matrix)

    % Inicialización
esurn2 = zeros(nnode + 1, 1);

% Extraer todos los nodos (columnas 1:4)
elemcontents = elem(:, 1:4);

% Máscara de nodos válidos
mask = elemcontents ~= 0;

% Recorrer solo nodos válidos (vector plano)
nodes = elemcontents(mask);

% ÚNICO bucle (equivalente al doble for)
for k = 1:numel(nodes)
    nodeval = nodes(k);
    esurn2(nodeval + 1) = esurn2(nodeval + 1) + 1;
end
%"esurnqnt" recives "esurn2" before this be recalculated in order to 
%represent an acumulate value (function "cumsum") 
esurnqnt = esurn2;    
%Create a vector with the amount accumulated from "esurn2" already 
%obtained (esurn2cum(i) = esurn2(i) + esurn2(i-1))
esurn2 = cumsum(esurn2);
esurn2aux = esurn2;

%Create and initialize "esurn1"
%The number of rows of "esurn1" is equal to amount accumulated of elements
%in each node
esurn1 = zeros(esurn2(nnode + 1),1);
%Initialize "nodestore" which control if the node evaluated 
%already was verifyed. If the node in the "elem" was already evaluated 
%inode is added in one.
nodestore = 0;
%"inodestore" is a counter of "nodestore"
inodestore = 1;

%Fill "esurn1" (list of elements surrounding each node)
    %Swept from 1 until the last element
    for ielem = 1:nelem  
        
        %Get "elemcontents"
        elemcontents = elem(ielem,1:4);
        %Get the amount of edges by element.
        amountedges = sum(elemcontents ~= 0);
        %Get a vector with the vertices for the element evaluated. 
        vertices = elemcontents(1:amountedges);
        %Get the number of vertices for this element.
        numvert = length(vertices);
        %Swept each node of row's "elem" evaluated
        for inode = 1:numvert   
            %"nodeval" receives the index of node on "elem" matrix
            nodeval = vertices(inode);
            %Just one element concors to node evaluated.
            %In general this occur in the domain corner.
            if esurnqnt(nodeval + 1) == 1
                %"store" obtain the value (of "esurn2aux" vector) which 
                %corresponds to position whose number is "nodeval" + 1.
                %Remenbering "esurn2" begins from 0 (for that +1)
                store = esurn2aux(nodeval + 1);
                %Put the value which corresponds to "Elem" row in the "esurn1"
                %whose the position is termined by "strore" index
                esurn1(store) = ielem;
                %Finaly, we have to diminish 1 from value of "esurn2aux" 
                %on the index evaluated. This is done in order to avoid use 
                %again the same index. 
                %(see basic example in the Marcio's notebook)
                esurn2aux(nodeval + 1) = esurn2aux(nodeval + 1) - 1;
            %Two or more elements concor to node evaluated and the node 
            %evaluated to be inside domain. In general is what happen in 
            %the domain.
            elseif esurnqnt(nodeval + 1) > 1 && ...
                    ismember(nodeval,nodestore) ~= 1
                %Initialize "elemorder". This parameter receives the
                %elements which try to node evaluated in increasing order.
                elemorder = zeros(1,esurnqnt(nodeval + 1));
                
                %----------------------------------------------------------
                %"nodeval" over boundary
                
                %Verify if the node evaluated ("nodeval") is over the 
                %boundary. In this case we have to find the first element
                %and the last element in order maintain the counter
                %clockwise way (little more complicated).
                rowsinbedge = ...
                    bedge(any(logical(bedge(:,1:2) == nodeval),2),1:3);
                %The node evaluated is inside domain. In this case is
                %easer than the earlier situation once any element may be 
                %the first one. The cycle finishes when the last element 
                %before the first is found.
                if isempty(rowsinbedge)  
                    %Attribute the first element to "elemorder". 
                    %Considering "ielem" contain one of elements evaluated.
                    elemorder(1) = ielem;

                %The node evaluated is over the boundary. The first element
                %in "elemorder" cames from "bedge"
                else
                    %Initialize "vertexout" and "firstelem"
                    vertexout = 0;
                    firstelem = 0;
                    %Initialize "elemcandidate". It stores the number of
                    %elements which can be the first element in "elemorder"
                    elemcandidate = rowsinbedge(:,3);
                    %"rowsvertices" receives just the number of vertices
                    rowsvertices = rowsinbedge(:,1:2);
                    %We need verify which other vertices belong to edge
                    %over the boundary.
                    vertoutnodeval = setdiff(reshape(rowsvertices',4,1),...
                        nodeval,'stable');
                    
                    %Evaluate which edge is the first:
                    %Get the vector (first edge)
                    edgevec = coord(vertoutnodeval,:) - ...
                        [coord(nodeval,:); coord(nodeval,:)];
                    %Get a vector which points to out of domain.
                    pointout = [coord(nodeval,:); coord(nodeval,:)] - ...
                        centelem(elemcandidate,:);
                    
                    %Attribute the other vertex (edge over the boundary).
                    vertexout = vertexout + vertoutnodeval(1)*...
                        (cross(pointout(1,:),edgevec(1,:)) > 0);
                    vertexout = vertexout + vertoutnodeval(2)*...
                        (cross(pointout(2,:),edgevec(2,:)) > 0);
                    %Attribute the first element sharing the edge over the 
                    %boundary.
                    firstelem = firstelem + elemcandidate(1)*...
                        (setdiff(cross(pointout(1,:),edgevec(1,:)),0) > 0);
                    firstelem = firstelem + elemcandidate(2)*...
                        (setdiff(cross(pointout(2,:),edgevec(2,:)),0) > 0);

                    %Attribute the first element to "elemorder(1)"
                    elemorder(1) = firstelem;
                end  %End of internal IF

                %Create a vector with number of "inedge" rows
                iinrow = 1:size(inedge,1);

                %Find the last element to be stored in the "elemorder"
                for iorder = 2:esurnqnt(nodeval + 1)
                   
                    %Initialize "rowpointer". This parameter must be
                    %initialized in order to avoid superposition of values.
                    rowpointer = 0;
                    %Obtain the row or rows of "inedge" whose "nodeval" 
                    %belongs to first two columns of "inedge" and, to same 
                    %time, the aforefind element belongs to two last 
                    %columns of "inedge".
                    %Evaluate the vertices and elements in "inedge":
                    
                    evalvertelem = all([any(ismember(inedge(:,1:2),...
                        nodeval),2) any(ismember(inedge(:,3:4),...
                        elemorder(iorder-1)),2)],2);
                    %Get the "inedge" row.
                    rowinedge = iinrow(evalvertelem)';

                    %This loop points to entity of "rowinedge" which
                    %satisfy the necessary conditions. That row which
                    %satisfies receive the value "1", while that one 
                    %which does not satisfy receive "0"
                    for irow = 1:length(rowinedge)
                        rowpointer(irow) = ...
                            (inedge(rowinedge(irow),2) == nodeval & ...
                            elemorder(iorder-1) == inedge(rowinedge(irow),3)) | ...
                            (inedge(rowinedge(irow),1) == nodeval & ...
                            elemorder(iorder-1) == inedge(rowinedge(irow),4));
                    end  %End of internal FOR

                    %"rowdef" is a definitive row which satisfy the
                    %necessary conditions aforesaid. That is the 
                    %position which have "1" insteady "0" 
                    rowdef = logical(rowpointer == 1);   
                    %"rowinedge" receives, finaly, the unic row which 
                    %satisfy the conditions above.
                    rowinedge = rowinedge(rowdef);
                    %A decision is made as function of "nodeval"'s 
                    %value. The commands written below avoids "IF"  
                    %"nodeval" is lower than other node which 
                    %constitute the common edge among two elements. 
                    %In this case, a left element is attributed to 
                    %"elemorder" just if "nodeval" is lower than other 
                    %element   
                    elemorder(iorder) = elemorder(iorder) + ...
                        inedge(rowinedge,3)*any((nodeval ~= ...
                        max(inedge(rowinedge,1:2))));
                    %In this case, a right element is attributed to 
                    %"elemorder" just if "nodeval" is major than other 
                    %element   
                    elemorder(iorder) = elemorder(iorder) + ...
                        inedge(rowinedge,4)*any((nodeval == ...
                        max(inedge(rowinedge,1:2))));
                end  %End of FOR
                
                %"nodestore" receives the node evaluated. Thus, when this
                %will evaluated in the future that will be by passed,
                %because was already (see line 1241)
                nodestore(inodestore) = nodeval;
                %Add 1 to "inodestore"
                inodestore = inodestore + 1;

                %Put the number in its places in the "esurn1"
                store = esurn2aux(nodeval + 1);
                %Put the value which corresponds to "elem" row in the 
                %"esurn1" whose the position is termined by "strore" index
                
                %"retrocount" is a retroative counter.
                retrocount = esurnqnt(nodeval + 1);
                for iretro = store:-1:store - (esurnqnt(nodeval+1) - 1)
                    esurn1(iretro) = elemorder(retrocount);
                    %backward the counter
                    retrocount = retrocount - 1;
                end  %End of FOR

                %Finaly, we have to diminish 1 from value of "esurn2aux" 
                %on the index evaluated. This is done in order to avoid use 
                %again the same index. 
                %(see basic example in the Marcio's notebook)
                esurn2aux(nodeval + 1) = esurn2aux(nodeval + 1) - ...
                    esurnqnt(nodeval + 1);
            end  %End of external IF
            %Jump a node position if any conditions above is atended
            inode = inode + 1;  %Next node
        end  %End of WHILE (read each column of each row of "elem")
    end  %End of FOR (read each row of "elem")

disp('"esurn" was generated!');
%--------------------------------------------------------------------------
%"nsurn" vectors - Report the amount of NODES surronding each NODE.
% ============================================================
% Construcción de nsurn1 y nsurn2 (CSR)
% ============================================================

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%"nsurn" vectors - Report the amount of NODES surronding each NODE.

%We will use the vectors "eserp1" and "esurn2" to swept the elements around
%each point. In the evaluated element we must find the near points to point
%who corresponds to position of "Esurn2" (the first position of this vector 
%corresponds to point 1 and consecutively)

%Initialization of variables
%Works as a counter to "nsurn1" vector
store = 0;
%"nsurn2" points to "nsurn1" and denote the begin and end of node 
%surrounding cycle of each node
nsurn2 = zeros(nnode + 1,1);

%Swept the vector "esurn2" to pass for each point (first for). Each one 
%index "i" is a position of aforesaid (abovementioned) vector.
%With the development of "i", points already visited will not any more  
    for i = 1:nnode
        %Swept from position "i" add 1 until position "i+1" (end of 
        %elements surrounding point cycle).
        %"iesurn" is a counter of thous elements surrounding each point

        %Initialization of variable
        %"nsurn1aux" is an auxiliar vector which avoid repeated node 
        %be stored in "nsurn1" (see below)
        nsurn1aux = 0;
        %"insurn" is a counter who will be used in "nsurn1aux" vector
        insurn = 0;
        %Swept all elements surround each node evaluated
        for iesurn = esurn2(i) + 1:esurn2(i + 1)
            %"ielem" will count which elements be in the cycle.
            %This parameter receives the elements 1,2 & 3 (whose positions 
            %are 7,8 and 9) - see basic example in Marcio's notebook
            ielem = esurn1(iesurn);   
            
            %Initializing "inode" (node counter in each element evaluated)
            inode = 1;
                %Verify in each element (use the "elem" matrix) which
                %nodes are linked with the one evaluated.
                %This happen while each node evaluated is different of zero 
                %and if that node does not belong to 5th column of "elem" 
                %matrix 
                while inode < 5 && elem(ielem,inode) > 0 
                    %The decision command store the value from "elem" in 
                    %"nsurn1" if and only if the node evaluated does not be
                    %the node "i" (becouse we are searching the neighbour 
                    %nodes to node "i", not the own node "i"). Another
                    %condition is that the neighbour node already stored
                    %does not be stored again
                    %This decision command is exclusive to TRIANGLES
                    if elem(ielem,4) == 0 & elem(ielem,inode) ~= i & ...
                            elem(ielem,inode) ~= nsurn1aux 
                        %"store" stores the position in "Nsurn1" vector
                        store = store + 1;
                        %Add 1 to counter "insurn". When "i" changes its 
                        %value this parameter is nulled
                        insurn = insurn + 1;
                        %"Nsurn1" receives the node who attends to above
                        %conditions
                        nsurn1(store) = elem(ielem,inode);
                        %"Nsurn1aux" receives the nodes in order to those
                        %do not be repeated
                        nsurn1aux(insurn) = elem(ielem,inode);
                    %The next decisions command is exclusive to QUADRAN.
                    %The neighbour to node evaluated will be considerate if
                    %those are either forward or backward positioned with 
                    %regard to evaluated node.  
                    elseif elem(ielem,4) ~= 0 && elem(ielem,inode) ~= i
                        %A special case where the first node evaluated into 
                        %row "ielem" is different of value "i" but does not 
                        %is neighbour of node "i". This happen only when
                        %the node evaluated is the first and the node who
                        %corresponds to "i" is the third.
                        if inode == 1 && elem(ielem,3) == i
                            inode = inode + 1;
                        end  %End of first internal IF
                            %To avoid which node already accounted be 
                            %stored again this decision command just 
                            %storage nodes who does not be in Nsurn1aux 
                            %(auxiliar vector)
                            if elem(ielem,inode) ~= nsurn1aux
                                store = store + 1;
                                insurn = insurn + 1;
                                nsurn1(store) = elem(ielem,inode);
                                nsurn1aux(insurn) = elem(ielem,inode);
                                %Avoid which the neighbour node be stored, 
                                %once it does not is straightly linked with 
                                %node evaluated
                            end  %End of second internal IF
                            
                            %If the node evaluated was considered neighbour
                            % of "i" (and stored) the counter "inode"
                            %must jump one node because the next one
                            %certaly does not will neighbour 
                            inode = inode + 1;
                    end  %End of external IF
                    inode = inode + 1;  %Evaluate the next node
                end  %End of WHILE
        end  %End of second FOR
        %Update "nsurn2". This variable receives the number of increments
        %that "store" had
        nsurn2(i+1) = store;
    end  %End of first FOR
%Reorder "nsurn1" (counterclockwise). It returns a column vector.
nsurn1 = reordernsurn(esurn1,esurn2,nsurn1,nsurn2,bedge,coord,elem);
%---------------------------------------------------------------------------
%Gives the information of "nsurn" generated
disp('"nsurn" was generated!');

%--------------------------------------------------------------------------
%Get the position of "bedge" and "inedge" rows.

%It has a structure similar to "nsurn1"
%rowposit = getrowposition(bedge,inedge,nsurn1,nsurn2);

%--------------------------------------------------------------------------
%Generate "esureface" (face neighbor element) and "esurefull" (vertices and
%face) neighbor

%[esureface1,esureface2,esurefull1,esurefull2] = getesure(elem,...
%    size(bedge,1),inedge,esurn1,esurn2,nsurn2,rowposit);

%Gives the information of "nsurn" generated
%disp('"esure" were generated!');

geometry.coord        = coord;
geometry.centelem     = centelem;
geometry.elem         = elem;

geometry.esurn1       = esurn1;
geometry.esurn2       = esurn2;

geometry.nsurn1       = nsurn1;
geometry.nsurn2       = nsurn2;

geometry.bedge        = bedge;
geometry.inedge       = inedge;

geometry.normals      = normals;
geometry.medEdge      = medEdge;
geometry.unitnormals  = unitnormals;

%geometry.esureface1   = esureface1;
%geometry.esureface2   = esureface2;
%geometry.esurefull1   = esurefull1;
%geometry.esurefull2   = esurefull2;

geometry.elemarea     = elemarea;
geometry.klb          = klb;

env.config=config;
env.geometry=geometry;
env.mainpathfolders=mainpathfolders;
end
%% ========================================================================
%% ========================================================================
%% ========================================================================
%--------------------------------------------------------------------------
%Generating or acessing the folder


function mainpathfolders = build_directories(config,isimu)

%==================================================
% PATH BASE
%==================================================

pathstore = char(config.pathstore);
complemfile = char(config.complemfile);

% Criar pasta principal se não existir
if ~exist(pathstore,'dir')
    mkdir(pathstore);
end

%==================================================
% DEFINIR MAINPATH
%==================================================

if strcmp(complemfile,'0')

    path = pathstore;

    resfolder  = 'Results';
    datafolder = 'Data';
    tabfolder  = 'Tables';

else

    path = fullfile(pathstore, complemfile);

    if ~exist(path,'dir')
        mkdir(path);
    end
    isimu = num2str(isimu);


    resfolder  = ['Results_' complemfile '_' isimu];
    datafolder = ['Data_' complemfile];
    tabfolder  = ['Tables_' complemfile '_' isimu];

end


%==================================================
% CRIAR SUBPASTAS DE FORMA VETORIZADA
%==================================================

folders = {
    fullfile(path,resfolder)
    fullfile(path,datafolder)
    fullfile(path,tabfolder)};

for k = 1:length(folders)
    if ~exist(folders{k},'dir')
        mkdir(folders{k});
    end
end
mainpathfolders.path=path;
mainpathfolders.resfolder=resfolder;
mainpathfolders.datafolder=datafolder;
mainpathfolders.tabfolder=tabfolder;

end

%--------------------------------------------------------------------------
%FUNCTION "getdatafile"
%--------------------------------------------------------------------------
function config = getdatafile()

fid = fopen('Start.dat');

if fid == -1
    error('Cannot open Start.dat');
end

%% ==================================================
% LER ARQUIVO
%% ==================================================

raw = textscan(fid,'%s','Delimiter','\n','Whitespace','');
raw = raw{1};
fclose(fid);

%% ==================================================
% LIMPAR COMENTÁRIOS
%% ==================================================

raw = regexprep(raw,'#.*$','');
raw = regexprep(raw,'%.*$','');
raw = regexprep(raw,'//.*$','');

raw = strtrim(raw);
raw = raw(~cellfun('isempty',raw));

lines = raw;

%% ==================================================
% FUNÇÕES AUXILIARES
%% ==================================================

i = 1;

%getValue = @() getLine();

    function val = getValue()
    line = strtrim(lines{i});
    i = i + 1;

    % tenta converter para número
    num = str2double(line);

    if ~isnan(num)
        val = num;       % retorna número
    else
        val = line;      % retorna texto
    end
   end

    function val = getNum()
    line = strtrim(lines{i});
    i = i + 1;
    val = sscanf(line, '%f')';
end

    function [M,i] = nextMatrixAuto(lines,i)

        rows = {};
        n = 0;

        while i <= length(lines)

            nums = sscanf(lines{i},'%f')';

            % parar se não for linha numérica
            if isempty(nums)
                break
            end

            % garantir que tem pelo menos 2 colunas
            if length(nums) >= 2
                n = n + 1;
                rows{n,1} = nums(1:2);   % pega só 2 colunas
            end

            i = i + 1;

        end

        if n == 0
            M = [];
        else
            M = vertcat(rows{:});
        end

    end
%% ==================================================
% INICIO PARSER
%% ==================================================

config = struct();

config.benchkey = getValue();

config.pathstore = getValue();

% extrair numero do caso
digits = regexp(config.pathstore, '\d+(?:[._]\d+)?', 'match');

if ~isempty(digits)
    token = digits{end};          % pega o último número encontrado
    token = strrep(token, '_', '.');
    config.numcase = str2double(token);
else
    config.numcase = [];
end
config.complemfile = getValue();
config.phasekey = getValue();
config.modflowcase=getValue();
config.xyrz = getValue();
config.r0=getValue();

config.symaxe=getValue();
config.keymsfv=getValue();
config.coarseratio=getValue();
config.pmethod=getValue();
config.auxcvfactor=getValue();

config.interptype=getValue();
config.nltol=getValue();
config.maxiter=getValue();
config.acel=getValue();
config.smethod=getValue();
config.multdopt=getValue();
config.additionalmult=getValue();
config.goefreeopt=getValue();
config.order=getValue();
config.timeorder=getValue();
config.recovtype=getValue();
config.lsneightype=getValue();
config.lsexp=getValue();
config.limiterflag=getValue();
config.limitertype=getValue();
config.factorgradientmethod=getValue();

config.patchlimiter=getValue();
config.MultiDlimiter=getValue();
config.factorDelta=getValue();
config.MUSCL=getValue();
config.greekcorrection=getValue();
config.MOOD=getValue();
config.MOODstrategy=getValue();
config.MUSCLposLimited=getValue();
config.MOODeMUSCL=getValue();
config.LPMUSCL=getValue();
config.ReimannSolver=getValue();
config.keypathmesh=getValue();
config.mesh=getValue();
config.dens = getValue();
config.visc = getValue();
config.satirrectuble = getValue();
i=i+1;
config.perm = getNum();
config.pormap=getValue();
config.keygravity=getValue();
config.g=getValue();
config.keycapil=getValue();
config.ncaplcorey=getValue();
config.courant=getValue();
i=i+1;
config.totaltime=getValue();
bb=getValue();
[M]= nextMatrixAuto(lines(1:i+(bb-1)),i);
config.bcflag=M;
end

%--------------------------------------------------------------------------
%FUNCTION "getcoord"
%--------------------------------------------------------------------------

function [coord, nnode] = getcoord(verifymshfile)

%==================================================
% ABRIR ARQUIVO
%==================================================

fid = fopen(verifymshfile,'r');
if fid == -1
    error('Cannot open file: %s', verifymshfile);
end

%==================================================
% LER HEADER E NÚMERO DE NÓS
% (funciona para formato padrão Gmsh)
%==================================================

% Procurar seção $Nodes (mais robusto que HeaderLines fixo)
line = fgetl(fid);
while ischar(line)
    if strcmp(strtrim(line),'$Nodes')
        break
    end
    line = fgetl(fid);
end

% Ler número de nós
nnode = fscanf(fid,'%d',1);

%==================================================
% LER COORDENADAS
% Formato padrão:
% node_id x y z
% (para 2D, z pode não existir)
%==================================================

% Ler bloco completo
data = fscanf(fid,'%f');

fclose(fid);

% Número de colunas por nó (3D padrão)
% Detectar automaticamente dimensão

if mod(length(data),4) == 0
    % Formato: id + 3 coordenadas
    data = reshape(data,4,nnode);
    coord = data(2:4,:)';

elseif mod(length(data),3) == 0
    % Formato: id + 2 coordenadas (2D)
    data = reshape(data,3,nnode);
    coord = data(2:3,:)';

    % adicionar terceira coluna zero para compatibilidade
    coord(:,3) = 0;

else
    error('Formato inesperado do arquivo .msh');
end

end
%--------------------------------------------------------------------------
%FUNCTION "getelem"
%--------------------------------------------------------------------------

function [elem,nbe,nelem,nodelim,flaglim] = ...
    getelem(verifymshfile,nnode)

%=========================================================
% ABRIR ARQUIVO UMA ÚNICA VEZ
%=========================================================
fid = fopen(verifymshfile,'r');
if fid == -1
    error('Cannot open mesh file');
end

%% =========================================================
%  IR ATÉ $Entities  (como no seu código original)
% ==========================================================
getmshdata = textscan(fid,'%u',1,'HeaderLines',7 + nnode);
nent = getmshdata{1};

% Leitura direta como matriz numérica (rápido)
auxmat = textscan(fid,'%*n %n %*n %n %*[^\n]', nent);
auxmat = [auxmat{:}];

entitype = auxmat(:,1);   % tipo da entidade física (1,2,3,15,...)

%% =========================================================
%  ESTATÍSTICAS (mantido exatamente como no seu código)
% =========================================================
nbe     = sum(entitype == 1);     % entidades tipo linha
nbeaux  = nbe;
nodelim = sum(entitype == 15);    % entidades tipo ponto físico
intnode = 0;                      % será atualizado depois

ntri  = sum(entitype == 2);       % entidades tipo superfície tri
nquad = sum(entitype == 3);       % entidades tipo superfície quad

nelem = ntri + nquad;             % total de elementos 2D

% flags externas (mantido)
flaglim = auxmat(entitype == 15 & auxmat(:,2) < 1000, 2);

fclose(fid);

%% =========================================================
%  LER ELEMENTOS (vetorizado)
% =========================================================
fid = fopen(verifymshfile);
C = textscan(fid,'%s','Delimiter','\n','Whitespace','');
fclose(fid);
C = C{1};

% Localiza $Elements
i1 = find(contains(C,'$Elements'),1);
nelem_file = str2double(strtrim(C{i1+1}));

% Linhas dos elementos
elem_lines = C(i1+2 : i1+1+nelem_file);

% Extrai todos os números de cada linha
tokens = regexp(elem_lines,'\d+','match');

% Tipo do elemento = segundo número
types = cellfun(@(t) str2double(t{2}), tokens);

%% =========================================================
%  TRIÂNGULOS (type = 2)
% =========================================================
tri_tokens = tokens(types == 2);
tri = cellfun(@(t) str2double(t(end-2:end)), tri_tokens, 'UniformOutput', false);
tri = vertcat(tri{:});   % ntri × 3

%% =========================================================
%  QUADRILÁTEROS (type = 3)
% =========================================================
quad_tokens = tokens(types == 3);
quad = cellfun(@(t) str2double(t(end-3:end)), quad_tokens, 'UniformOutput', false);
quad = vertcat(quad{:}); % nquad × 4

%% =========================================================
%  MATRIZ elem FINAL
% =========================================================
ntri  = size(tri,1);
nquad = size(quad,1);

elem = zeros(ntri+nquad,5);   % 5 colunas como no seu código original

elem(1:ntri,1:3) = tri;       % triângulos
elem(ntri+1:end,1:4) = quad;  % quadriláteros

%% =========================================================
%  GARANTIR MATERIAL != 0
% =========================================================
elem(:,5) = 1:size(elem,1);
end
%--------------------------------------------------------------------------
%FUNCTION "getinfoedge"
%--------------------------------------------------------------------------
function [bedge, inedge, klb,auxpar] = getinfoedge(coord, elem, nnode, nbe, nelem, ...
    nodelim, flaglim, verifymshfile, bcflag)

%--------------------------------------------------------------
% 1) Ler dados de borda e montar bedgeaux
%--------------------------------------------------------------
bedgeaux = zeros(nbe,5);
inedge   = zeros(0,4);   % será preenchida depois

% Abrir arquivo .msh
readbound = fopen(verifymshfile);
intnode   = 0;

% Ler: flag de contorno + dois nós da aresta
getboundata = textscan(readbound, '%*u %*u %*u %u %*u %u %u', nbe, ...
    'HeaderLines', 8 + nnode + nodelim + intnode);

fclose(readbound);

bedgeaux(:,5)   = getboundata{1};          % flag geométrica / BC
bedgeaux(:,1:2) = [getboundata{2} getboundata{3}];

%--------------------------------------------------------------
% 2) Filtrar bedge pelas flags válidas em bcflag
%--------------------------------------------------------------
%--------------------------------------------------------------
% 2) Filtrar bedge pelas flags válidas em bcflag
%--------------------------------------------------------------
pntcorrectrow = ismember(bedgeaux(:,5), bcflag(:,1));
bedge = bedgeaux(pntcorrectrow,:);

%--------------------------------------------------------------
% 3) Garantir orientação anti-horária (CCW)
%--------------------------------------------------------------
P1 = coord(bedge(:,1), :);
P2 = coord(bedge(:,2), :);
C  = mean(coord, 1);

v = P2 - P1;      % vetor da aresta
r = P1 - C;       % vetor do centro até o nó inicial

% Produto vetorial 2D: cross(r, v)
cross2D = r(:,1).*v(:,2) - r(:,2).*v(:,1);

flipmask = cross2D < 0;
bedge(flipmask,1:2) = bedge(flipmask,[2 1]);



%--------------------------------------------------------------
% 3) Encontrar linhas com Neumann não nulo → klb
%--------------------------------------------------------------
pointflag = bcflag(:,1) > 200 & bcflag(:,2) ~= 0;

if any(pointflag)
    flagval = bcflag(pointflag,1);
    mask    = ismember(bedge(:,5), flagval);
    klb     = find(mask).';
else
    klb = [];
end

%--------------------------------------------------------------
% 4) Ajustar 4ª coluna de bedge com flaglim
%--------------------------------------------------------------
bedge(:,4) = bedge(:,5);   % inicialmente copia a 5ª

[tf, loc] = ismember(bedge(:,1), 1:nodelim);
bedge(tf,4) = flaglim(loc(tf));

%--------------------------------------------------------------
% 5) Construir arestas de todos os elementos
%    (assumindo até 4 nós por elemento, zeros = inexistentes)
%--------------------------------------------------------------
if max(elem(:,4))~=0
    nodes = elem(:,1:4);

    % pares (1-2, 2-3, 3-4, 4-1)
    tmp = [ ...
        nodes(:,1) nodes(:,2) ...  % a b
        nodes(:,2) nodes(:,3) ...  % b c
        nodes(:,3) nodes(:,4) ...  % c d
        nodes(:,4) nodes(:,1) ...  % d a
        ];
    auxpar=4;
else
    nodes = elem(:,1:3);

    % pares (1-2, 2-3, 3-4, 4-1)
    tmp = [ ...
        nodes(:,1) nodes(:,2) ...  % a b
        nodes(:,2) nodes(:,3) ...  % b c
        nodes(:,3) nodes(:,1) ...  % c d
        ];
    auxpar=3;
end

edges_all = reshape(tmp.', 2, []).';



%edges_all = [e1; e2; e3; e4];
elem_all  = repelem((1:nelem)',auxpar);

% remover arestas com nó zero ou repetido
invalid = edges_all(:,1)==0 | edges_all(:,2)==0 | edges_all(:,1)==edges_all(:,2);
edges_all(invalid,:) = [];
elem_all(invalid)    = [];

%--------------------------------------------------------------
% 6) Associar arestas de bedge aos elementos (preencher col. 3)
%   (assumindo mesma orientação de nós entre elem e bedge)
%--------------------------------------------------------------
% Ordenar arestas do contorno
edges_bedge_sorted = sort(bedge(:,1:2), 2);

% Ordenar arestas internas
edges_all_sorted = sort(edges_all, 2);

% Comparar sempre arestas ordenadas
[tfB, locB] = ismember(edges_bedge_sorted, edges_all_sorted, 'rows');

% Preencher coluna 3 com o elemento vizinho
bedge(tfB,3) = elem_all(locB(tfB));



%--------------------------------------------------------------
% 7) Arestas internas: aquelas de edges_all que não estão em bedge
%   (também assumindo mesma orientação)
%--------------------------------------------------------------
[tfInt, ~] = ismember(edges_all, edges_bedge_sorted, 'rows');
internalEdges = edges_all(~tfInt,:);
internalElem  = elem_all(~tfInt);

% ============================================================
% 1) Normalizar arestas internas (sem sort)
% ============================================================
internalNorm = [ ...
    min(internalEdges(:,1), internalEdges(:,2)), ...
    max(internalEdges(:,1), internalEdges(:,2)) ...
    ];

[U, ~, ic] = unique(internalNorm,'rows','stable');
counts = accumarray(ic,1);
idxTwo = find(counts == 2);

inedge = zeros(numel(idxTwo),4);
inedge(:,1:2) = U(idxTwo,:);

% ============================================================
% 2) Encontrar elementos vizinhos (loop NECESSÁRIO)
% ============================================================

A = inedge(:,1);
B = inedge(:,2);

leftElem  = zeros(numel(idxTwo),1);
rightElem = zeros(numel(idxTwo),1);

% lista de pares (nó, elemento)
[row, ~] = find(elem(:,1:auxpar) ~= 0);
nodeList = elem(elem(:,1:auxpar) ~= 0);
elemList = row;

% node2elem{k} = lista de elementos que usam o nó k
node2elem = accumarray(nodeList, elemList, [nnode 1], @(x){x});

for k = 1:numel(idxTwo)

    a = A(k);
    b = B(k);

    % encontrar elementos que contenham ambos os nós
    %elemsA = any(elem == a, 2);
    %elemsB = any(elem == b, 2);
    elems = intersect(node2elem{a}, node2elem{b});

    %elems = find(elemsA & elemsB);   % deve ter exatamente 2 elementos

    if numel(elems) ~= 2
        error('Aresta interna não encontrada em exatamente dois elementos.');
    end

    el1 = elems(1);
    el2 = elems(2);

    % ============================================================
    % 3) Determinar left/right via produto vetorial
    % ============================================================

    p1 = coord(a,:);
    p2 = coord(b,:);
    v1 = p2 - p1;

    % terceiro nó do primeiro elemento
    nodes1 = elem(el1, elem(el1,:)~=0);
    nodes1(nodes1==a | nodes1==b) = [];
    p3L = coord(nodes1(1),:);

    % terceiro nó do segundo elemento
    nodes2 = elem(el2, elem(el2,:)~=0);
    nodes2(nodes2==a | nodes2==b) = [];
    p3R = coord(nodes2(1),:);

    vL = p3L - p1;
    crossL = cross(v1, vL);

    if crossL(3) >= 0
        leftElem(k)  = el1;
        rightElem(k) = el2;
    else
        leftElem(k)  = el2;
        rightElem(k) = el1;
    end
end

inedge(:,3) = leftElem;
inedge(:,4) = rightElem;
end
%----------------------------
%Function "calcnormals"
%--------------------------------------------------------------------------

function [normals,unitnormals,medEdge] = calcnormals(coord, centelem, bedge, inedge)

%% ============================================================
% 1) MATRIZ DE ROTAÇÃO
% =============================================================
R = [0 1 0; -1 0 0; 0 0 0];

%% ============================================================
% 2) NORMAIS DAS ARESTAS DE FRONTEIRA (bedge)
% =============================================================

nb = size(bedge,1);

% Vetores das arestas
vB = coord(bedge(:,2),:) - coord(bedge(:,1),:);
medB=0.5.*(coord(bedge(:,2),:) + coord(bedge(:,1),:));
% Normais rotacionadas
normB = (R * vB')';

% Ponteiro para verificar orientação
pointer = coord(bedge(:,1),:) - centelem(bedge(:,3),:);

% Correção de orientação
signFix = sign(sum(pointer .* normB, 2));
normB = normB .* signFix;


%% ============================================================
% 3) NORMAIS DAS ARESTAS INTERNAS (inedge)
% =============================================================

ni = size(inedge,1);

vI = coord(inedge(:,2),:) - coord(inedge(:,1),:);
normI = (R * vI')';
medI=0.5.*(coord(inedge(:,2),:) + coord(inedge(:,1),:));
%% ============================================================
% 4) CONCATENAÇÃO FINAL
% =============================================================

normals = [normB; normI];
unitnormals = [normB./vecnorm(normB,2,2); normI./vecnorm(normI,2,2)];
medEdge=[medB;medI];
end
%--------------------------------------------------------------------------
%FUNCTION "getsurnode"
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%FUNCTION "shiftchoosen"
%--------------------------------------------------------------------------

function [vecout] = shiftchoosen(vecin,numelem,letter)
%Initialize "vecout". It is the vector reordered according vector element
%choosen.
vecout(1:length(vecin),1) = vecin;

%It finds the position of element choosen in vector "vecin"
if strcmp(letter,'val')
    %Define "i" (auxiliary conuter 1:length(vecin))
    i = 1:length(vecin);
    %Poits its position
    pointpos = i(logical(vecin == numelem));
    %"pointpos" receives the position
else
    pointpos = numelem;
end  %End of IF

%It uses "circshift" to reorder the vector
if pointpos > 1
    %From "pointpos" on
    auxvecout = vecout(pointpos:length(vecout));
    %From "pointpos" on
    auxvecout(length(auxvecout) + 1:length(vecout)) = ...
        vecout(1:pointpos - 1);
    %Update "vecout"
    vecout = auxvecout;
end  %End of IF
end
%--------------------------------------------------------------------------
%FUNCTION "getcentelem"
%--------------------------------------------------------------------------
function [centelem] = getcentelem(coord,elem)

nelem = size(elem,1);
dim   = size(coord,2);

% Apenas nós válidos
mask = elem(:,1:4) ~= 0;

% Expandir coordenadas
nodes = elem(:,1:4);
nodes(~mask) = 1;

% Somar usando indexação vetorizada
x = coord(nodes(:),:);

x = reshape(x, nelem, 4, dim);

sumcoords = squeeze(sum(x,2));

% Número real de nós
numnodes = sum(mask,2);

% Centróide
centelem = sumcoords ./ numnodes;

end
%--------------------------------------------------------------------------
%Function "calcelemarea"
%--------------------------------------------------------------------------

%This function calculate the area of each element. The element is recognaze
%due parameter "numelem"
function elemarea = calcelemarea(coord, elem)

% Número de elementos
ne = size(elem,1);

% Inicializa vetor de áreas
elemarea = zeros(ne,1);

% --- Identifica elementos triangulares e quadrangulares ---
isQuad = elem(:,4) > 0;   % quadrangulares
isTri  = elem(:,4) == 0;  % triangulares

%% ============================
%   ÁREA DOS TRIÂNGULOS
% =============================
tri = elem(isTri, 1:3);

if any(isTri)
    A = coord(tri(:,1),:);
    B = coord(tri(:,2),:);
    C = coord(tri(:,3),:);

    v1 = B - A;
    v2 = B - C;

    elemarea(isTri) = 0.5 * vecnorm(cross(v1, v2, 2), 2, 2);
end

%% ============================
%   ÁREA DOS QUADRILÁTEROS
% =============================
quad = elem(isQuad, :);

if any(isQuad)
    A = coord(quad(:,1),:);
    B = coord(quad(:,2),:);
    C = coord(quad(:,3),:);
    D = coord(quad(:,4),:);

    % Triângulo 1: A-B-C
    v1 = B - A;
    v2 = B - C;
    A1 = 0.5 * vecnorm(cross(v1, v2, 2), 2, 2);

    % Triângulo 2: C-D-A
    v3 = D - C;
    v4 = D - A;
    A2 = 0.5 * vecnorm(cross(v3, v4, 2), 2, 2);

    elemarea(isQuad) = A1 + A2;
end

end
%--------------------------------------------------------------------------
%Function "getesure"
%--------------------------------------------------------------------------
%end
function [esureface1,esureface2,esurefull1,esurefull2] = getesure(elem,...
    bedgesize,inedge,esurn1,esurn2,nsurn2,rowposit)

% Inicializa ponteiros
esureface2 = 0;
esurefull2 = 0;

% Pré-alocação grosseira (pode ajustar depois se quiser otimizar memória)
nelem = size(elem,1);
maxdeg = 4; % elem tem até 4 nós
esureface1 = zeros(nelem*maxdeg,1);
esurefull1 = zeros(nelem*maxdeg*2,1); % full pode ter mais

mface = 0;
mfull = 0;

for ielem = 1:nelem
    neighborelemface = 0;
    neighborelemfull = 0;

    amountedge = sum(elem(ielem,1:4) ~= 0);

    m = 1;
    cface = 1;
    cfull = 1;
    boundflagface = 0;
    boundflagfull = 0;

    elemvertices = elem(ielem,1:amountedge);

    for j = 1:amountedge + 1
        vertex = elemvertices(1);
        next   = elemvertices(2);

        % linhas de rowposit associadas a cada vértice
        survtx1 = rowposit(nsurn2(vertex) + 1:nsurn2(vertex + 1));
        survtx2 = rowposit(nsurn2(next)   + 1:nsurn2(next   + 1));

        % interseção vectorizada (ignorando zeros)
        common = survtx1( survtx1 ~= 0 & ismember(survtx1,survtx2) );
        if isempty(common)
            pointrow = 0;
        else
            pointrow = common(1); % assume uma só aresta
        end

        if pointrow > bedgesize
            possiblelem = 1:2;
            pointactualelem = inedge(pointrow - bedgesize,3:4) ~= ielem;

            neigh = inedge(pointrow - bedgesize,2 + possiblelem(pointactualelem));

            neighborelemfull(cfull) = neigh;
            neighborelemface(cface) = neigh;

            if m == 2
                esurn = esurn1(esurn2(vertex) + 1:esurn2(vertex + 1));

                % reordena esurn
                newesurn = shiftchoosen(esurn,...
                    neighborelemfull(length(neighborelemfull) - 1),'val');

                elemvector = [ielem ...
                    neighborelemfull(length(neighborelemfull) - 1:...
                    length(neighborelemfull))];

                % setdiff vectorizado
                extraelem = setdiff(newesurn,elemvector,'stable');

                if ~isempty(extraelem)
                    neighboraux = zeros(1,length(extraelem)+2);
                    neighboraux(1) = neighborelemfull(length(neighborelemfull) - 1);
                    neighboraux(2:length(extraelem)+1) = extraelem;
                    neighboraux(length(extraelem)+2) = ...
                        neighborelemfull(length(neighborelemfull));

                    neighborelemfull(cfull - 1:cfull + length(neighboraux) - 2) = neighboraux;

                    m = 1;
                    cfull = cfull + length(extraelem);
                end
            end

            m = m + 1;
            cface = cface + 1;
            cfull = cfull + 1;

        else
            m = 1;
            boundflagface = cface;
            boundflagfull = cfull;
        end

        if amountedge == 3
            elemvertices = circshift(elemvertices,1:3);
        else
            elemvertices = circshift(elemvertices,2:4);
        end
    end

    neighborelemface = unique(neighborelemface,'stable');
    neighborelemfull = unique(neighborelemfull,'stable');

    if boundflagface ~= 0 && boundflagface <= length(neighborelemface)
        neighborelemface = shiftchoosen(neighborelemface,...
            neighborelemface(boundflagface),'val');
    end

    if boundflagfull ~= 0 && boundflagfull <= length(neighborelemfull)
        neighborelemfull = shiftchoosen(neighborelemfull,...
            neighborelemfull(boundflagfull),'val');
    end

    lenf = length(neighborelemface);
    lenfull = length(neighborelemfull);

    esureface1(mface + 1:mface + lenf) = neighborelemface;
    esurefull1(mfull + 1:mfull + lenfull) = neighborelemfull;

    esureface2 = vertcat(esureface2,mface + lenf);
    esurefull2 = vertcat(esurefull2,mfull + lenfull);

    mface = mface + lenf;
    mfull = mfull + lenfull;
end

% corta excesso de pré-alocação
esureface1 = esureface1(1:mface);
esurefull1 = esurefull1(1:mfull);
end
%--------------------------------------------------------------------------
%FUNCTION "getrowposition"
function rowposit = getrowposition(bedge, inedge, nsurn1, nsurn2)

rowposit = zeros(length(nsurn1),1);

% ============================================================
% PROCESSA BOUNDARY EDGES
% ============================================================

v1 = bedge(:,1);
v2 = bedge(:,2);
nb = size(bedge,1);

for k = 1:nb
    % vizinhos do nó v1(k)
    s = nsurn2(v1(k)) + 1;
    e = nsurn2(v1(k)+1);
    pos = find(nsurn1(s:e) == v2(k));
    rowposit(s + pos - 1) = k;

    % vizinhos do nó v2(k)
    s = nsurn2(v2(k)) + 1;
    e = nsurn2(v2(k)+1);
    pos = find(nsurn1(s:e) == v1(k));
    rowposit(s + pos - 1) = k;
end

% ============================================================
% PROCESSA INTERNAL EDGES
% ============================================================

v1 = inedge(:,1);
v2 = inedge(:,2);
ni = size(inedge,1);

for k = 1:ni
    idx = k + nb;  % deslocamento

    % vizinhos do nó v1(k)
    s = nsurn2(v1(k)) + 1;
    e = nsurn2(v1(k)+1);
    pos = find(nsurn1(s:e) == v2(k));
    rowposit(s + pos - 1) = idx;

    % vizinhos do nó v2(k)
    s = nsurn2(v2(k)) + 1;
    e = nsurn2(v2(k)+1);
    pos = find(nsurn1(s:e) == v1(k));
    rowposit(s + pos - 1) = idx;
end

end

function [node2edge, node2elem] = build_connectivities(bedge, inedge, elem, nNodes)

node2edge = cell(nNodes,1);
node2elem = cell(nNodes,1);

% --- Elementos na vizinhança do nó ---
for e = 1:size(elem,1)
    nd = elem(e, elem(e,:) > 0);
    for n = nd
        node2elem{n} = [node2elem{n}, e];
    end
end

% --- Nós conectados via aresta interna ---
for k = 1:size(inedge,1)
    i = inedge(k,1); j = inedge(k,2);
    node2edge{i} = [node2edge{i}, j];
    node2edge{j} = [node2edge{j}, i];
end

% --- Nós conectados via aresta de contorno ---
for k = 1:size(bedge,1)
    i = bedge(k,1); j = bedge(k,2);
    node2edge{i} = [node2edge{i}, j];
    node2edge{j} = [node2edge{j}, i];
end

node2edge = cellfun(@unique, node2edge, 'UniformOutput', false);

end

function [nsurn1] = reordernsurn(esurn1,esurn2,nsurn1,nsurn2,bedge,coord,...
    elem)
%Put "nsurn1" in a vector column
nsurn1 = nsurn1'; 

%Fernando's modification

% a aloca��o do vetor serve para identificar todos os nos
% do contorno e damos flag os demais que est�o no interior � zero
is_bound = zeros(size(coord,1),1);
for i = 1:size(bedge,1),
   is_bound(bedge(i,1)) = 1;
   % procura se tem algum no se repitidno no coluna 1 do bedge
   % quando passa isso que dizer que ordena��o do bedge esta 
   % sentido horario e/ou antihorario
   auxc = find(bedge(:,1) == bedge(i,1));
   if ~isempty(auxc) && length(auxc) > 1
       is_bound(bedge(auxc(1),1)) = 1;
   end
   % procura se o n� esta se repitindo na coluna 2 bedge
  
    auxc1 = find(bedge(:,2) == bedge(i,2));
    if ~isempty(auxc1) && length(auxc1) > 1
        is_bound(bedge(auxc1(1),2))=1;  
    end
end  %End of FOR

% is_bound = zeros(size(coord,1),1);
% 
% for i=1:size(bedge,1),
%    is_bound(bedge(i,1))=1; 
% end

for i=1:size(coord,1),
    if is_bound(i)==0
        nnsn=zeros(nsurn2(i+1)-nsurn2(i),1);
        prim=esurn1(esurn2(i)+1);
        ult=esurn1(esurn2(i+1));
        if elem(prim,4)==0
            bp=3;
        else
            bp=4;
        end
        if elem(ult,4)==0
            bu=3;
        else
            bu=4;
        end
        for j=1:bp,
            for u=1:bu,
                if ((elem(ult,u)==elem(prim,j))&&(elem(prim,j)~=i))
                    nnsn(1)=elem(prim,j);
                end
            end
        end
        for t=2:(esurn2(i+1)-esurn2(i)),
            atual=esurn1(esurn2(i)+t);
            ant=esurn1(esurn2(i)+t-1);
            if elem(atual,4)==0
                bp=3;
            else
                bp=4;
            end
            if elem(ant,4)==0
                bu=3;
            else
                bu=4;
            end
            for j=1:bp,
                for u=1:bu,
                    if ((elem(ant,u)==elem(atual,j))&&(elem(atual,j)~=i))
                        nnsn(t)=elem(atual,j);
                    end
                end
            end
        end
        for t=1:(nsurn2(i+1)-nsurn2(i)),
            nsurn1(nsurn2(i)+t)=nnsn(t);
        end     
    else 
        nnsn=zeros(nsurn2(i+1)-nsurn2(i),1);
        prim=esurn1(esurn2(i)+1);
        ult=esurn1(esurn2(i+1));
        if elem(prim,4)==0
            bp=3;
        else
            bp=4;
        end
        if elem(ult,4)==0
            bu=3;
        else
            bu=4;
        end
        for j=1:size(bedge,1)
            for p=1:bp,
                if (i==elem(prim,p))
                    if p==bp
                        n1=elem(prim,1);
                    else
                        n1=elem(prim,p+1);
                    end
                end
            end
            nnsn(1)=n1;            
            for u=1:bu,
                if (i==elem(ult,u))
                    if u==1
                        n1=elem(ult,bu);
                    else
                        n1=elem(ult,u-1);
                    end
                end
            end
            nnsn((nsurn2(i+1)-nsurn2(i)))=n1;
        end
        for t=2:(esurn2(i+1)-esurn2(i)),
            atual=esurn1(esurn2(i)+t);
            ant=esurn1(esurn2(i)+t-1);
            if elem(atual,4)==0
                bp=3;
            else
                bp=4;
            end
            if elem(ant,4)==0
                bu=3;
            else
                bu=4;
            end
            for j=1:bp,
                for u=1:bu,
                    if ((elem(ant,u)==elem(atual,j))&&(elem(atual,j)~=i))
                        nnsn(t)=elem(atual,j);
                    end
                end
            end
        end
        for t=1:(nsurn2(i+1)-nsurn2(i)),
            nsurn1(nsurn2(i)+t)=nnsn(t);
        end
    end            
end




end