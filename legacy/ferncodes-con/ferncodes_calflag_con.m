%It is called by "preMPFA.m"

function [nflagnoc,nflagfacec] = ferncodes_calflag_con(a)
% determinar o flag do nó interior e fronteira de Neumann
global coord bedge bcflagc inedge numcase;
nflagnoc = 5000*ones(size(coord,1),3);
nflagnoc(:,2:3)=0;

nflagfacec=5000*ones(size(bedge,1)+size(inedge,1),3);
nflagfacec(:,2:3)=0;

for ifacont = 1:size(bedge,1)

    %----------------------------------------------------------------------
    % flag no vertice
    x = logical(bcflagc(:,1) == bedge(ifacont,6));
    vertex = bedge(ifacont,1);
    %First column receives the boundary condition flag.

     nflagnoc(vertex,1) = bcflagc(x,1);

    if numcase==248
        % avalia o valor da concentração no contorno utilizando a solução
        % exata

        x1 = coord(vertex,1);
        y1 = coord(vertex,2);
        x= x1+y1;
        %Second column receives the boundary condition value.
        nflagnoc(vertex,2) = PLUG_bcfunction_con(x,a);
    else
        %Second column receives the boundary condition value.

        nflagnoc(vertex,2) = PLUG_bcfunction_con(x,a);

    end
    % primeira coluna flag das condicoes de contorno de Dirichlet
    % segunda coluna valores das concentrações impostos sobre as condicoes
    % de contorno
    % terceira coluna, é flag de no vertice.

    if nflagnoc(vertex,1)<100

        nflagnoc(vertex,3)=1;
    else
        nflagnoc(vertex,3)=0;
    end

    %% --------------------------------------------------------------------
    % flag sobre a face

    y = logical(bcflagc(:,1) == bedge(ifacont,7));
    %First column receives the boundary condition flag.
    nflagfacec(ifacont,1) = bcflagc(y,1);
    if numcase==248
        % avalia o valor da concentração no contorno utilizando a solução
        vertex2 = bedge(ifacont,2);
        auxvertex= 0.5*(coord(vertex,:)+coord(vertex2,:));
        x1 = auxvertex(1,1);
        y1 = auxvertex(1,2);
        y= x1+y1;
        %Second column receives the boundary condition value.
        nflagfacec(ifacont,2) = PLUG_bcfunction_con(y,a);
    else
        %Second column receives the boundary condition value.
        nflagfacec(ifacont,2) = PLUG_bcfunction_con(y,a);
    end
    if nflagfacec(ifacont,1)<100
        nflagfacec(ifacont,3)=1;
    else
        nflagfacec(ifacont,3)=0;
    end

end  %End of FOR