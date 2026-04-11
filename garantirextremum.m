function [Sleft,Sright]=garantirextremum(Sleft,Sright,Sw,leftelem,...
    rightelem,dotvn,i,vertices,verticescoord,taylorterms,...
    limiterflag,constraint,flagknownvert,satonvertices,countinter,mlplimiter)
global centelem

% verificação dos extremos
%==========================================================================
if dotvn<0 && abs(dotvn)>1e-5
    if Sw(leftelem)<Sleft || Sleft<Sright
        % Sleft=Sw(leftelem);
        %==================================================================
        %Define the elements that share the edge evaluated
        elemeval = [leftelem rightelem];
        %Get the saturation value recovered on each quadrature point ("on_q")
        %Get the statment regarding to mlp limiter:
        boolmlp = (length(mlplimiter) == 1);
        mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
        limiterflagaux=limiterflag;
        limiterflagaux{1}='on';
        limiterflagaux{4}='wf';
        limiterflagaux{7}=-1;
        faceorderauxlef=2;
        Sleft = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,...
            limiterflagaux,faceorderauxlef,constraint,flagknownvert,satonvertices,...
            mlpbyelem,centelem(elemeval,1:2),countinter);
        if Sw(leftelem)<Sleft || Sleft<Sright
            Sleft=Sw(leftelem);
        end
        %==================================================================
    end
    if Sright>Sleft || Sright<Sw(rightelem)
        %            Sright=Sw(rightelem);
        %==================================================================
        % Define the elements that share the edge evaluated
        elemeval = [rightelem leftelem];
        % Get the saturation value recovered on each quadrature point ("on_q")
        %     Get the statment regarding to mlp limiter:
        boolmlp = (length(mlplimiter) == 1);
        %     Get the statment regarding to mlp limiter:
        mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
        limiterflagaux=limiterflag;
        limiterflagaux{1}='on';
        limiterflagaux{4}='wf';
        limiterflagaux{7}=-1;
        faceorderauxright=2;
        Sright = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,...
            limiterflagaux,faceorderauxright,constraint,flagknownvert,satonvertices,...
            mlpbyelem,centelem(elemeval,1:2),countinter);
        if Sright>Sleft || Sright<Sw(rightelem)
            Sright=Sw(rightelem);
        end
        %==================================================================
    end
%     if Sw(leftelem)<=Sleft && Sleft<=Sright && Sright<=Sw(rightelem)
%     else
%         i
%         disp('viola extremo aqui1')
%         pause
%     end
elseif dotvn>0 &&  abs(dotvn)>1e-5
    if Sw(leftelem)<Sleft || Sleft<Sright
        %             Sleft=Sw(leftelem);
        %==================================================================
        %Define the elements that share the edge evaluated
        elemeval = [leftelem rightelem];
        %Get the saturation value recovered on each quadrature point ("on_q")
        %Get the statment regarding to mlp limiter:
        boolmlp = (length(mlplimiter) == 1);
        mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
        limiterflagaux=limiterflag;
        limiterflagaux{1}='on';
        limiterflagaux{4}='wf';
        limiterflagaux{7}=-1;
        faceorderauxlef=2;
        
        Sleft = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,...
            limiterflagaux,faceorderauxlef,constraint,flagknownvert,satonvertices,...
            mlpbyelem,centelem(elemeval,1:2),countinter);
        if Sw(leftelem)<Sleft || Sleft<Sright
            Sleft=Sw(leftelem);
        end
        %==================================================================
    end
    if Sleft< Sright || Sright<Sw(rightelem)
        %Sright=Sw(rightelem);
        %==================================================================
        % Define the elements that share the edge evaluated
        elemeval = [rightelem leftelem];
        % Get the saturation value recovered on each quadrature point ("on_q")
        %     Get the statment regarding to mlp limiter:
        boolmlp = (length(mlplimiter) == 1);
        %     Get the statment regarding to mlp limiter:
        mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
        limiterflagaux=limiterflag;
        limiterflagaux{1}='on';
        limiterflagaux{4}='wf';
        limiterflagaux{7}=-1;
        faceorderauxright=2;
        Sright = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,...
            limiterflagaux,faceorderauxright,constraint,flagknownvert,satonvertices,...
            mlpbyelem,centelem(elemeval,1:2),countinter);
        if Sleft< Sright || Sright<Sw(rightelem)
            Sright=Sw(rightelem);
        end
        %==================================================================
    end
%     if Sw(leftelem)>=Sleft && Sleft>=Sright && Sright>=Sw(rightelem)
%     else
%         i
%         disp('viola extremo aqui2')
%         pause
%     end
elseif abs(dotvn)<1e-5
    Sleft=Sleft;%Sw(leftelem);
    Sright=Sright;% Sw(rightelem);
end


end