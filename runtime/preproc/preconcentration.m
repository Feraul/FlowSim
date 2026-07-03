
%Modified: Fernando Contreras, 2021


function [env,paramconcentra]= ...
    preconcentration(env,wells)
global numcase centelem


auxpar=0;
wellsc=zeros(2,6);
%% media velocities
if numcase==233
    velmedio=0.5; % quando a pessão no contorno é 5
elseif numcase==235
    velmedio=2;    % quando a pressão no contorno é 15
elseif numcase==231 || numcase==232 || numcase==243 || ...
        numcase==236 || numcase==237 || numcase==238 ||...
        numcase==239 || numcase==241 || numcase==242 ||...
        numcase==244  || numcase==234|| numcase==233 ||...
        numcase==248 || numcase==251 || numcase==380 || numcase==380.1
    velmedio=1;
end
%% calculate diffusion tensor
[dmap,Dmedio,gamma] = PLUG_dfunction;

%% flags adicionais para o problema 5.3.1-5.3.2
switch numcase
    case 231
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % flag de concentracao tipo Dirichlet
        bcflagc(2,2)=10; % valor concentracao
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor do fluxo de concentracao
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % bcflag(6,1)=50;  % dicichlet concen
        % bcflag(7,1)=51; % dirich concen
        % bcflag(8,1)=250;  % dicichlet concen
        % bcflag(9,1)=251; % dirich concen
        % bcflag(6,2)=10; % valor dirich
        % bcflag(7,2)=10;  % valor concent
        % bcflag(8,2)=0; % valor dirich
        % bcflag(9,2)=0;  % valor concent
        %
        % c=logical(bedge(:,4)==201);
        % bedge([c],6)=201;
        % d=logical(bedge(:,4)==101);
        % bedge([d],6)=50;
        % e=logical(bedge(:,4)==102);
        % bedge([e],6)=51;
        % f=logical(bedge(:,4)==103);
        % bedge([f],6)=250;
        % g=logical(bedge(:,4)==104);
        % bedge([g],6)=251;
        %
        % c1=logical(bedge(:,5)==201);
        % bedge([c1],7)=201;
        % d1=logical(bedge(:,5)==101);
        % bedge([d1],7)=50;
        % e1=logical(bedge(:,5)==102);
        % bedge([e1],7)=51;
        % f1=logical(bedge(:,5)==103);
        % bedge([f1],7)=250;
        % g1=logical(bedge(:,5)==104);
        % bedge([g1],7)=251;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 232
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
    case 233
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
        
    case 234
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
    case 235
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
    case 236
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
    case 237
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
    case 238
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
    case 239
        bcflagc(1,1)=50;% dicichlet concen
        bcflagc(1,2)=0;% valor dirich
        bcflagc(2,1)=51; % dicichlet concen
        bcflagc(2,2)=1; % valor dirich
        bcflagc(3,1)=250; % neumann boundary
        bcflagc(3,2)=0;  % valor
        bcflagc(4,1)=251;  % neumann boundary
        bcflagc(4,2)=0; % valor
        
        c=logical(bedge(:,4)==101);
        bedge([c],6)=50;
        d=logical(bedge(:,4)==102);
        bedge([d],6)=51;
        e=logical(bedge(:,4)==103);
        bedge([e],6)=250;
        f=logical(bedge(:,4)==201);
        bedge([f],6)=251;
        
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=50;
        d1=logical(bedge(:,5)==102);
        bedge([d1],7)=51;
        e1=logical(bedge(:,5)==103);
        bedge([e1],7)=250;
        f1=logical(bedge(:,5)==201);
        bedge([f1],7)=251;
        wells=0;
    case 241
        bcflagc(1,1)=50;% dicichlet concen
        bcflagc(1,2)=0;% valor dirich
        bcflagc(2,1)=51; % dicichlet concen
        bcflagc(2,2)=1; % valor dirich
        bcflagc(3,1)=250; % neumann boundary
        bcflagc(3,2)=0;  % valor
        bcflagc(4,1)=251;  % neumann boundary
        bcflagc(4,2)=0; % valor
        
        c=logical(bedge(:,4)==101);
        bedge([c],6)=50;
        d=logical(bedge(:,4)==102);
        bedge([d],6)=51;
        e=logical(bedge(:,4)==103);
        bedge([e],6)=250;
        f=logical(bedge(:,4)==201);
        bedge([f],6)=251;
        
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=50;
        d1=logical(bedge(:,5)==102);
        bedge([d1],7)=51;
        e1=logical(bedge(:,5)==103);
        bedge([e1],7)=250;
        f1=logical(bedge(:,5)==201);
        bedge([f1],7)=251;
        wells=0;
    case 242
        bcflagc(1,1)=50;% dicichlet concen
        bcflagc(1,2)=0;% valor dirich
        bcflagc(2,1)=51; % dicichlet concen
        bcflagc(2,2)=1; % valor dirich
        bcflagc(3,1)=250; % neumann boundary
        bcflagc(3,2)=0;  % valor
        bcflagc(4,1)=251;  % neumann boundary
        bcflagc(4,2)=0; % valor
        
        c=logical(bedge(:,4)==101);
        bedge([c],6)=50;
        d=logical(bedge(:,4)==102);
        bedge([d],6)=51;
        e=logical(bedge(:,4)==103);
        bedge([e],6)=250;
        f=logical(bedge(:,4)==201);
        bedge([f],6)=251;
        
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=50;
        d1=logical(bedge(:,5)==102);
        bedge([d1],7)=51;
        e1=logical(bedge(:,5)==103);
        bedge([e1],7)=250;
        f1=logical(bedge(:,5)==201);
        bedge([f1],7)=251;
        wells=0;
    case 243
        bcflagc(1,1)=50;% dicichlet concen
        bcflagc(1,2)=0;% valor dirich
        bcflagc(2,1)=51; % dicichlet concen
        bcflagc(2,2)=1; % valor dirich
        bcflagc(3,1)=250; % neumann boundary
        bcflagc(3,2)=0;  % valor
        bcflagc(4,1)=251;  % neumann boundary
        bcflagc(4,2)=0; % valor
        
        c=logical(bedge(:,4)==101);
        bedge([c],6)=50;
        d=logical(bedge(:,4)==102);
        bedge([d],6)=51;
        e=logical(bedge(:,4)==103);
        bedge([e],6)=250;
        f=logical(bedge(:,4)==201);
        bedge([f],6)=251;
        
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=50;
        d1=logical(bedge(:,5)==102);
        bedge([d1],7)=51;
        e1=logical(bedge(:,5)==103);
        bedge([e1],7)=250;
        f1=logical(bedge(:,5)==201);
        bedge([f1],7)=251;
        %         bcflagc(1,1)=251;
        %         bcflagc(1,2)=0;
        %         bcflagc(2,1)=50;  % dicichlet concen
        %         bcflagc(2,2)=10; % valor dirich
        %         bcflagc(3,1)=250; % dirich neumann
        %         bcflagc(3,2)=0;  % valor concent
        %
        %         c=logical(bedge(:,4)==201);
        %         bedge([c],6)=251;
        %         d=logical(bedge(:,4)==101);
        %         bedge([d],6)=50;
        %         e=logical(bedge(:,4)==102);
        %         bedge([e],6)=250;
        %
        %         c1=logical(bedge(:,5)==201);
        %         bedge([c1],7)=251;
        %         d1=logical(bedge(:,5)==101);
        %         bedge([d1],7)=50;
        %         e1=logical(bedge(:,5)==102);
        %         bedge([e1],7)=250;
        wells=0;
    case 244
        bcflagc(1,1)=250;
        bcflagc(1,2)=0;
        
        c=logical(bedge(:,4)==101);
        bedge([c],6)=250;
        
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=250;
        wells (1,1:6)=[1201 1 301 1 0 1];
    case 245
        bcflagc(1,1)=250;
        bcflagc(1,2)=0;
        c=logical(bedge(:,4)==201);
        bedge([c],6)=250;
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=250;
    case 247
        bcflagc(1,1)=250;
        bcflagc(1,2)=0;
        c=logical(bedge(:,4)==201);
        bedge([c],6)=250;
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=250;
        auxpar=1;
        %         wellsc(1,1)=1620;
        %         wellsc(1,2)=1;
        %         wellsc(1,3)=301;
        %         wellsc(1,4)=1;
        %         wellsc(1,5)=401;
        %         wellsc(1,6)=14.2146;
        wellsc(2,:)=wells(2,:);
        wellsc(1,:)=wells(1,:);
        
    case 248
        bcflagc(1,1)=50;
        bcflagc(1,2)=0;
        c=logical(bedge(:,4)==101);
        bedge([c],6)=50;
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=50;
    case 249
        
        
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=1; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        j=1;
        for i=1:size(bedge,1)
            if bedge(i,7)==250
                wells(j,1)=bedge(i,3);
                wells(j,2)=2;
                wells(j,3)=0;
                wells(j,4)=0;
                wells(j,5)=501;
                wells(j,6)=0;
                j=j+1;
            end
        end
        auxpar=1;
    case 250
        
        bcflagc(1,1)=250;
        bcflagc(1,2)=0;
        c=logical(bedge(:,4)==201);
        bedge([c],6)=250;
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=250;
        
    case 246
        bcflagc(1,1)=50;% dicichlet concen
        bcflagc(1,2)=0;% valor dirich
        bcflagc(2,1)=51; % dicichlet concen
        bcflagc(2,2)=1; % valor dirich
        bcflagc(3,1)=250; % neumann boundary
        bcflagc(3,2)=0;  % valor
        bcflagc(4,1)=251;  % neumann boundary
        bcflagc(4,2)=0; % valor
        
        c=logical(bedge(:,4)==101);
        bedge([c],6)=50;
        d=logical(bedge(:,4)==102);
        bedge([d],6)=51;
        e=logical(bedge(:,4)==103);
        bedge([e],6)=250;
        f=logical(bedge(:,4)==201);
        bedge([f],6)=251;
        
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=50;
        d1=logical(bedge(:,5)==102);
        bedge([d1],7)=51;
        e1=logical(bedge(:,5)==103);
        bedge([e1],7)=250;
        f1=logical(bedge(:,5)==201);
        bedge([f1],7)=251;
    case 251
        bcflagc(1,1)=50;% dicichlet concen
        bcflagc(1,2)=1;% valor dirich
        bcflagc(2,1)=251; % Neumann
        bcflagc(2,2)=0; % valor dirich
        bcflagc(3,1)=250; % neumann boundary
        bcflagc(3,2)=0;  % valor
        
        
        c=logical(bedge(:,4)==101);
        bedge([c],6)=50;
        d=logical(bedge(:,4)==102);
        bedge([d],6)=251;
        f=logical(bedge(:,4)==201);
        bedge([f],6)=250;
        
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=50;
        d1=logical(bedge(:,5)==102);
        bedge([d1],7)=251;
        f1=logical(bedge(:,5)==201);
        bedge([f1],7)=250;
    case 380
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        %bcflagc(2,1)=251;  % dicichlet concen
        %bcflagc(2,2)=0; % valor dirich
        %bcflagc(3,1)=251; % dirich neumann
        %bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=251;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=251;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=251;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=251;
        
        b=find((10.85<centelem(:,1) & centelem(:,1)<10.95) & (0<centelem(:,2) & centelem(:,2)<1));
        b1=find((10.85<centelem(:,1) & centelem(:,1)<10.95) & (4<centelem(:,2) & centelem(:,2)<5));
        %bb=logical;
        wellsc(1:10,1)=b;
        wellsc(11:20,1)=b1;
        wellsc(:,2)=2;
        wellsc(:,3)=380;
        wellsc(:,4)=100;
        wellsc(:,5)=0;
        wellsc(:,6)=0;
    case 380.1
        bcflagc(1,1)=50;% dicichlet concen
        bcflagc(1,2)=0;% valor dirich
        bcflagc(2,1)=51; % dicichlet concen
        bcflagc(2,2)=0; % valor dirich
        bcflagc(3,1)=52; % dirichlet boundary
        bcflagc(3,2)=100;  % valor
        bcflagc(4,1)=251;  % neumann boundary
        bcflagc(4,2)=0; % valor
        
        c=logical(bedge(:,4)==101);
        bedge([c],6)=50;
        d=logical(bedge(:,4)==102);
        bedge([d],6)=51;
        e=logical(bedge(:,4)==103);
        bedge([e],6)=52;
        f=logical(bedge(:,4)==201);
        bedge([f],6)=251;
        
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=50;
        d1=logical(bedge(:,5)==102);
        bedge([d1],7)=51;
        e1=logical(bedge(:,5)==103);
        bedge([e1],7)=52;
        f1=logical(bedge(:,5)==201);
        bedge([f1],7)=251;


end

env.geometry.bedge=bedge;
% end flags
%% flags adicionais para o problema 5.3.3 Darlan
% if strcmp(Prob,'2D')
%     bcflag(8,1)=50; % flag aux para a concentracáo no médio
%     bcflag(9,1)=51; % flag aux para a concentracáo acima e embaixo
%     bcflag(10,1)=52; % flag aux para a concentracáo acima e embaixo
%     bcflag(11,1)=250; % flag aux para a concentracáo face 2
%     bcflag(12,1)=251; % flag aux para a concentracáo face 2
%     bcflag(13,1)=252; % flag aux para a concentracáo face 2
%     bcflag(8,2)=0; % valor da Concentração em bcflag(9,1)
%     bcflag(9,2)=0; % valor da Concentração em bcflag(10,1)
%
%     % if t_old < 0.2
%     bcflag(10,2)=1; % valor da Concentração em bcflag(11,1)
%     %  else
%     %      bcflag(10,2)=0; % valor da Concentração em bcflag(11,1)
%     %  end
%     bcflag(11,2)=0; % valor da Concentração em bcflag(11,1)
%     bcflag(12,2)=0; % valor da Concentração em bcflag(11,1)
%     bcflag(13,2)=0; % valor da Concentração em bcflag(11,1)
%
%     c=logical(bedge(:,4)==101);
%     bedge([c],6)=50;
%     d=logical(bedge(:,4)==102);
%     bedge([d],6)=51;
%     e=logical(bedge(:,4)==103);
%     bedge([e],6)=52;
%     f=logical(bedge(:,4)==104);
%     bedge([f],6)=250;
%     g=logical(bedge(:,4)==105);
%     bedge([g],6)=251;
%     h=logical(bedge(:,4)==106);
%     bedge([h],6)=252;
%     i=logical(bedge(:,4)==201);
%     bedge([i],6)=201;
%
%     c1=logical(bedge(:,5)==101);
%     bedge([c1],7)=50;
%     d1=logical(bedge(:,5)==102);
%     bedge([d1],7)=51;
%     e1=logical(bedge(:,5)==103);
%     bedge([e1],7)=52;
%     f1=logical(bedge(:,5)==104);
%     bedge([f1],7)=250;
%     g1=logical(bedge(:,5)==105);
%     bedge([g1],7)=251;
%     h1=logical(bedge(:,5)==106);
%     bedge([h1],7)=252;
%     i1=logical(bedge(:,5)==201);
%     bedge([i1],7)=201;
% end
% % end flags
% %%
%
% nflagp=5000*ones(size(coord,1),2); % flags 7.15934 é uma constante
% nflagc=5000*ones(size(coord,1),2); % flags 7.15934 é uma constante
%
% for ifacont=1:size(bedge,1)
%
%     x=bcflag(:,1)==bedge(ifacont,4);
%     r=find(x==1);
%     nflagp(bedge(ifacont,1),2)=bcflag(r,2);
%     nflagp(bedge(ifacont,1),1)=bcflag(r,1);
%
%     x=bcflag(:,1)==bedge(ifacont,6);
%     r=find(x==1);
%     nflagc(bedge(ifacont,1),2)=bcflag(r,2);
%     nflagc(bedge(ifacont,1),1)=bcflag(r,1);
% end
paramconcentra.bcflagc=bcflag;
paramconcentra.wellsc=wellsc;
paramconcentra.auxpar=auxpar;
paramconcentra.velmedio=velmedio;
paramconcentra.dmap=dmap;
paramconcentra.Dmedio=Dmedio;
paramconcentra.gamma=gamma;
end