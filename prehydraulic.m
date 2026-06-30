
% Esta funcao determina os parametros iniciais como adequacao dos pocos
% artesanais, condicao inicial hydrologica e alguns outros
% parametro fisicos necessarios para rodar os problemas, alem disso, pode
% colocar os parametros para os novos casos
function [parmgroundwater,source_wells]=prehydraulic(env,source_wells)

numcase=env.config.numcase;
elem=env.geometry.elem;
centelem=env.geometry.centelem;
% inicialize os parametros
SS=0; h_init=0; MM=0; wells=0; dt=0; P=0;

switch numcase
    % The cases 330-333 were obtained of the article -- A local grid-refined
    % numerical groundwater model based on the vertex centered finite
    % volume method
    case 330
        % Case 1: single pumping well in a confined aquifer (Y. Qian et al)
        % initially hydraulic charge
        h_init=100*ones(size(elem,1),1);
        % the espeficied storage
        SS=0.001;
        % aquifer thickness
        MM=3;
        % step time
        dt=0.01;
        % find the well element
        %b=find(abs(centelem(:,1)-500)<1e-9 & abs(centelem(:,2)-500)<1e-9);
        b=find((abs(centelem(:,1)-500)/500)<1e-2 & (abs(centelem(:,2)-500)/500)<1e-2);
        wells(1,1)=b;
        wells(1,2)=2;
        wells(1,3)=0;
        wells(1,4)=0;
        wells(1,5)=0;
        wells(1,6)=-10000;
    case 331
        % Case 2: single pumping well in a unconfined aquifer (Y. Qian et al)

        % initially hydraulic charge
        h_init=90*ones(size(elem,1),1);
        % the specified yield
        SS=0.1;
        % aquifer thickness
        MM=100;
        % step time
        dt=0.01;
        % find the well element
        b=find((abs(centelem(:,1)-500)/500)<1e-2 & (abs(centelem(:,2)-500)/500)<1e-2);
        wells(1,1)=b;
        wells(1,2)=2;
        wells(1,3)=0;
        wells(1,4)=0;
        wells(1,5)=0;
        wells(1,6)=-40000;
        h_old=1*ones(size(elem,1),1);
    case 332
        % Case 3: multiple pumping wells in a confined aquifer (Y. Qian et al)

        % initially hydraulic charge
        h_init=100*ones(size(elem,1),1);
        % coeficiente de armazenamento especifico
        SS=0.001;
        % espesura do aquifero
        MM=3;
        % step time
        dt=0.01;
        % number of divisions
        const=25;
        pumpingrate=2500;
        % find the well element

        % the flow rate value in the well is divided by number wells

        %b1=find(abs(centelem(:,1)-250)<1e-9 & abs(centelem(:,2)-250)<1e-9);
        b1=find((const*10<centelem(:,1)& centelem(:,1)<const*11) & (const*10<centelem(:,2)& centelem(:,2)<const*11));
        %b1=find((abs(centelem(:,1)-250)/250)<1e-2 & (abs(centelem(:,2)-250)/250)<1e-2);
        wells(1,1)=b1;
        wells(1,2)=2;
        wells(1,3)=0;
        wells(1,4)=0;
        wells(1,5)=0;
        wells(1,6)=-pumpingrate;
        %-----------------------------------------------------------------
        b2=find((const*29<centelem(:,1)& centelem(:,1)<const*30) & (const*10<centelem(:,2)& centelem(:,2)<const*11));
        %b2=find((abs(centelem(:,1)-750)/750)<1e-2 & (abs(centelem(:,2)-250)/250)<1e-2);
        %b2=find(abs(centelem(:,1)-750)<1e-9 & abs(centelem(:,2)-250)<1e-9);
        wells(2,1)=b2;
        wells(2,2)=2;
        wells(2,3)=0;
        wells(2,4)=0;
        wells(2,5)=0;
        wells(2,6)=-pumpingrate;
        %------------------------------------------------------------------
        b3=find((const*29<centelem(:,1)& centelem(:,1)<const*30) & (const*29<centelem(:,2)& centelem(:,2)<const*30));
        %b3=find((abs(centelem(:,1)-750)/750)<1e-2 & (abs(centelem(:,2)-750)/750)<1e-2);

        wells(3,1)=b3;
        wells(3,2)=2;
        wells(3,3)=0;
        wells(3,4)=0;
        wells(3,5)=0;
        wells(3,6)=-pumpingrate;
        %------------------------------------------------------------------
        b4=find((const*10<centelem(:,1)& centelem(:,1)<const*11) & (const*29<centelem(:,2)& centelem(:,2)<const*30));
        %b4=find((abs(centelem(:,1)-250)/250)<1e-2 & (abs(centelem(:,2)-750)/750)<1e-2);

        wells(4,1)=b4;
        wells(4,2)=2;
        wells(4,3)=0;
        wells(4,4)=0;
        wells(4,5)=0;
        wells(4,6)=-pumpingrate;
    case 333
        % Case 4:Two parallel canals (Qian et al)
        elem(:,5)=1:size(elem,1);
        % initially hydraulic charge
        h_init=2*ones(size(elem,1),1);
        % coeficiente de armazenamento especifico
        SS=0.1;
        % espesura do aquifero
        MM=5;
        % step time no considered
        dt=20;
        % precipitation infiltration
        P=0.002;
    case 334
        % Case flow between two rivers (Mark Bakker), pag. 10
        % aquifer confined flow
        % aquifer thickness
        MM=10;
    case 335
        % Case areal recharge between two rivers (Mark Bakker), pag. 13
        % aquifer unconfined flow
        % aquifer thickness
        MM=10;
        % precipitation infiltration
        P=0.001;
    case 337
        % Case areal recharge between an impermeable noundary and a river
        % (Mark Bakker), pag. 19
        % aquifer unconfined flow
        % aquifer thickness
        MM=10;
        % precipitation infiltration
        P=0.001;
    case 338
        % Case areal recharge between an impermeable noundary and a river
        % with river bed resistance
        % (Mark Bakker), pag. 21
        % aquifer unconfined flow
        % aquifer thickness
        MM=10;
        % precipitation infiltration
        P=0.001;
    case 339
        % Case flow through two zones of different transmissiviteis
        % (Mark Bakker), pag. 23
        % aquifer unconfined flow
        % aquifer thickness
        MM=10;
    case 342
        % book: introduction to groundwater: Herbert Wang, page 80

        % initially hydraulic charge
        h_init=0*ones(size(elem,1),1);
        % the espeficied storage
        SS=3.28*10^-3;
        % aquifer thickness
        MM= 3;
        % step time, change for each time
        dt=1;
        % find the well element
        %wells(1,1)=1;
        %wells(1,2)=2;
        %wells(1,3)=0;
        %wells(1,4)=0;
        %wells(1,5)=0;
        %wells(1,6)=-2000; % pumping rate
    case 343
        h_init=0*ones(size(elem,1),1);
        % the espeficied storage
        SS=0.01;
        % aquifer thickness
        MM= 1;
        % step time, change for each time
        dt=10; % 1, 2, 5, 10, 20 and 50
    case 347
        % initially hydraulic charge
        h_init=1034.5*ones(size(elem,1),1);
        % coeficiente de armazenamento especifico
        SS=0.04;
        % espesura do aquifero
        MM=136;
        % step time no considered
        dt=1;
        % precipitation or infiltration
        P=0.042; %(4.2mm)
        %pumping rate
        pumpingrate=-261; %261


        b1 = 1; % Índice para wells
        %wells(b1, 1) = 125;            % Define o índice do poço
        wells(b1, 1) = 38;   
        wells(b1, 2) = 2;    % Referencia ao elemento correspondente
        wells(b1, 3) = 0;             % Valor placeholder (exemplo: coordenada x)
        wells(b1, 4) = 0;             % Valor placeholder (exemplo: coordenada y)
        wells(b1, 5) = 0;             % Valor placeholder (exemplo: coordenada z)
        wells(b1, 6) = pumpingrate;   % Define a taxa de bombeamento

        b2 = 2; % Índice para wells
        %wells(b2, 1) = 230;
        wells(b2, 1) = 5;
        wells(b2, 2) = 2;
        wells(b2, 3) = 0;
        wells(b2, 4) = 0;
        wells(b2, 5) = 501;
        wells(b2, 6) = 88;

        b3 = 3; % Índice para wells
        %wells(b3, 1) = 238;
        wells(b3, 1)  =15;
        wells(b3, 2) = 2;
        wells(b3, 3) = 0;
        wells(b3, 4) = 0;
        wells(b3, 5) = 0;
        wells(b3, 6) = pumpingrate;

        b4 = 4; % Índice para wells
        %wells(b4, 1) = 371;
        wells(b4, 1) = 45;
        wells(b4, 2) = 2;
        wells(b4, 3) = 0;
        wells(b4, 4) = 0;
        wells(b4, 5) = 501;
        wells(b4, 6) = 88; % Neste caso, a taxa de bombeamento é negativa

        b5 = 5; % Índice para wells
        %wells(b5, 1) = 194;
        wells(b5, 1) = 50;
        wells(b5, 2) = 2;
        wells(b5, 3) = 0;
        wells(b5, 4) = 0;
        wells(b5, 5) = 0;
        wells(b5, 6) = pumpingrate;
% 
%         b6 = 6; % Índice para wells GAMELEIRA
%         wells(b6, 1) = 347;
%         wells(b6, 2) = 2;
%         wells(b6, 3) = 0;
%         wells(b6, 4) = 0;
%         wells(b6, 5) = 0;
%         wells(b6, 6) = pumpingrate;
% 
%         b7 = 7; % Índice para wells LAGOINHA
%         wells(b7, 1) = 377;
%         wells(b7, 2) = 2;
%         wells(b7, 3) = 0;
%         wells(b7, 4) = 0;
%         wells(b7, 5) = 0;
%         wells(b7, 6) = pumpingrate;
% 
%         b8 = 8; % Índice para wells PAROIS
%         wells(b8, 1) = 396;
%         wells(b8, 2) = 2;
%         wells(b8, 3) = 0;
%         wells(b8, 4) = 0;
%         wells(b8, 5) = 0;
%         wells(b8, 6) = pumpingrate;
% 
%         b9 = 9;
%         wells(b9, 1) = 147;
%         wells(b9, 2) = 2;
%         wells(b9, 3) = 0;
%         wells(b9, 4) = 0;
%         wells(b9, 5) = 0;
%         wells(b9, 6) = pumpingrate;
% 
%         b10 = 10;
%         wells(b10, 1) = 401;
%         wells(b10, 2) = 2;
%         wells(b10, 3) = 0;
%         wells(b10, 4) = 0;
%         wells(b10, 5) = 0;
%         wells(b10, 6) = pumpingrate;
% 
%         b11 = 11;
%         wells(b11, 1) = 196;
%         wells(b11, 2) = 2;
%         wells(b11, 3) = 0;
%         wells(b11, 4) = 0;
%         wells(b11, 5) = 0;
%         wells(b11, 6) = pumpingrate;
% 
%         b12 = 12;
%         wells(b12, 1) = 434;
%         wells(b12, 2) = 2;
%         wells(b12, 3) = 0;
%         wells(b12, 4) = 0;
%         wells(b12, 5) = 0;
%         wells(b12, 6) = pumpingrate;
% 
%         b13 = 13;
%         wells(b13, 1) = 247;
%         wells(b13, 2) = 2;
%         wells(b13, 3) = 0;
%         wells(b13, 4) = 0;
%         wells(b13, 5) = 0;
%         wells(b13, 6) = pumpingrate;


        %-------------------------------------------------------------------------


        %         % initially hydraulic charge
        %         h_old=1*ones(size(elem,1),1);
        %         % coeficiente de armazenamento especifico
        %         SS=0.001*1000;
        %         % espesura do aquifero
        %         MM=100;
        %         % step time no considered
        %         dt=1;
        %         % precipitation infiltration
        %         P=0.0009593;
        %         %P=0;
        %         %pumping rate
        %         %pumpingrate=-0;
        %         pumpingrate=-0.0038475;
        %
        %
        %         b1 = 1; % Índice para wells
        %         wells(b1, 1) = 5688;            % Define o índice do poço
        %         wells(b1, 2) = 2;    % Referencia ao elemento correspondente
        %         wells(b1, 3) = 0;             % Valor placeholder (exemplo: coordenada x)
        %         wells(b1, 4) = 0;             % Valor placeholder (exemplo: coordenada y)
        %         wells(b1, 5) = 0;             % Valor placeholder (exemplo: coordenada z)
        %         wells(b1, 6) = pumpingrate;   % Define a taxa de bombeamento
        %
        %         b2 = 2; % Índice para wells
        %         wells(b2, 1) = 9466;
        %         wells(b2, 2) = 2;
        %         wells(b2, 3) = 0;
        %         wells(b2, 4) = 0;
        %         wells(b2, 5) = 0;
        %         wells(b2, 6) = pumpingrate;
        %
        %         b3 = 3; % Índice para wells
        %         wells(b3, 1) = 8343;
        %         wells(b3, 2) = 2;
        %         wells(b3, 3) = 0;
        %         wells(b3, 4) = 0;
        %         wells(b3, 5) = 0;
        %         wells(b3, 6) = pumpingrate;
        %
        %         b4 = 4; % Índice para wells
        %         wells(b4, 1) = 7461;
        %         wells(b4, 2) = 2;
        %         wells(b4, 3) = 0;
        %         wells(b4, 4) = 0;
        %         wells(b4, 5) = 0;
        %         wells(b4, 6) = pumpingrate; % Neste caso, a taxa de bombeamento é negativa
        %
        %         b5 = 5; % Índice para wells
        %         wells(b5, 1) = 7462;
        %         wells(b5, 2) = 2;
        %         wells(b5, 3) = 0;
        %         wells(b5, 4) = 0;
        %         wells(b5, 5) = 0;
        %         wells(b5, 6) = pumpingrate;
        %
        %         b6 = 6; % Índice para wells GAMELEIRA
        %         wells(b6, 1) = 3955;
        %         wells(b6, 2) = 2;
        %         wells(b6, 3) = 0;
        %         wells(b6, 4) = 0;
        %         wells(b6, 5) = 0;
        %         wells(b6, 6) = pumpingrate;
        %
        %         b7 = 7; % Índice para wells LAGOINHA
        %         wells(b7, 1) = 1280;
        %         wells(b7, 2) = 2;
        %         wells(b7, 3) = 0;
        %         wells(b7, 4) = 0;
        %         wells(b7, 5) = 0;
        %         wells(b7, 6) = pumpingrate;
        %
        %         b8 = 8; % Índice para wells PAROIS
        %         wells(b8, 1) = 3510;
        %         wells(b8, 2) = 2;
        %         wells(b8, 3) = 0;
        %         wells(b8, 4) = 0;
        %         wells(b8, 5) = 0;
        %         wells(b8, 6) = pumpingrate;
    case 380
        % article: Numerical modeling of contaminant transport in a stratified
        % heterogeneous aquifer with dipping anisotropy, author QIN 2013
        % coeficiente de armazenamento especifico
        SS=1;
        % espesura do aquifero
        MM=1;
        % step time
    case 380.1
        % well 1
        %wells(1, 1) = 6071;
        wells(1, 1) = 98;
        wells(1, 2) = 2;
        wells(1, 3) = 0;
        wells(1, 4) = 0;
        wells(1, 5) = 0;
        wells(1, 6) = -0.1;
        % well 2
        %wells(2, 1) = 7002;
        wells(2, 1) = 413;
        wells(2, 2) = 2;
        wells(2, 3) = 0;
        wells(2, 4) = 0;
        wells(2, 5) = 0;
        wells(2, 6) = -0.1;
    
end
parmgroundwater.SS=SS;
parmgroundwater.h_init=h_init;
parmgroundwater.MM=MM;
parmgroundwater.dt=dt;
parmgroundwater.P=P;
parmgroundwater.h_hold=h_hold;
%--------------------------------------------------------------------------
source_wells.wells=wells;
end