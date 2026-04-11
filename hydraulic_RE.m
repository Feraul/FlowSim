%--------------------------------------------------------------------------
%Subject: obtain the hydraulic head from the Richards' equation
%Type of file: FUNCTION
%Programer: Fernando Contreras,
%--------------------------------------------------------------------------
%This routine receives geometry and physical data.
function hydraulic_RE(env,preTPFA,preMPFAD,preMPFAH,parmRichardEq,source_wells)
% inicializando variaveis globais
centelem=env.geometry.centelem;
numcase=env.config.numcase;
elem=env.geometry.elem;
inedge=env.geometry.inedge;
bedge=env.geometry.bedge;
dt= parmRichardEq.dt;
finaltime=env.config.totaltime;

% determina a altura 75 onde esta localizado a face
afaceaux=find(abs(env.geometry.coord(env.geometry.inedge(:,1),2)-60)<1e-9);
facestop=afaceaux+size(env.geometry.bedge,1);
% os elementos que conforman ate a acima da altura 75
elemaux=find(env.geometry.centelem(:,2)>75);


%---------------------------------------------------------------------------
% h: representa a carga hidraulica
% h_new: carga hidraulica atualizado a cada passo de tempo
% h_old: carga hidraulica inicial (a partir da condicao inicial)
%--------------------------------------------------------------------------
%% ================== Inicialização de parâmetros ==================
hanalit          = 0;
haux             = 0;
time             = 0;              % tempo (dimensional ou adimensional)
stopcriteria     = 0;
orderintimestep  = zeros(size(env.geometry.elem,1),1);
timew            = 0;
count            =1;
zero             =zeros(size(env.geometry.elem,1),1);
x = centelem(:,1);
y = centelem(:,2);

h                = parmRichardEq.h_init;
count_aux = 1; MBE = 0; MBE2 = 0;

h_storage      = [env.geometry.centelem(:,2) h];
time_storage   = 0;
jx             = 1;
sizebedgeinedge= size(env.geometry.bedge,1)+size(env.geometry.inedge,1);
% Armazena VTK no tempo 0
postprocessor(h,0*ones(sizebedgeinedge,1),0*parmRichardEq.h_init,...
    time_storage,env,1,parmRichardEq)

%% ================== Casos especiais por numcase ==================
if env.config.numcase == 432
    auxelem22 = env.goemetry.bedge(57:80,3);
    plot(parmRichardEq.h_init(auxelem22), env.config.perm(auxelem22,2))
    xlabel('h(m)')
    ylabel('Hidraulic conductivity')

elseif env.config.numcase == 431
    mmm   = (size(env.geometry.bedge,1)/2) + 1;
    sumax = 0;
    sumax1 = 0;
end

dtaux        = parmRichardEq.dt;
iterinicial  = 1;
theta_init     = thetafunction(h,parmRichardEq,env);
theta_storage  = [env.geometry.centelem(:,2) theta_init];
kmap_storage   = env.config.perm(:,2);
sum2           = 0;
sum1           = 0;

tic
u_exact(:,count) = -0 .* x.*(1-x) .* y.*(1-y) - 1;
%% ================== Loop temporal principal ==================
while stopcriteria < 100

    %% --------- Cálculo da carga hidráulica (h_new) e vazão (flowrate) -----
    if strcmp(env.config.pmethod,'tpfa')
        % TPFA
        [h_new, flowrate] = ...
            ferncodes_solvePressure_TPFA(env,preTPFA,parmgroundwater,parmRichardEq);

    elseif strcmp(env.config.pmethod,'mpfad')
        % MPFA-D (diamond patch)
        [h_new, flowrate, flowresult, flowratedif,faceaux,parmRichardEq] = ...
            ferncodes_solverpressure(env,preMPFAD,parmRichardEq,dtaux,...
            source_wells,time );

        if count == 1
            sum1 = sum(flowresult);
        end
    elseif strcmp(env.config.pmethod,'mpfah')
        % MPFA-H (harmonic points)
        [h_new, flowrate] = ...
            ferncodes_solverpressureMPFAH(env,preMPFAH,parmgroundwater,parmRichardEq);
    end

    % --------- Atualização da permeabilidade kmap_h -----------------------
    kmap_storage(:,2*count-1:2*count) = [centelem(:,2) parmRichardEq.auxperm(:,2)];

    % --------- Atualização do tempo e critério de parada ------------------
    disp('>> Time evolution:');
    time = time + dt;
    concluded     = time*100/finaltime;
    stopcriteria  = concluded;

    status = [num2str(concluded) '% concluded']

    count = count + 1;

    % Atualiza carga hidráulica
    h     = h_new;
    dtaux = dt;

    %% --------- Pós-processamento e armazenamento --------------------------


    h_storage(:,2*count-1:2*count) = [centelem(:,2) h];

    theta_n = thetafunction(h,parmRichardEq,env);
    theta_storage(:,2*count-1:2*count) = [centelem(:,2) theta_n];

    time_storage(count,1) = time;

    %% ================== Tratamento por numcase ==================
    if numcase == 431

        postprocessor(h,0*parmRichardEq.h_init,0*parmRichardEq.h_init,...
            time_storage,env,count,parmRichardEq);
        % chute inicial para o método iterativo
        p_oldaux1 = (h <= 0);
        p_oldaux2 = (h > 0);
        p_old     = 100*p_oldaux2 - 100*p_oldaux1;
        parmRichardEq.h_old=p_old;
        parmRichardEq.h_init=h;
        % cálculo do balanço de massa
        pqp = 100/size(elem,1);

        % Conservação de massa em altura definida
        sum2 = sum2 + (-flowrate(mmm,1) + flowrate(facestop,1));
        MBE1(count_aux,1)  = time;
        MBE1(count_aux,2)  = ( (sum(theta_n(elemaux)-theta_init(elemaux))*pqp ...
            - sum2*dt) / (sum(theta_init(elemaux))*pqp) );

        MBER1(count_aux,1) = time;
        MBER1(count_aux,2) = (sum(theta_n(elemaux)-theta_init(elemaux))*pqp) ...
            /(abs(sum2)*dt);

        % Conservação de massa global
        sum1 = sum1 + (-flowrate(mmm,1) + flowrate(1,1));
        MBE2(count_aux,1)  = time;
        MBE2(count_aux,2)  = (sum(theta_n-theta_init)*pqp - sum1*dt) ...
            /(sum(theta_init)*pqp);

        MBER2(count_aux,1) = time;
        MBER2(count_aux,2) = (sum(theta_n-theta_init)*pqp)/(sum1*dt);

        count_aux = count_aux + 1;

        % Fluxos de contorno
        flux_contour_inter(jx,1) = time;
        flux_contour_inter(jx,2) = abs(flowrate(mmm));
        flux_contour_inter(jx,3) = abs(flowrate(facestop));

        sumax = sumax + abs(flux_contour_inter(jx,2))*dt;
        kmap_h=parmRichardEq.auxperm;
        kmapaux=env.config.perm;
        % Verifica se atingiu altura desejada
        if  kmap_h(inedge(afaceaux,3),2) == kmapaux(1,2) && ...
                kmap_h(inedge(afaceaux,4),2) == kmapaux(1,2)

            store_influx_contourflux_and_time(1,1) = time;
            store_influx_contourflux_and_time(1,2) = abs(flowrate(facestop));
            store_influx_contourflux_and_time(1,3) = abs(flowrate(mmm));

            store_acumulativeflux_time(1,1) = time;
            store_acumulativeflux_time(1,2) = sumax;

            MBEunico  = MBE2;
            MBERunico = MBER2;
            break

        elseif stopcriteria >= 100
            store_influx_contourflux_and_time(1,1) = time;
            store_influx_contourflux_and_time(1,2) = abs(flowrate(facestop));
            store_influx_contourflux_and_time(1,3) = abs(flowrate(mmm));

            store_acumulativeflux_time(1,1) = time;
            store_acumulativeflux_time(1,2) = sumax;

            MBEunico  = MBE2;
            MBERunico = MBER2;
        end

        jx = jx + 1;

    elseif numcase == 433
        p_oldaux1 = (h <= 0);
        p_oldaux2 = (h > 0);
        p_old     = 25*p_oldaux2 - 25*p_oldaux1;
        parmRichardEq.h_old=p_old;

    elseif numcase == 432
        if time <= 1      % Beit Netofa clay
            bcflag(3,2) = -2 + (2.2*(time/1));
        else
            bcflag(3,2) = 0.2;
        end

        % Atualiza condição de Dirichlet
        [nflag, nflagface] = ferncodes_calflag(0);

        p_oldaux1 = (h <= 0);
        p_oldaux2 = (h > 0);
        p_old     = 2*p_oldaux2 - 2*p_oldaux1;
        parmRichardEq.h_old=p_old;

    elseif numcase == 434
        p_oldaux1 = (h <= 0);
        p_oldaux2 = (h > 0);
        p_old     = 50*p_oldaux2 - 50*p_oldaux1;
        parmRichardEq.h_old=p_old;

    elseif numcase == 435
        p_oldaux1 = (h <= 0);
        p_oldaux2 = (h > 0);
        p_old     = 5*p_oldaux2 - 5*p_oldaux1;
        parmRichardEq.h_old=p_old;

    elseif numcase == 436
        p_oldaux1 = (h <= 0);
        p_oldaux2 = (h > 0);
        p_old     = 10*p_oldaux2 - 5*p_oldaux1;
        
        %------------------------------------------------------------------
        x = centelem(:,1);
        y = centelem(:,2);
        t = time;
        u_exact = -3*t .* x.*(1-x) .* y.*(1-y) - 1;

        postprocessor(h,u_exact,0*parmRichardEq.h_init,...
            time_storage,env,count,parmRichardEq);

        parmRichardEq.h_old=p_old;

        parmRichardEq.h_init=u_exact;
        [source] = PLUG_sourcefunction(0,env,time,parmRichardEq);
        source_wells.source=source;
    end

    %% --------- Atualizações específicas para MPFA-D ---------------------
    if strcmp(env.config.pmethod,'mpfad')
        [env,parmRichardEq] = PLUG_kfunction(env, parmRichardEq, time);

        [preMPFAD] = ferncodes_Kde_Ded_Kt_Kn(env,parmRichardEq,preMPFAD,time);

        [preMPFAD,~,~] = ferncodes_Pre_LPEW_2_vect(zero,preMPFAD,parmRichardEq,env);
        if numcase==436
            [nflag, nflagface] = ferncodes_calflag(env,time);
            preMPFAD.nflag=nflag;
            preMPFAD.nflagface=nflagface;
        end
    end

    iterinicial = iterinicial + 1;

    %% --------- Pós-processamento extra e gráficos -------------------------
    %postprocessor(h, flowrate, Con, 1-Con, count, ...
    %              overedgecoord, orderintimestep, ...
    %              'i', 1, kmap_h(:,2), time);
    if numcase == 432 && count == 2
        auxelem22 = bedge(49:80,3);
        plot(h(auxelem22), kmap(auxelem22,2))
        xlabel('h(m)')
        ylabel('hidraulic conductivity')
    end

    if stopcriteria < 100
        faceaux(:,1) = 0;
    end

end  % while
errorateconv(u_exact, h, 0,0,env)
if numcase==435
    % theta=thetafunction(h,theta_s,theta_r,alpha,pp,q);
    % theta_init=thetafunction(h_init,theta_s,theta_r,alpha,pp,q);
    % figure(1)
    % plot(centelem(:,2),theta)
    % xlabel('Elevation (cm)')
    % ylabel('\theta(h)')
    %
    % hold on
    % plot(centelem(:,2),theta_init )
    % grid
    % figure(2)
    % plot(centelem(:,2),h)
    % grid
    % xlabel('Elevantion (cm)')
    % ylabel('h(cm)')
    % hold on
elseif numcase==431
    % theta_n=thetafunction(h,theta_s,theta_r,alpha,pp,q);
    % theta_init=thetafunction(h_init,theta_s,theta_r,alpha,pp,q);

    %MBE1=(sum(theta_n-theta_init)*0.25)/(sum2*dt);
    % figure(1)
    % plot(auxflow(:,1), auxflow(:,2))
    % hold on
    % plot(auxflow(:,1), auxflow(:,3),'o')
    % figure(2)
    % plot(auxflowresult(:,1), auxflowresult(:,2))
    % figure(1)
    % theta=thetafunction(h,theta_s,theta_r,alpha,pp,q);
    % plot(theta,centelem(:,2))
    % xlabel('\theta (cm)')
    % ylabel('Elevation')
    % grid
    %
    % figure(2)
    % plot(h,centelem(:,2))
    % hold on

    % Warrick solution 11700
    % BBB=[-343.291	73.1783
    % -170.127	74.5736
    % -78.9873	77.5194];

    % Warrick solution 23400
    % BBB=[-345.291	60.7752
    % -170.127	62.3256
    % -78.9873	65.8915];
    % Warrick solution 46800

    % BBB=[-343.291	39.8450
    % -170.127	41.7054
    % -78.9873	45.8915];
    % plot(BBB(:,1),BBB(:,2), 's')
    % xlabel('h (cm)')
    % ylabel('Elevation(cm)')
    % hold on
    % grid
    % figure(1)
    % plot(MBE(:,1), MBE(:,2))
    % xlabel('Time step (min)')
    % ylabel('Massa Balance Error')
    % grid
end
toc
% plotagem dos graficos em determinados regioes do dominio
%plotandwrite(producelem,Con,h,satonvertices,0,0,0,time,overedgecoord,...
%    hanalit,haux);
%--------------------------------------------------------------------------
% activate if you want to visualize the hydraulic field in the final time

%h

%--------------------------------------------------------------------------
% profile off
% profsave(profile('info'),'myprofile_results')

%Mesage for the user:
disp('------------------------------------------------');
disp('>> Global Hydraulic head extrema values [hmax hmin]:');
max_hyval = max(h)
min_hyval = min(h)
% if numcase==342
%     erro=norm(haux(:,2)-hanalit(:,2))
% elseif numcase==343
%    plot(h_time(:,1),h_time(:,2))
%    hold on
%    grid
% end
%It deletes the "restart.dat" file
%command = ['del ' char(filepath) '\' 'restart.dat'];
if numcase==431
    filepath=env.mainpathfolders.path;
    tabfolder=env.mainpathfolders.tabfolder;

    fname = [filepath '\' tabfolder '\'];

    writematrix(h_storage,[fname 'h_steptime3.txt'])

    writematrix(theta_storage, [fname 'WaterContent_steptime3.txt'])

    writematrix(centelem,[fname 'centrocell3.txt'])

    writematrix(time_storage,[fname 'time_step3.txt'])

    writematrix(kmap_storage,[fname 'condhydraulic_steptime3.txt'])

    writematrix(flux_contour_inter(:,1:3),[fname 'flux_contour_inter_steptime3.txt'])
end
%It calls system
%system(command);
