%--------------------------------------------------------------------------
%Subject: obtain the hydraulic head from the Richards' equation
%Type of file: FUNCTION
%Programer: Fernando Contreras,
%--------------------------------------------------------------------------
%This routine receives geometry and physical data.
function hydraulic_RE(env,preTPFA,preMPFAD,preMPFAH,parmRichardEq,source_wells)
% inicializando variaveis globais
centelem=env.geometry.centelem;numcase=env.config.numcase;elem=env.geometry.elem;
inedge=env.geometry.inedge;bedge=env.geometry.bedge;dt= parmRichardEq.dt;
alpha=parmRichardEq.alpha;nvg=parmRichardEq.nvg;normals=env.geometry.normals;
coord=env.geometry.coord;finaltime=env.config.totaltime;elemarea=env.geometry.elemarea;
F=preMPFAD.F;
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
count            =1;
auxcount         =1;
count1           =1;
x                = centelem(:,1);
y                = centelem(:,2);
hnsum=0;

h                = parmRichardEq.h_init;
count_aux = 1; MBE = 0; MBE2 = 0;

h_storage      = [env.geometry.centelem(:,2) h];
time_storage   = 0;
jx             = 1;
sizebedge=size(env.geometry.bedge,1);
sizeinedge=size(env.geometry.inedge,1);
sizebedgeinedge= sizebedge+sizeinedge;

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
elseif env.config.numcase == 436
u_exact(:,count) = -0 .* x.*(1-x) .* y.*(1-y) - 1;

elseif env.config.numcase == 437
    wA=elemarea;
    h_init_exact=h;%carga_hidraulica(centelem(:,1),centelem(:,2),0,parmRichardEq);
    theta_init_exact=thetafunction(h_init_exact,parmRichardEq,env);
    theta_init_num=thetafunction(h,parmRichardEq,env);
    

    totalmassa(auxcount,1)=0;
    totalmassa(auxcount,2)=sum(wA.*theta_init_exact);
    totalmassa(auxcount,3)=sum(wA.*theta_init_num);
    DT_num=theta_init_num;
    DT_exac=theta_init_exact;

    %relative_mass=norm(DT_exac-DT_num)/norm(DT_exac);

    %totalmassa(auxcount,4)=relative_mass;
    %----------------------------------------
end

dtaux        = parmRichardEq.dt;
iterinicial  = 1;
theta_init     = thetafunction(h,parmRichardEq,env);
theta_storage  = [env.geometry.centelem(:,2) theta_init];
kmap_storage   = env.config.perm(:,2);
sum2           = 0;
sum1           = 0;
vnsum          =0;

figure(1); clf;
hold on;
hplot = plot(NaN, NaN, 'o');   % marcador vazio
xlabel('tempo');
ylabel('flowrate(4584)');

tic
%% ================== Loop temporal principal ==================
while stopcriteria < 100

    %% --------- Cálculo da carga hidráulica (h_new) e vazão (flowrate) -----
    if strcmp(env.config.pmethod,'tpfa')
        % TPFA
        [h_new, flowrate] = ...
            ferncodes_solvePressure_TPFA(env,preTPFA,parmgroundwater,parmRichardEq);

    elseif strcmp(env.config.pmethod,'mpfad')
        % MPFA-D (diamond patch)
        [h_new, flowrate, flowresult, flowratedif,faceaux,parmRichardEq,preMPFAD] = ...
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
    time = time + dt
    concluded     = time*100/finaltime;
    stopcriteria  = concluded;

    status = [num2str(concluded) '% concluded']

    count = count + 1

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
        %------------------------------------------------------------------
        v1=bedge(:,1); v2=bedge(:,2);
        m_bedge=0.5*(coord(v1,1:2)+coord(v2,1:2));
        areanormal(1:sizebedge,1)=norm(coord(v1,1:2)-coord(v2,1:2));
        xbedge=m_bedge(:,1); ybedge=m_bedge(:,2);
        u_exact_bedge=-3*t .* xbedge.*(1-xbedge) .* ybedge.*(1-ybedge) - 1;
        % gradiente de p
        dpdx = -3*t .* (1 - 2*xbedge) .* ybedge.*(1-ybedge);
        dpdy = -3*t .* xbedge.*(1-xbedge) .* (1 - 2*ybedge);

        % k(theta(p)) por partes
        ktheta = ones(size(u_exact_bedge));                 % caso p > 0
        mask = u_exact_bedge <= 0;                        % caso p <= 0
        ktheta(mask) = (1 + (-alpha*u_exact_bedge(mask)).^nvg).^(-(nvg-1)/nvg);

        % velocidade v = -k * grad(p)
        vx1 = -ktheta .* dpdx;
        vy1 = -ktheta .* dpdy;

        vn(1:sizebedge,1) = vx1 .* normals(1:sizebedge,1) + vy1 .* normals(1:sizebedge,2);
        %------------------------------------------------------------------
        v11=inedge(:,1); v21=inedge(:,2);
        m_inedge=0.5*(coord(v11,1:2)+coord(v21,1:2));
        areanormal(sizebedge+1:sizebedgeinedge,1)=norm(coord(v11,1:2)-coord(v21,1:2));
        xinedge=m_inedge(:,1); yinedge=m_inedge(:,2);
        % gradiente de p
        dpdx = -3*t .* (1 - 2*xinedge) .* yinedge.*(1-yinedge);
        dpdy = -3*t .* xinedge.*(1-xinedge) .* (1 - 2*yinedge);
        u_exact_inedge = -3*t .* xinedge.*(1-xinedge) .* yinedge.*(1-yinedge) - 1;

        % k(theta(p)) por partes
        ktheta = ones(size(u_exact_inedge));                 % caso p > 0
        mask = u_exact_inedge <= 0;                        % caso p <= 0
        ktheta(mask) = (1 + (-alpha*u_exact_inedge(mask)).^nvg).^(-(nvg-1)/nvg);

        % velocidade v = -k * grad(p)
        vx = -ktheta .* dpdx;
        vy = -ktheta .* dpdy;

        vn(sizebedge+1:sizebedgeinedge,1) = vx .* normals(sizebedge+1:sizebedgeinedge,1) + vy .* normals(sizebedge+1:sizebedgeinedge,2);
        %--------------------------------------------------------------
        vn_numerico=flowrate./areanormal;

        %----------------------------------------------------------------


        wA = elemarea;                         % pesos por elemento
        diffp = u_exact-h;                   % diferença analítica - numérica
        hnsum=hnsum+ dt*(sum(wA .* diffp.^2) / sum(wA) );

        % erro de velocidade
        % pesos Q para cada aresta
        Q = zeros(size(inedge,1) + size(bedge,1), 1);

        Q(1:size(bedge,1)) = elemarea(bedge(:,3));
        % modifique aqui colocando 0.5 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Q(size(bedge,1)+1:end) = 0.5.*(elemarea(inedge(:,3)) + elemarea(inedge(:,4)));

        e=vn-vn_numerico;
        er=e.^2;
        vnsum=vnsum+dt*(Q'*er)/sum(Q');


        postprocessor(h,u_exact,0*parmRichardEq.h_init,...
            time_storage,env,count,parmRichardEq);

        parmRichardEq.h_old=p_old;

        parmRichardEq.h_init=u_exact;
        [source] = PLUG_sourcefunction(0,env,time,parmRichardEq);
        source_wells.source=source;
    elseif numcase == 437

        p_oldaux1 = (h <= 0);
        p_oldaux2 = (h > 0);
        p_old     = 1*p_oldaux2 - 10*p_oldaux1;
        parmRichardEq.h_old=p_old;
        parmRichardEq.h_init=h;
        auxcount=auxcount+1;

        flowresultZ=preMPFAD.flowresultZ;
        %------------------------------------------------------------------
        [h_exact]=carga_hidraulica(centelem(:,1),centelem(:,2),time,parmRichardEq);
        wA = elemarea;                         % pesos por elemento
        diffp = h_exact-h;                     % diferença analítica - numérica
        hnsum=hnsum+ dt*(sum(wA .* diffp.^2) / sum(wA) );
        %------------------------------------------------------------------
        theta_exact=thetafunction(h_exact,parmRichardEq,env);
        theta_num=thetafunction(h,parmRichardEq,env);

        totalmassa(auxcount,1)=time;
        totalmassa(auxcount,2)=sum(wA.*theta_exact);
        totalmassa(auxcount,3)=sum(wA .*theta_num);

        if 0.1<=time
            relativemassa(count1,1)=time;
            DT_num=wA.*(theta_num-theta_init_num);
            DT_exac=wA.*(theta_exact-theta_init_exact);

            relative_mass=norm(DT_exac-DT_num)/norm(DT_exac);

            relativemassa(count1,2)=relative_mass;
            count1=count1+1;
        end
        %------------------------------------------------------------------
        postprocessor(h,h_exact,flowresultZ,time_storage,env,count,...
            parmRichardEq);
        %------------------------------------------------------------------

        [q,q_numerico] = fluxo_Richards(env.geometry.medEdge(:,1),...
            env.geometry.medEdge(:,2), time, parmRichardEq,env,flowrate);

        Q = zeros(size(inedge,1) + size(bedge,1), 1);

        Q(1:size(bedge,1)) = elemarea(bedge(:,3));
        Q(size(bedge,1)+1:end) = 0.5*(elemarea(inedge(:,3)) + elemarea(inedge(:,4)));

        e=q-q_numerico;
        er=e.^2;
        vnsum=vnsum+dt*(Q'*er)/sum(Q');
        %------------------------------------------------------------------
        % Conservacao de massa
        sum2 = sum2+sum(flowrate(1:size(bedge,1),1));
        MBE1(count_aux,1)  = time;
        MBE1(count_aux,2)  = (sum(wA.*theta_num-wA.*theta_init_num) ...
            - abs(sum2)*dt) / (sum(wA.*theta_init_num)) ;

        MBER1(count_aux,1) = time;
        MBER1(count_aux,2) = (sum(wA.*theta_num-wA.*theta_init_num)) ...
            /(abs(sum2)*dt);

        count_aux=count_aux+1;
    elseif numcase==438
        p_oldaux1 = (h <1);
        p_oldaux2 = (h >= 1);
        p_old     = 10*p_oldaux2 - 10*p_oldaux1;
        parmRichardEq.h_old=p_old;
        parmRichardEq.h_init=h;
        flowresultZ=preMPFAD.flowresultZ;
         %------------------------------------------------------------------
        postprocessor(h,0*h,0*flowresultZ,time_storage,env,count,...
            parmRichardEq);
    elseif numcase==439
        p_oldaux1 = (h <0);
        p_oldaux2 = (h >= 0);
        p_old     = 20*p_oldaux2 - 30*p_oldaux1;
        parmRichardEq.h_old=p_old;
        parmRichardEq.h_init=h;
        flowresultZ=preMPFAD.flowresultZ;
         %------------------------------------------------------------------
         h1=(h-(200-centelem(:,2)));
        postprocessor(h,h1,0*flowresultZ,time_storage,env,count,...
            parmRichardEq);
       %-------------------------------------------------------------------
       mmmm=bedge(:,5)==202;

       facescontor=find(mmmm==1);

       elementosedge=bedge(facescontor,3);

       faceselem=F(elementosedge,:);
       
       mask1=faceselem<size(bedge,1);
       
       sss=max(abs(flowrate(faceselem(~mask1))));

        %env.config.bcflag(3,2)=min(14.8, sss);
        %------------------------------------------------------------------
        
        set(hplot, 'XData', [get(hplot,'XData') time], ...
           'YData', [get(hplot,'YData') sss]);
        drawnow;
    end

    %% --------- Atualizações específicas para MPFA-D ---------------------
    if strcmp(env.config.pmethod,'mpfad')
        [env,parmRichardEq] = PLUG_kfunction(env, parmRichardEq, time);

        [preMPFAD] = ferncodes_Kde_Ded_Kt_Kn(env,parmRichardEq,preMPFAD);

        [preMPFAD,~,~] = ferncodes_Pre_LPEW_2_vect(preMPFAD,parmRichardEq,env);
        if numcase==436
            [env] = ferncodes_calflag(env,parmRichardEq,time);
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

    faceaux(:,1) = faceaux(:,1) .* (stopcriteria >= 100);

end  % while
if numcase==436
    wA = elemarea;                         % pesos por elemento
    diffp = u_exact-h;                   % diferença analítica - numérica
    hnsum1=(sum(wA .* diffp.^2) / sum(wA) );

    errortotal1=(hnsum+vnsum) % Analysis of an Euler Implicit-Mixed Finite Element
    % Scheme for Reactive Solute Transport in Porous Media, veja o paragrafo da
    % equacao 5.5, dh=dt.

    errortotal2=(hnsum1+vnsum) % artigo:Convergence analysis for a mixed finite element scheme for flow in
    %strictly unsaturated porous media, ultimo paragrafo da equacao 37, dh=dt.
    errorateconv(u_exact, h, vn,vn_numerico,env)
elseif numcase==437

    wA = elemarea;                         % pesos por elemento
    diffp = h_exact-h;                   % diferença analítica - numérica
    hnsum1=(sum(wA .* diffp.^2) / sum(wA) );

    errortotal1=(hnsum+vnsum) % Analysis of an Euler Implicit-Mixed Finite Element
    % Scheme for Reactive Solute Transport in Porous Media, veja o paragrafo da
    % equacao 5.5, dh=dt.

    errortotal2=(hnsum1+vnsum) % artigo:Convergence analysis for a mixed finite element scheme for flow in
    %strictly unsaturated porous media, ultimo paragrafo da equacao 37, dh=dt.

    figure(1)
    plot(totalmassa(:,1),totalmassa(:,2),'k-')

    hold on
    plot(totalmassa(:,1),totalmassa(:,3),'b-')
    legend('Exact solution','Numerical solution')
    xlabel('t')
    ylabel('Total mass')
    figure(2)
    plot(relativemassa(:,1),relativemassa(:,2))

    figure(3)
    plot(MBE1(:,1), MBE1(:,2),'g-')

    figure(4)
    plot(MBER1(:,1), MBER1(:,2),'r-')

    wA = elemarea;                         % pesos por elemento

    difftheta=theta_exact-theta_num;

    thetasum=(sum(wA .* difftheta.^2) / sum(wA) );

    diffp = h_exact-h;                     % diferença analítica - numérica
    hnsum1=sqrt(sum(wA .* diffp.^2) / sum(h_exact.^2.*wA) );

    Q = zeros(size(inedge,1) + size(bedge,1), 1);

    Q(1:size(bedge,1)) = elemarea(bedge(:,3));
    Q(size(bedge,1)+1:end) = 0.5*(elemarea(inedge(:,3)) + elemarea(inedge(:,4)));

    e=q-q_numerico;
    er=e.^2;
    vnsum=sqrt((Q'*er)/sum(Q));
end
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
elseif numcase==439
inedge=env.geometry.inedge;

Lef = inedge(:,3);
Rel = inedge(:,4);

% faces onde há mudança de sinal
mask = h(Lef).*h(Rel) < 0;

% elementos candidatos
Lef_neg = h(Lef) > 0;

% vetor final de elementos
elemento = Lef;
elemento(~Lef_neg) = Rel(~Lef_neg);

% aplicar a máscara
elemento = elemento(mask);

centros=centelem(elemento,:);
A = sortrows(centros,1);
 % remover duplicados em x
[xu, ia] = unique(A(:,1));
zu = A(ia,2);

% interpolação
xq = linspace(min(xu), max(xu), 500);
yq = interp1(xu, zu, xq, 'spline');

% plot
plot(A(:,1),A(:,2),'o')
hold on
plot(xq,yq,'r','LineWidth',2)
grid on
legend('Datos','Curva interpolada')
figure(2)

        %t=2
    B=[0.00000	78.7
5.45455	78.5
14.1176	78.3307
20.2139	78.3307
27.9144	78.0096
34.0107	77.6886
39.7861	77.0465
49.7326	76.0835
60.0000	75.4414
65.7754	74.7994
73.7968	74.1573
81.4973	73.1942
88.2353	72.5522
95.9358	71.9101
104.920	70.9470
114.545	69.6629
123.529	68.3788
133.155	67.4157
144.706	66.7737
159.465	66.4527
174.866	65.8106
199.251	65.0
222.674	65.0
247.059	65.0
300.321	65.0
];

    plot(A(:,1), A(:,2),'k')
    hold on
    plot(B(:,1),B(:,2),'b')

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
