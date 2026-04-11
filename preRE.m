function [parmRichardEq,env]=preRE(env)
% inicializacao de parametros globais
elem=env.geometry.elem;
hs=0;
theta_s=0;
theta_r=0;
nvg=0;
pp=0;
q=0;
env.config.numcase=env.config.numcase+30;
switch env.config.numcase
    case 431
        % Dados obtidos de Caviedes, malha z=100, 400 ok
        %ks=[1 0.0001 0 0 0.0001];
        theta_s=0.363;
        theta_r=0.186;
        alpha=0.01;
        pp=1.53;
        q=1-(1/pp);
        dt=195 % segundo a tese de Caviedes
        h_init=-800*ones(size(elem,1),1);
        h_old(1:size(elem,1),1)=-100;
        %------------------------------------------------------------------
        % theta_s=0.5;
        % theta_r=0.05;
        % alpha=0.01;
        % pp=1.67;
        % q=1-(1/pp);
        % h_init=-211.3413*ones(size(env.geometry.elem,1),1);
        % dt=1; 
        % h_old(1:size(env.geometry.elem,1),1)=-100;
    case 432
        % silt loam
        % theta_s=0.396;
        % theta_r=0.131;
        % alpha=0.423;
        % pp=2.06;
        % q=1-(1/pp);
        % dt=0.0208;       
        %------------------------------------------------------------------
        % Beit Netofa clay
        theta_s=0.446;
        theta_r=0.0;
        alpha=0.152;
        pp=1.17;
        q=1-(1/pp);
        dt=0.3;
        % valor incial
        for i=1:size(env.geometry.elem,1)
            h_init(i,1) =1-env.geometry.centelem(i,2);
        end
       
        % chute inicial
        for i=1:size(env.geometry.centelem,1)
            if env.geometry.centelem(i,2)<1
                h_old(i,1)=2.5;
            else
                h_old(i,1)=-2.5;
            end

        end

    case 433
        theta_s=0.43;
        theta_r=0.078;
        alpha=0.036;
        pp=1.56;
        q=1-(1/pp);
        % K= 24.96
        %flux 20  CC. Neumann
        h_init=-51.3949*ones(size(env.geometry.elem,1),1);
        dt=0.05;
    case 434
        theta_s=0.3;
        theta_r=0.01;
        alpha=0.033;
        pp=4.1;
        q=1-(1/pp);
        for i=1:size(env.geometry.elem,1)
            h_init(i,1) =65-env.geometry.centelem(i,2);
        end
        dt=5;
    case 435
        %------------------------------------------------------------------
        theta_s=0.43;
        theta_r=0.078;
        alpha=0.036;
        pp=1.56;
        % K=0.0001
        % h=20  CC. Dirichlet
        q=1-(1/pp);
        for i=1:size(env.geometry.elem,1)
            a =10-env.geometry.centelem(i,2);
            if a>0
                h_init(i,1)=10;
                h_old(i,1)=1;
            else
                h_init(i,1)=-90;
                h_old(i,1)=-1;
            end
        end
        dt=0.03;
        % chute inicial
    case 436
       % adaptado de
       %MPFA methods for Richards' equation
       %Truls Moholt
        alpha=0.1844;
        nvg=3;
        
        h_init=-1*ones(size(env.geometry.elem,1),1);
        h_old=-5*ones(size(env.geometry.elem,1),1);
        dt= (1/32)^2;
        % chute inicial

end
parmRichardEq.theta_s=theta_s;
parmRichardEq.theta_r=theta_r;
parmRichardEq.alpha=alpha;
parmRichardEq.pp=pp;
parmRichardEq.q=q;
parmRichardEq.h_init=h_init;
parmRichardEq.dt=dt;
parmRichardEq.h_old=h_old;
parmRichardEq.hs=hs;
parmRichardEq.nvg=nvg;
end