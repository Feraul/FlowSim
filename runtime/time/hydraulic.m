%--------------------------------------------------------------------------
%Subject: obtain the hydraulic head
%Type of file: FUNCTION
%Programer: Fernando Contreras,
%--------------------------------------------------------------------------
%This routine receives geometry and physical data.
function hydraulic(wells,overedgecoord,V,N,Hesq,Kde,Kn,Kt,Ded,kmap,nflag,...
    parameter,h_init,contnorm,SS,MM,weight,s,dt,gravrate,...
    nflagface,weightDMP,P,weightDMPc,nflagfacec,weightc,h_old,source)
%Define global parameters:
global timew  totaltime coord pmethod filepath elem numcase inedge ...
    bedge interptype centelem ;

%---------------------------------------------------------------------------
% h: representa a carga hidraulica
% h_new: carga hidraulica atualizado a cada passo de tempo
% h_old: carga hidraulica inicial (a partir da condicao inicial)
%--------------------------------------------------------------------------
%Initialize parameters:
hanalit=0;
haux=0;
%"time" is a parameter which add "dt" in each looping (dimentional or adm.)
theta_s=0;theta_r=0; alpha=0;pp=0; qq=0;hs=0;iterinicial=0;gravresult=0;
flowrateZ=0;flowresultZ=0;
time = 0;
h_kickoff=h_old;
stopcriteria = 0;
orderintimestep = ones(size(elem,1),1)*0;
%Attribute to time limit ("finaltime") the value put in "Start.dat".
finaltime = totaltime(2);
timew = 0;
% inicialization paramenters 
satonvertices=0;producelem=0;h=h_init;Con=0;Kdec=0;Knc=0;nflagc=0;viscosity=1;
count=1;auxkmap=0;mobility=1;Ktc=0;Dedc=0;wightc=0;sc=0;dparameter=0;
% storage the file vtk s in the time 0
postprocessor(h_init,zeros(size(inedge,1)+size(bedge,1),1),Con,1-Con,...
    count,overedgecoord,orderintimestep,'i',1,auxkmap,0);
if numcase==342
    %----------------------------------------------------------------------
    vx=63; % malha quadrilateral ortogonal e distorcida
    %vx=212;% malha triangular nao-estruturada
    haux(1,1)=time;
    haux(1,2)=h_init(vx);
    hanalit(1,1)=time;
    hanalit(1,2)=3*erfc(centelem(vx,1)/(2*sqrt(30.5*time/(3.28*10^-3))));
end
dtaux=dt;
tic
%Get "hydraulic head" and "flowrate"
while stopcriteria < 100
    % utiliza o metodo TPFA para aproximar a carga hidraulica
    if strcmp(pmethod,'tpfa')
        [h_new,flowrate,] = ferncodes_solvePressure_TPFA(Kde, Kn,...
            nflagface, Hesq,wells,viscosity, Kdec, Knc,nflagc,Con,SS,dt,h,...
            MM,P,time);
        % utiliza o metodo MPFA-D para aproximar a carga hidraulica
    elseif strcmp(pmethod,'mpfad')
        % kickoff for non-linear system
        
        % Calculate hydraulic head and flowrate using the MPFA with diamond pacth
        [h_new,flowrate,] = ferncodes_solverpressure(mobility,wells,Hesq,...
            Kde,Kn,Kt,Ded,nflag,nflagface,weight,s,Con,Kdec,Knc,Ktc,Dedc,...
            nflagc,wightc,sc,SS,dt,h,MM,gravrate,P,kmap,time,N,h_kickoff,...
            source,theta_s,theta_r,alpha,pp,qq,iterinicial,gravresult,...
            flowrateZ,flowresultZ);

        % utiliza o metodo MPFA-H para aproximar a carga hidraulica
    elseif strcmp(pmethod,'mpfah')
        % Calculate hydraulic head and flowrate using the MPFA with harmonic
        % points
        [h_new,flowrate,]=ferncodes_solverpressureMPFAH(nflagface,...
            parameter,weightDMP,wells,SS,dtaux,h,MM,gravrate,viscosity,P,time);
        % utiliza o metodo NL-TPFA para aproximar a carga hidraulica
    elseif strcmp(pmethod,'nlfvpp')
        [h_new,flowrate,]=...
            ferncodes_solverpressureNLFVPP(nflag,parameter,kmap,wells,...
            mobility,V,N,p_old,contnorm,weight,s,Con,nflagc,weightc,...
            sc,weightDMPc,nflagfacec,dparameter,SS,dt,h,MM,gravrate,source);
    end
    %% time step calculation
    disp('>> Time evolution:');
    time = time + dt
    % calcula a evoluicao da simulacao em porcentagem
    concluded = time*100/finaltime;
    stopcriteria = concluded;
    concluded = num2str(concluded);
    status = [concluded '% concluded']
    % contador
    count=count+1
    % update the hydraulic head
    %h_init=h;
    h=h_new;
    
    %p_old=h_init;
    % delta t auxiliar
    dtaux=dt;
    %----------------------------------------------------------------------
    % case unconfined aquifer
    % case 1 and 4 of the article Qian, et al 2023
    if numcase==333 || numcase==331 || numcase==347
        if strcmp(pmethod,'mpfah')
            [facelement]=ferncodes_elementfacempfaH;
                kmap = PLUG_kfunction(kmap,h_kickoff,MM,theta_s,theta_r,...
                alpha,pp,qq,hs);
            [pointarmonic]=ferncodes_harmonicopoint(kmap);
            [parameter,]=ferncodes_coefficientmpfaH(facelement,pointarmonic,kmap);
            [weightDMP]=ferncodes_weightnlfvDMP(kmap,elem);
        elseif strcmp(pmethod,'mpfad')

            kmap = PLUG_kfunction(kmap,h_kickoff,MM,theta_s,theta_r,...
                alpha,pp,qq,hs);
            %Get preprocessed terms:
            [Hesq,Kde,Kn,Kt,Ded] = ferncodes_Kde_Ded_Kt_Kn(kmap,elem,theta_r,...
                theta_s,pp,alpha,qq);
            
            %It switches according to "interptype"
            switch char(interptype)
                %LPEW 1
                case 'lpew1'
                    % calculo dos pesos que correspondem ao LPEW1
                    [weight,s] = ferncodes_Pre_LPEW_1(kmap,N);
                    %LPEW 2
                case 'lpew2'
                    % calculo dos pesos que correspondem ao LPEW2
                    [weight,s] = ferncodes_Pre_LPEW_2(kmap,N,zeros(size(elem,1),1));
            end  %End of SWITCH
        else % TPFA
            [~,kmap] = calcnormk(kmap,MM,h);
            %Get preprocessed terms:
            [Hesq,Kde,Kn,Kt,Ded] = ferncodes_Kde_Ded_Kt_Kn(kmap,elem);
        end
    elseif numcase==342 && strcmp(pmethod,'mpfah')
        %----------------------------------------------------------------------
        vx=63; % malha quadrilateral ortogonal distorcida
        %vx=212; % malha triangular nao estruturada
        haux(count+1,1)=log10(time);
        haux(count+1,2)=h_new(vx);
        hanalit(count+1,1)=log10(time);
        hanalit(count+1,2)=3*erfc(centelem(vx,1)/(2*sqrt(30.5*(time)/(3.28*10^-3))));
        % a cada passo de tempo atualiza a condicao de contorno do outro da
        % face a direita
        vv=find(bedge(:,5)==102);
        gg=(coord(bedge(vv,1),:) +coord(bedge(vv,2),:))*0.5;
        %hbound=h(gg);
        hbound=3*erfc(gg(:,1)/(2*sqrt(30.5*(time)/(3.28*10^-3))));
        nflagface(vv,2)=hbound;
    elseif numcase==343
        % calculate the element with coordinate (0.5, 0.5) 
         b2=find((0.4<centelem(:,1)& centelem(:,1)<0.56) & ...
            (0.4<centelem(:,2)& centelem(:,2)<0.56));
         h_time(count,1)=time;
         h_time(count,2)=h(b2);
    end
    
    % storage the vtks and calculate errors
    postprocessor(h,flowrate,Con,1-Con,count,overedgecoord,orderintimestep,...
        'i',1,auxkmap,time);
end

toc
% plotagem dos graficos em determinados regioes do dominio
plotandwrite(producelem,Con,h,satonvertices,0,0,0,time,overedgecoord,...
    hanalit,haux);
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
command = ['del ' char(filepath) '\' 'restart.dat'];
%It calls system
system(command);
