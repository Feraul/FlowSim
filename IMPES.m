%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter
%Type of file: FUNCTION
%Criate date: 10/01/2012
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%Modified: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals:
%Determinate the saturation and presure fields (2D) in a eithe homogen or
%heterogen domain such as isotropic and anisotropic media for each time
%step or in the steady state when will be important.

%--------------------------------------------------------------------------
%This routine receives geometry and physical data.

%--------------------------------------------------------------------------

function IMPES(Sw,injecelem,producelem,satinbound,wells,klb,satonvertices,...
    satonedges,flagknownvert,flagknownedge,wvector,wmap,constraint,lsw,...
    transmvecleft,transmvecright,knownvecleft,knownvecright,mapinv,...
    maptransm,mapknownvec,pointedge,storeinv,Bleft,Bright,Fg,overedgecoord,...
    bodyterm,normk,limiterflag,massweigmap,othervertexmap,V,N,Hesq,Kde,Kn,...
    Kt,Ded,kmap,nflag,swsequence,ntriang,areatriang,lastimelevel,...
    lastimeval,prodwellbedg,prodwellinedg,mwmaprodelem,vtxmaprodelem,...
    coordmaprodelem,amountofneigvec,rtmd_storepos,rtmd_storeleft,...
    rtmd_storeright,isonbound,elemsize,bedgesize,inedgesize,parameter,...
    weightDMP,nflagface,p_old,contnorm,Kdec,Knc,Ktc,Dedc,weight,s,Con,nflagc,...
    weightc,sc,dparameter,SS,dth,h,MM,gravrate,Dmedio,velmedio,nflagfacec,weightDMPc)
%Define global parameters:
global timew elemarea totaltime timelevel pormap numcase pmethod smethod ...
    filepath benchkey resfolder order inedge;

%--------------------------------------------------------------------------
%Initialize parameters:
c = 0;
gamma=0;
%"time" is a parameter which add "dt" in each looping (dimentional or adm.)
time = lastimeval;
stopcriteria = 0;
%"timelevel" is a parameter used to plot each time step result. Image
%1,2,...,n. This parameter is incremented and sent to "postprocessor"
timelevel = lastimelevel + 1;
%Attribute to time limit ("finaltime") the value put in "Start.dat".
finaltime = totaltime(2);
%"earlyswonedge" is a key used for define the strategy to calculate the
%mobility. See the "getmobility" function to more detail.
earlysw = 0;
timew = 0;
%Initialize "flagtoplot".  Initialy it is "0" and when is "1", plot the
%vtk. It avoids too much vtk files.
flagtoplot = 0;
%"contiterplot" is the number of the vtk created.
contiterplot = 1;

%Parameters to plot
%It is an auxiliary counter. When it get "10", the production parameters
%and saturation field are stored in a file *.dat
countstore = 0;
oilaccum = 0;
S_extrema_old = zeros(1,2);

Sleft = 0;
Sright = 0;

%--------------------------------------------------------------------------
%Verify if there exists a restart

%Verify if the "lastimelevel" is bigger than zero
%In this case, the production parameters must be cought
if lastimelevel > 0
    oilaccum = getrestartdata(producelem);
    %Generate the *.vtk file for the first timestep (initial condition)
else
    %Call the postprocessor:
    postprocessor(ones(elemsize,1),0,Sw,1 - Sw,contiterplot - 1,...
        overedgecoord,order*ones(elemsize,1),'i',1,normk);
end  %End of IF
if numcase==31.1
    [pressure,flowrate,flowresult] = getknownflowrate(elemsize,producelem);
end
tic
while stopcriteria < 100
    
    %User message:
    %Jump a row (in prompt)
    disp(' ');
    disp('---------------------------------------------------');
    disp('>> Show timelevel:')
    timelevel
    %Define Mobility and Saturation on VERTICES and MID-EDGES.
        mobility = getmobility(satinbound,injecelem,Sw,earlysw,smethod,...
            timelevel,numcase,Sleft,Sright,c,overedgecoord);
    %Chose the type of MPFA according "pmethod"
    %Traditional Two-Point Flux Approximation (TPFA), Aziz and Set. (1979)
    if strcmp(pmethod,'tpfa') && numcase~=31.1
        
        %Get "pressure" and "flowrate"
        [pressure,flowrate,flowresult] = solvePressure_TPFA(transmvecleft,...
            knownvecleft,mobility,wells,Fg,bodyterm);
        %MPFA-D (Gao and Wu, 2010)
    elseif strcmp(pmethod,'mpfad') &&  numcase~=31.1
        
        %Calculate "pressure", "flowrate" and "flowresult"
        [pressure,flowrate,flowresult,] = ferncodes_solverpressure(...
            mobility,wells,Hesq,Kde,Kn,Kt,Ded,nflag,weight,s,Con,...
            Kdec,Knc,Ktc,Dedc,nflagc, weightc,sc,SS,dth,h,MM,gravrate);
    elseif strcmp(pmethod,'mpfaql')&& numcase~=31.1
        
        [pressure,flowrate,flowresult,]=ferncodes_solverpressureMPFAQL(nflag,...
            parameter,kmap,weightDMP,wells,mobility,V,Sw,N);
    elseif strcmp(pmethod,'mpfah') && numcase~=31.1% Revisar com cuidado
        
        [pressure,flowrate,flowresult]=ferncodes_solverpressureMPFAH(nflagface,...
            parameter,weightDMP,wells,SS,dth,h,MM,gravrate,mobility);
        
    elseif strcmp(pmethod,'nlfvpp')&& numcase~=31.1
        
    [pressure,flowrate,flowresult,]=ferncodes_solverpressureNLFVPP(nflag,...
    parameter,kmap,wells,mobility,V,Sw,N,p_old,contnorm,weight,s,Con,nflagc,...
                           weightc,sc,0,0,dparameter,SS,dth,h,MM,gravrate);
                                            
        %Any other type of scheme to solve the Pressure Equation
    elseif strcmp(pmethod,'nlfvh')&& numcase~=31.1 % revisar com cuidado
        
        [pressure,flowrate,flowresult]=ferncodes_solverpressureNLFVH(nflagface,...
            parameter,wells,mobility,weightDMP,p_old,0,0,0,contnorm,...
            weightDMPc,Con,nflagfacec,dparameter,weightc,sc,SS,dth,h,MM,gravrate);

    elseif strcmp(pmethod,'nlfvdmp')&& numcase~=31.1
       
        [pressure,flowrate,flowresult]=ferncodes_solverpressureDMP(nflagface,...
            parameter,wells,mobility,weightDMP,p_old,0,0,0,contnorm);
    elseif numcase~=31.1
        
        %Calculate the PRESSURE field (Two-Phase context):
        [pressure,flowrate,flowresult] = solvePressure(transmvecleft,...
            transmvecright,knownvecleft,knownvecright,storeinv,Bleft,...
            Bright,wells,mapinv,maptransm,mapknownvec,pointedge,mobility,...
            bodyterm);
    end  %End of IF (type of pressure solver)
    
    %----------------------------------------------------------------------
    %Calculate "dt" using the function "calctimestep"
    
    %This function obtains the time step using the "Courant" number.
    %The necessity of calculate the time step is ensure the stability of
    %explicit saturation formulation.
    
    dt = calctimestepsat(flowrate,gravrate,Sw,satinbound,injecelem,klb)
    
    %Verify if the "dt" is the last one
    domainvol = sum(elemarea);
    %Get the total flowrate
    if any(producelem)
        totalflowrate = abs(sum(flowresult(producelem)));
    end  %End of IF
    
    %Non-dimensional case:
    if timelevel > 1 && totaltime(1) ~= 0
        booleancond = ((time + dt*totalflowrate/(domainvol*pormap)) > ...
            finaltime);
        %Define again "dt". In this case it will be lower that that
        %calculated.
        dt = booleancond*(finaltime - time)*...
            (domainvol*pormap/totalflowrate) + (1 - booleancond)*dt;
        %Dimensional case:
    elseif totaltime(1) == 0 && (time + dt > finaltime)
        dt = finaltime - time;
    end  %End of IF
    
    %----------------------------------------------------------------------
    
    %Calculate the SATURATION field (choose saturation method):
    [newSw,orderintimestep,waterflowrate,oilflowrate,earlysw,Sleft,Sright] = ...
        solveSaturation(Sw,flowrate,0,dt,injecelem,producelem,satinbound,...
        gravrate,flagknownvert,satonvertices,flagknownedge,satonedges,flowresult,...
        wvector,wmap,constraint,lsw,limiterflag,mobility,massweigmap,...
        othervertexmap,swsequence,ntriang,areatriang,prodwellbedg,...
        prodwellinedg,mwmaprodelem,vtxmaprodelem,coordmaprodelem,...
        amountofneigvec,rtmd_storepos,rtmd_storeleft,rtmd_storeright,...
        isonbound,elemsize,bedgesize,inedgesize,gamma,time);
    
    %Update the saturation field
    Sw = newSw;
    
    %Calculate the OIL SATURATION field:
    %The OIL saturation is obtained using a restriction equation
    So = 1 - newSw;
    
    %----------------------------------------------------------------------
    %Define PVI or DIMENTIONAL TIME
    
    %Dimentional (s, h, day, etc)
    if totaltime(1) == 0
        time = time + dt;
        concluded = time*100/finaltime;
        %It is used for restart activation
        percentdt = dt*100/finaltime;
        stopcriteria = concluded;
        %Define a flag to plot the vtk file
        flagtoplot = flagtoplot + dt*100/finaltime;
        concluded = num2str(concluded);
        status = [concluded '% concluded']
        %VPI (non-dimentional)
    else
        %Increment the parameter "time"
        admtime = ((oilflowrate + waterflowrate)*dt)/(domainvol*pormap);
        time = time + admtime;
        concluded = time*100/finaltime;
        %It is used for restart activation
        percentdt = admtime*100/finaltime;
        stopcriteria = concluded;
        %Define a flag to plot the vtk file
        flagtoplot = flagtoplot + admtime*100/finaltime;
        timew = time/finaltime;
        concluded = num2str(concluded);
        status = [concluded '% concluded']
    end  %End of IF
    
    %----------------------------------------------------------------------
    %Define the WATER, OIL and TIME accumulated
    
    %Update values when there exists producer well(s)
    if any(producelem)
        %Acumulative oil
        oilaccum = oilaccum + (oilflowrate*dt);
        %Water cut for ALL PRODUCER ELEMENT
        watercut = Sw(producelem(1:length(producelem)));
        
        %Write table (Time (VPI), Oil Flow rate, Accumulated Oil and Water Cut)
        %Create the file name
        prfilename = [resfolder '_' 'ProdutionReport.dat'];
        
        %Select the "letter" for "fopen"
        if timelevel == 1
            letter = 'w';
        else
            letter = 'a';
        end  %End of IF
        
        %Open the file
        writeproductionreport = fopen([filepath '\' prfilename],letter);
        
        %Write "productionreport" according to "producelem" size
        %There is one producer well
        if length(producelem) == 1
            fprintf(writeproductionreport,...
                '%26.16E %26.16E %26.16E %26.16E\r\n',...
                [time oilflowrate oilaccum watercut]);
            %There are more than one producer wells (but just write two)
        else
            fprintf(writeproductionreport,...
                '%26.16E %26.16E %26.16E %26.16E %26.16E\r\n',...
                [time oilflowrate oilaccum watercut(1:2)']);
        end  %End of IF
        
        %Close the file "writeproductionreport.dat"
        fclose(writeproductionreport);
    end  %End of IF (when there exists producer well)
    
    %----------------------------------------------------------------------
    %Call the "postprocessor" (plot results in each time step)
    
    %Just create the vtk file if "flagtoplot" reaches 0.1.
    %if flagtoplot >= 1
        %This function create the "*.vtk" file used in VISIT to posprocessing
        %the results
        postprocessor(pressure,flowrate,Sw,So,contiterplot,overedgecoord,...
            orderintimestep,'i',1,normk);
        %Update "flagtoplot"
        flagtoplot = 0;
        %Update "contiterplot"
        contiterplot = contiterplot + 1;
 %   end  %End of IF
    
    %User mesage
    disp('>> Saturation field calculated with success!');
    disp('>> Saturation extrema values [Smax Smin]:');
    %Show extrema values
    S_extrema = [max(Sw); min(Sw)]
    
    %Store maximum and minimum saturation values:
    maxminsatval = [max(S_extrema(1),S_extrema_old(1)) ...
        min(S_extrema(2),S_extrema_old(2))];
    %Update "S_extrema"
    S_extrema_old = S_extrema;
    
    %----------------------------------------------------------------------
    %Write the "restart" file:
    
    %first value --> "timelevel";
    %second value --> "time"
    %from this on --> the saturation field
    
    %Get a boolean condition
    boonwrtfile = strcmp(benchkey,'r')*1 + (1 - strcmp(benchkey,'r'))*100;
    
    %It is actived for each 1% of simulation time
    if countstore > boonwrtfile || stopcriteria == 100
        %Open the file for store the SATURATION FIELD
        writefield = fopen([filepath '\' 'restart.dat'],'w');
        %Write the file
        fprintf(writefield,'%26.16E\r\n',[timelevel; time; Sw]);
        %Close the file "writeresult.dat"
        fclose(writefield);
        
        %Turn "countstore" null
        countstore = 0;
    end  %End of IF
    
    %Increment the parameters "timelevel" and "countstore"
    timelevel = timelevel + 1;
    countstore = countstore + percentdt;
    c=c+1;
    %It gives the time spent per "timelevel"
 
end  %End of While
toc
%Write data file ("ProdutionReport.dat" and others)

plotandwrite(producelem,Sw,pressure,overedgecoord(:,1),injecelem,Dmedio,velmedio,time);

%--------------------------------------------------------------------------
%Write data file ("ProdutionReport.dat" and others)

% plotandwrite(producelem,Sw,pressure,overedgecoord(:,1),injecelem);

% profile off
% profsave(profile('info'),'myprofile_results')

%Mesage for the user:
disp('------------------------------------------------');
disp('>> Global Saturation extrema values [Smax Smin]:');
max_satval = max(maxminsatval(:,1))
min_satval = min(maxminsatval(:,2))
%It deletes the "restart.dat" file
command = ['del ' char(filepath) '\' 'restart.dat'];
%It calls system
system(command);
