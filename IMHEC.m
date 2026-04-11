
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter
%Type of file: FUNCTION
%Criate date: 10/01/2012
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

function IMHEC(Con,injecelem,producelem,satinbound,wells,klb,satonvertices,...
    satonedges,flagknownvert,flagknownedge,wvector,wmap,constraint,lsw,...
    transmvecleft,transmvecright,knownvecleft,knownvecright,mapinv,...
    maptransm,mapknownvec,pointedge,storeinv,Bleft,Bright,Fg,overedgecoord,...
    bodyterm,normk,limiterflag,massweigmap,othervertexmap,V,N,Hesq,Kde,Kn,...
    Kt,Ded,kmap,nflag,swsequence,ntriang,areatriang,lastimelevel,...
    lastimeval,prodwellbedg,prodwellinedg,mwmaprodelem,vtxmaprodelem,...
    coordmaprodelem,amountofneigvec,rtmd_storepos,rtmd_storeleft,...
    rtmd_storeright,isonbound,elemsize,bedgesize,inedgesize,parameter,...
    weightDMP,nflagface,p_old,contnorm,dparameter,nflagc,gamma,...
    Dmedio,Kdec,Knc,Ktc,Dedc,wightc,sc,nflagfacec,weightDMPc,wellsc,...
    weight,s,velmedio,transmvecleftcon,transmvecrightcon,...
    knownvecleftcon,knownvecrightcon,storeinvcon,...
    Bleftcon,Brightcon,Fgcon,mapinvcon,maptransmcon,mapknownveccon,...
    pointedgecon, bodytermcon,gravrate,SS,MM,P,tempo)
%Define global parameters:
global timew  totaltime timelevel numcase pmethod filepath  resfolder ...
    order bcflagc ;

%--------------------------------------------------------------------------
%Initialize parameters:
c = 0;
Conaux=Con;
aux=1;
aux1=1;
mm=0;
nn=0;
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
timew = 0;
%Initialize "flagtoplot".  Initialy it is "0" and when is "1", plot the
%vtk. It avoids too much vtk files.
flagtoplot = 0;
%"contiterplot" is the number of the vtk created.
contiterplot = 1;
earlysw = 0;
%--------------------------------------------------------------------------
%Verify if there exists a restart

%Verify if the "lastimelevel" is bigger than zero
%In this case, the production parameters must be cought
% Alocate sleft and Sright
if lastimelevel~=0
    %Open the restart.dat file:
    command = [char(filepath) '\' 'Results_teste_SleftReport.dat'];
    readfile = fopen(command);

    %"getgeodata" reads each row of *.geo and store it as a string array.
    getdata = textscan(readfile,'%f');
    %Attribute the data to "getgeodata"
    getvecdata = cell2mat(getdata);
    %Fill the variables:

    Sleft = getvecdata(3:length(getvecdata));

    %Open the restart.dat file:
    command = [char(filepath) '\' 'Results_teste_RightReport.dat'];
    readfile = fopen(command);

    %"getgeodata" reads each row of *.geo and store it as a string array.
    getdata = textscan(readfile,'%f');
    %Attribute the data to "getgeodata"
    getvecdata = cell2mat(getdata);
    %Fill the variables:

    Sright = getvecdata(3:length(getvecdata));
else
    Sleft = 0;
    Sright = 0;
end

%Call the postprocessor:
auxkmap=logical(numcase==247)*log(kmap(:,2))+logical(numcase~=247)*normk;
postprocessor_con(ones(elemsize,1),Con,contiterplot - 1,auxkmap);

% mobility
mobility=1;
% calculate the pressure, flow rate advective and dispersive terms
if numcase<379 || numcase==380.1
[hydraulic,flowrateadvec,flowresult,flowratedif]=auxiliarysolverhydraulic(Hesq,Kde,Kn,Kt,Ded,nflag,...
    kmap,V,N,contnorm,weight,s,Con,nflagc,wightc,sc,dparameter,...
    nflagface,parameter,weightDMP,p_old,transmvecleft,...
    transmvecright,knownvecleft,knownvecright,storeinv,Bleft,...
    Bright,wells,mapinv,maptransm,mapknownvec,pointedge,mobility,...
    bodyterm,Kdec,Knc,Ktc,Dedc,weightDMPc,nflagfacec,gravrate,SS,MM,P,tempo);
end
if numcase==241
    %----------------------------------------------------------------------
    %Calculate "dt" using the function "calctimestep"

    %This function obtains the time step using the "Courant" number.
    %The necessity of calculate the time step is ensure the stability of
    %explicit concentration formulation.

    dt = calctimestep(flowrateadvec,satinbound,gamma,Dmedio)
end
%
while stopcriteria < 100

    % while time< finaltime
    %User message:
    %Jump a row (in prompt)
    disp(' ');
    disp('---------------------------------------------------');
    disp('>> Show timelevel:')
    timelevel

    if numcase==246 || numcase==245 || numcase==247 || numcase==248 ||...
            numcase==249 || numcase==250 || numcase==251 || numcase==380 

        if numcase==380 
            %[viscosity] = ferncodes_getviscosity(satinbound,injecelem,Con,earlysw,...
            %    Sleft,Sright,c,overedgecoord,nflagc,nflagfacec);
            viscosity=1;
        else
            [viscosity]=calc_viscosity(nflagfacec,Sleft,Sright,timelevel,...
                earlysw,Con,lastimelevel);
        end
        % To calculate the hydraulic head, advective and dispersive flux
        [hydraulic,flowrateadvec,flowresult,flowratedif]=...
            auxiliarysolverhydraulic(Hesq,Kde,Kn,Kt,Ded,nflag,...
            kmap,V,N,contnorm,weight,s,Con,nflagc,wightc,sc,dparameter,...
            nflagface,parameter,weightDMP,p_old,transmvecleft,...
            transmvecright,knownvecleft,knownvecright,storeinv,Bleft,...
            Bright,wells,mapinv,maptransm,mapknownvec,pointedge,viscosity,...
            bodyterm,Kdec,Knc,Ktc,Dedc,weightDMPc,nflagfacec,gravrate,SS,MM,P,tempo,source);
        %----------------------------------------------------------------------
        %Calculate "dt" using the function "calctimestep"

        %This function obtains the time step using the "Courant" number.
        %The necessity of calculate the time step is ensure the stability of
        %explicit concentration formulation.
        %dt=1;
        %dt = calctimestep(flowrateadvec,satinbound,gamma,Dmedio)
    elseif numcase==241 || numcase==242 || numcase==231 || numcase==232 || numcase==380.1
        viscosity=1;

        % calcula fluxo dispersivo
        [~,~,flowratedif]=auxiliarysolverflux(hydraulic,...
            Con,transmvecleftcon,transmvecrightcon,...
            knownvecleftcon,knownvecrightcon,storeinvcon,Bleftcon,...
            Brightcon,Fgcon,mapinvcon,maptransmcon,mapknownveccon,...
            pointedgecon, bodytermcon,transmvecleft,transmvecright,...
            knownvecleft,knownvecright,storeinv,Bleft,Bright,Fg,...
            mapinv,maptransm,mapknownvec,pointedge, bodyterm,nflag,...
            weight,s,nflagc,wightc,sc,nflagface,weightDMPc,nflagfacec,...
            dparameter,Kde,Ded,Kn,Kt,Hesq,Kdec,Knc,Ktc,Dedc,time,...
            viscosity,parameter,weightDMP,gravrate);
    end
    dt=1;
    %----------------------------------------------------------------------

    %Calculate the CONCENTRATION field (choose concentration method):
    [newC,~,~,~,earlysw,Sleft,Sright] = ...
        solveSaturation(Con,flowrateadvec,flowratedif,dt,injecelem,producelem,satinbound,...
        Fg,flagknownvert,satonvertices,flagknownedge,satonedges,flowresult,...
        wvector,wmap,constraint,lsw,limiterflag,1,massweigmap,...
        othervertexmap,swsequence,ntriang,areatriang,prodwellbedg,...
        prodwellinedg,mwmaprodelem,vtxmaprodelem,coordmaprodelem,...
        amountofneigvec,rtmd_storepos,rtmd_storeleft,rtmd_storeright,...
        isonbound,elemsize,bedgesize,inedgesize,gamma,time);

    %Update the concentration field
    Con = newC;
    %----------------------------------------------------------------------
    if numcase==248
        %Write table (Time , conecntration)
        %Create the file name
        prfilename = [resfolder '_' 'ConReport.dat'];

        %Select the "letter" for "fopen"
        if timelevel == 1
            letter = 'w';
        end  %End of IF

        %Open the file
        writereport = fopen([filepath '\' prfilename],letter);

        A=[timelevel;time+dt; Con];
        fprintf(writereport,'%26.16E \r\n',A);

        %------------------------------------------------------------------
        % Pressure field storage
        %Create the file name
        prfilename = [resfolder '_' 'PresReport.dat'];

        %Open the file
        writereport = fopen([filepath '\' prfilename],letter);
        B=[timelevel;time+dt; hydraulic];
        fprintf(writereport,'%26.16E \r\n',B);
        %Close the file "writeproductionreport.dat"
        fclose(writereport);
        %------------------------------------------------------------------
        % Concentration storage in the cell left
        prfilename = [resfolder '_' 'SleftReport.dat'];

        %Open the file
        writereport = fopen([filepath '\' prfilename],letter);
        CC=[timelevel;time+dt; Sleft];
        fprintf(writereport,'%26.16E \r\n',CC);
        %Close the file "writeproductionreport.dat"
        fclose(writereport);
        %-----------------------------------------------------------------------
        % Concentration storage in the cell right
        prfilename = [resfolder '_' 'RightReport.dat'];

        %Open the file
        writereport = fopen([filepath '\' prfilename],letter);
        D=[timelevel;time+dt; Sright];
        fprintf(writereport,'%26.16E \r\n',D);
        %Close the file "writeproductionreport.dat"
        fclose(writereport);
    end

    %Dimentional (s, h, day, etc)
    disp('>> Time evolution:');
    time = time + dt

    concluded = time*100/finaltime;
    %It is used for restart activation
    percentdt = dt*100/finaltime;
    stopcriteria = concluded;
    %Define a flag to plot the vtk file
    flagtoplot = flagtoplot + dt*100/finaltime;
    concluded = num2str(concluded);
    status = [concluded '% concluded']

    %----------------------------------------------------------------------
    %Call the "postprocessor" (plot results in each time step)

    %Just create the vtk file if "flagtoplot" reaches 0.1.
    if numcase~=380 && numcase~=380.1
        if flagtoplot >= 1
            %This function create the "*.vtk" file used in VISIT to posprocessing
            %the results

            auxkmap=logical(numcase==247)*log(kmap(:,2))+logical(numcase~=247)*normk;

            postprocessor_con(hydraulic,Con,contiterplot,auxkmap);

            flagtoplot = 0;
            %Update "contiterplot"
            contiterplot = contiterplot + 1;
        end  %End of IF
    else
        auxkmap=logical(numcase==247)*log(kmap(:,2))+logical(numcase~=247)*normk;

        postprocessor_con(hydraulic,Con,contiterplot,kmap(:,2));

        flagtoplot = 0;
        %Update "contiterplot"
        contiterplot = contiterplot + 1;

    end

    %User mesage
    disp('>> Concentration field calculated with success!');
    disp('>> Concentration extrema values [Con_max con_min]:');
    %Show extrema values
    C_extrema = [max(Con); min(Con)]

    %Increment the parameters "timelevel" and "countstore"
    timelevel = timelevel + 1;

    c=c+1;
    %It gives the time spent per "timelevel"
    %if time>1000 && (numcase==242 || numcase==243 || numcase==249 || numcase==250)
    if (time>1 && (numcase==242 || numcase==243 || numcase==245))...
            || (time>1 && numcase==380) || (time>200 && numcase==380.1)
        if numcase==245 || numcase==247
            wellsc(1:2,4)=0;
            Con(wells(1,1),1)=wells(1,4);
            Con(wells(2,1),1)=wells(2,4);
            %Initialize and preprocess the parameters:
            [nflagc,nflagfacec] = ferncodes_calflag_con(0);
            %Define elements associated to INJECTOR and PRODUCER wells.
            [injecelem,producelem,satinbound,Conaux,wells] = wellsparameter(wells,...
                Conaux,klb);

            %Define flags and known saturation on the vertices and edges.
            [satonvertices,satonedges,] = ...
                getsatandflag(satinbound,injecelem,Conaux,nflagc,nflagfacec,1);
        elseif time>1 && (numcase==380)
            wellsc(1:20,4)=0;
            %Initialize and preprocess the parameters:
            [nflagc,nflagfacec] = ferncodes_calflag_con(0);
            %Define elements associated to INJECTOR and PRODUCER wells.
            [injecelem,producelem,satinbound,Conaux,wells] = wellsparameter(wells,...
                Conaux,klb);

            %Define flags and known saturation on the vertices and edges.
            [satonvertices,satonedges,] = ...
                getsatandflag(satinbound,injecelem,Conaux,nflagc,nflagfacec,1);
        elseif numcase==380.1
            if time>200
                bcflagc(3,2)=0; % Dirichlet value
                wellsc(1:20,4)=0;
                %Initialize and preprocess the parameters:
                [nflagc,nflagfacec] = ferncodes_calflag_con(0);
                %Define elements associated to INJECTOR and PRODUCER wells.
                [injecelem,producelem,satinbound,Conaux,wells] = wellsparameter(wells,...
                    Conaux,klb);
            end
            
            if Con(wells(1,1))>0 && aux==1
                mm=time;
                aux=aux+1;
            end
            if Con(wells(2,1))>0 && aux1==1
                nn=time;
                aux1=aux1+1;
            end
            

            %Define flags and known saturation on the vertices and edges.
            [satonvertices,satonedges,] = ...
                getsatandflag(satinbound,injecelem,Conaux,nflagc,nflagfacec,1);
        else
            bcflagc(2,2)=0; % Dirichlet value
            %Initialize and preprocess the parameters:
            [nflagc,nflagfacec] = ferncodes_calflag_con(0);
            %Define elements associated to INJECTOR and PRODUCER wells.
            [injecelem,producelem,satinbound,Conaux,wells] = wellsparameter(wells,...
                Conaux,klb);

            %Define flags and known saturation on the vertices and edges.
            [satonvertices,satonedges,] = ...
                getsatandflag(satinbound,injecelem,Conaux,nflagc,nflagfacec,0);

        end
    end

    if numcase==248
        %Initialize and preprocess the parameters:
        nflag = ferncodes_calflag(time);
        [nflagc,nflagfacec] = ferncodes_calflag_con(time);
        %Define flags and known saturation on the vertices and edges.
        [satonvertices,satonedges,flagknownvert,flagknownedge] = ...
            getsatandflag(satinbound,injecelem,Con,nflagc,nflagfacec);
    end
    %[analsol]=anasolaux(velmedio,Dmedio,time);
end  %End of While
toc
%Write data file ("ProdutionReport.dat" and others)

plotandwrite(producelem,Con,hydraulic,satonvertices,Dmedio,velmedio,gamma,time-dt);

%--------------------------------------------------------------------------

% profile off
% profsave(profile('info'),'myprofile_results')

%Mesage for the user:
disp('------------------------------------------------');
disp('>> Global Concentration extrema values [Con_max Con_min]:');
max_conval = max(Con)
min_conval = min(Con)
auxCon=Con(Con<0);
auxCon1=Con(Con>1);
mm
nn
percentCon=((length(auxCon)+length(auxCon1))/length(Con))*100;
%Mesage for the user:
disp('------------------------------------------------');
sprintf('>> Percentagem negative values for the concentratio: %s', num2str(percentCon))

%It deletes the "restart.dat" file
%command = ['del ' char(filepath) '\' 'restart.dat'];
%It calls system
%system(command);
end
%==========================================================================
% we calculate the hydraulic head, flow rate advective and dispersive
function [pressure,flowrateadvec,flowresult,flowratedif]=...
    auxiliarysolverhydraulic(Hesq,Kde,Kn,Kt,Ded,nflag,kmap,V,N,contnorm,...
    weight,s,Con, nflagc,wightc,sc,dparameter,nflagface,parameter,...
    weightDMP,p_old,transmvecleft,transmvecright,knownvecleft,...
    knownvecright,storeinv,Bleft,Bright,wells,mapinv,maptransm,...
    mapknownvec,pointedge,mobility,bodyterm,Kdec,Knc,Ktc,Dedc,...
    weightDMPc,nflagfacec,gravrate,SS,MM,P,tempo,source)

global numcase pmethod
pressure=0;
flowrateadvec=0;
flowresult=0;
flowratedif=0;
dt=0; h=0;

if numcase~=246 & numcase~=246 & numcase~=247 & numcase~=248 & ...
        numcase~=249 & numcase~=250 & numcase~=251

    if strcmp(pmethod,'tpfa')

        %Get "pressure" and "flowrate"
        [pressure,flowrateadvec,flowresult] = solvePressure_TPFA(Kde, Kn,...
            nflag, Hesq,wells,gravrate,source);
    elseif (strcmp(pmethod,'mpfao') || strcmp(pmethod,'fps')) &&  numcase~=31.1
        %Calculate the PRESSURE field (Two-Phase context):
        [pressure,flowrateadvec,flowresult] = solvePressure(transmvecleft,...
            transmvecright,knownvecleft,knownvecright,storeinv,Bleft,...
            Bright,wells,mapinv,maptransm,mapknownvec,pointedge,mobility,...
            bodyterm);
        %MPFA-D (Gao and Wu, 2010)
    elseif strcmp(pmethod,'mpfad')

        %Calculate "pressure", "flowrate" and "flowresult"
        [pressure,flowrateadvec,flowresult,flowratedif] = ferncodes_solverpressure(...
            mobility,wells,Hesq,Kde,Kn,Kt,Ded,nflag,nflagface,weight,s,Con,Kdec,Knc,...
            Ktc,Dedc,nflagc,wightc,sc,SS,dt,h,MM,gravrate,P,kmap,tempo,source);

    elseif strcmp(pmethod,'mpfaql')
        [pressure,flowrateadvec,flowresult]=ferncodes_solverpressureMPFAQL(nflag,...
            parameter,kmap,weightDMP,wells,mobility,V,Con,N,weight,s);

    elseif strcmp(pmethod,'mpfah')

        [pressure,flowrate,]=ferncodes_solverpressureMPFAH(nflagface,...
            parameter,weightDMP,wells,SS,dt,h_old,MM,gravrate,1,P,tempo,source);

    elseif strcmp(pmethod,'nlfvpp')
        [pressure,flowrateadvec,flowresult,flowratedif]=...
            ferncodes_solverpressureNLFVPP(nflag,parameter,kmap,wells,...
            mobility,V,N,p_old,contnorm,weight,s,Con,nflagc,wightc,...
            sc,weightDMPc,nflagfacec,dparameter,SS,dt,h,MM,gravrate,source);

        %Any other type of scheme to solve the Pressure Equation
    elseif strcmp(pmethod,'nlfvh') % revisar com cuidado

        [pressure,flowrateadvec,flowresult,flowratedif]=...
            ferncodes_solverpressureNLFVH(nflagface,...
            parameter,wells,mobility,weightDMP,p_old,weight,s,nflag,...
            contnorm,weightDMPc,Con,nflagfacec,dparameter,wightc,sc);

    elseif strcmp(pmethod,'nlfvdmp')

        [pressure,flowrate,flowresult]=ferncodes_solverpressureDMP(nflagface,...
            parameter,wells,mobility,weightDMP,p_old,0,0,0,contnorm);
    elseif numcase~=31.1

        %Calculate the PRESSURE field (Two-Phase context):
        [pressure,flowrate,flowresult] = solvePressure(transmvecleft,...
            transmvecright,knownvecleft,knownvecright,storeinv,Bleft,...
            Bright,wells,mapinv,maptransm,mapknownvec,pointedge,mobility,...
            bodyterm);
    end  %End of IF (type of pressure solver)
end
end

function [flowrate,flowresult,flowratedif]=auxiliarysolverflux(pressure,...
    Con,transmvecleftcon,transmvecrightcon,...
    knownvecleftcon,knownvecrightcon,storeinvcon,Bleftcon,Brightcon,...
    Fgcon,mapinvcon,maptransmcon,mapknownveccon,pointedgecon, ...
    bodytermcon,transmvecleft,transmvecright,knownvecleft,knownvecright,...
    storeinv,Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge, bodyterm,nflag,...
    weight,s,nflagc,wightc,sc,nflagface,weightDMPc,nflagfacec,...
    dparameter,Kde,Ded,Kn,Kt,Hesq,Kdec,Knc,Ktc,Dedc,time,viscosity,...
    parameter,weightDMP,gravrate)
global  pmethod
if strcmp(pmethod,'nlfvpp')
    % pressure and concentration interpolation
    [pinterp,cinterp]=ferncodes_pressureinterpNLFVPP(pressure,nflag,...
        weight,s,Con,nflagc,wightc,sc);
    % calculate dispersive flux
    [flowrate,flowresult,flowratedif]=ferncodes_flowrateNLFVPP(pressure,...
        pinterp, parameter,viscosity,Con,nflagc,wightc,sc,dparameter,cinterp,gravrate,weightDMPc);
elseif strcmp(pmethod,'nlfvh')
    % pressure and concentration interpoltion on the harmonic points
    [pinterp,cinterp]=ferncodes_pressureinterpHP(pressure,nflagface,parameter,...
        weightDMP,weightDMPc,Con,nflagfacec,dparameter);
    % calculate dispersive flux
    [flowrate,flowresult,flowratedif,]=ferncodes_flowrateNLFVH(pressure,...
        pinterp, parameter,viscosity,Con,cinterp,dparameter);
elseif strcmp(pmethod,'mpfad')
    % pressure and concentration interpolation
    [pinterp,cinterp]=ferncodes_pressureinterpNLFVPP(pressure,nflag,...
        weight,s,Con,nflagc,wightc,sc);
    % calculate dispersive flux
    [flowrate,flowresult,flowratedif] = ferncodes_flowrate(pressure,...
        pinterp,cinterp,Kde,Ded,Kn,Kt,Hesq,viscosity,nflag,Con,...
        Kdec,Knc,Ktc,Dedc,nflagc);
elseif strcmp(pmethod,'mpfao') || strcmp(pmethod,'fps')
    % precisa investir
    [flowrate,flowresult,flowratedif] = calcflowrateMPFAcon(pressure,...
        Con,transmvecleftcon,transmvecrightcon,...
        knownvecleftcon,knownvecrightcon,storeinvcon,Bleftcon,Brightcon,...
        Fgcon,mapinvcon,maptransmcon,mapknownveccon,pointedgecon, ...
        bodytermcon,transmvecleft,transmvecright,knownvecleft,knownvecright,...
        storeinv,Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,viscosity);
else
    % calculate dispersive flux
    [flowrate,flowresult,flowratedif] = calcflowrateTPFA(transmvecleft,...
        pressure,Con,transmvecleftcon,time);
end


end
