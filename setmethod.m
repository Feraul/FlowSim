
%Programer: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals: this FUNCTION maneger the kind of pressure, headr hydraulic,
% solute concentration, water saturation solver
%--------------------------------------------------------------------------
%Additional comments:
%--------------------------------------------------------------------------

function setmethod(source_wells,keywrite,invh,env,parmconcentra, parmgroundwater, parmRichardEq)

%Get a preprocessment of pressure scheme (used in One-phase and Two-Phase).
if env.config.phasekey ~= 0
    %Call "preMPFA". This function calculate some important parameters to
    %be used in any PRESSURE solver.
    [env,preTPFA,preMPFAD,preMPFAQL,preNLTPFA, preMPFAH, preconcentraMPFAD,...
        preconcentraNLTPFA, preGravity,parmRichardEq] = preMPFA(env,parmgroundwater,...
        parmRichardEq);

end  %End of IF (execute "preMPFA")
env.config.phasekey = env.config.phasekey;

%According "phasekey" the "One-phase" or "Two-phase" procedures are choose.
switch env.config.phasekey
    case 1
        pmethod=env.config.pmethod;
        solvers = struct( ...
            'tpfa',   @( ) solvePressure_TPFA(preMPFAD), ...
            'mpfad',  @( ) ferncodes_solverpressure(env,preMPFAD,parmRichardEq,0,source_wells,0), ...
            'mpfaql', @( ) ferncodes_solverpressureMPFAQL(nflag,parameter,kmap,weightDMP,wells,1,V,1,N,weight,s), ...
            'mpfah',  @( ) ferncodes_solverpressureMPFAH(nflagface,parameter,weightDMP,wells,SS,dt,h_init,MM,gravrate,P), ...
            'nlfvpp', @( ) ferncodes_solverpressureNLFVPP(nflagno,parameter,kmap,wells,1,V,1,N,p_old,contnorm), ...
            'nlfvh',  @( ) ferncodes_solverpressureNLFVH(nflagface,parameter,kmap,wells,1,weightDMP,p_old,contnorm), ...
            'nlfvdmp',@( ) ferncodes_solverpressureDMP(nflagface,parameter,wells,1,weightDMP,p_old,0,0,0) ...
            );

        if isfield(solvers, pmethod)
            [pressure, flowrate] = solvers.(pmethod)();
        else
            [pressure, flowrate] = solvePressure(transmvecleft,transmvecright,knownvecleft, ...
                knownvecright,storeinv,Bleft,Bright,wells,mapinv,maptransm, ...
                mapknownvec,pointedge,1,bodyterm);
        end
        max(pressure)
        min(pressure) 
        postprocessor(pressure,flowrate,0, 0,env,1,parmRichardEq)
    case 4 % hydraulic head simulation
        %% ===========================================================
        % steady-state problem
        if ismember(env.config.numcase, [336, 334, 335,337,338,339,340,341,347,341.1])
            solvers = struct( ...
                'tpfa',   @( ) solvePressure_TPFA(preTPFA), ...
                'mpfad',  @( ) ferncodes_solverpressure(env,preMPFAD, preconcentraMPFAD, parmconcentra, parmgroundwater, preGravity), ...
                'mpfaql', @( ) ferncodes_solverpressureMPFAQL(env,preMPFAQL), ...
                'mpfah',  @( ) ferncodes_solverpressureMPFAH(env,preMPFAH, parmgroundwater, preGravity), ...
                'nlfvpp', @( ) ferncodes_solverpressureNLFVPP(env, preNLTPFA));
            time=0;
            % Verifica se o método existe na tabela
            if isfield(solvers, env.config.pmethod)
                [pressure, flowrate] = solvers.(env.config.pmethod)();
            end
            %Plot the fields (pressure, normal velocity, etc)
            %This function create the "*.vtk" file used in VISIT to
            %postprocessing the results
            postprocessor(pressure,flowrate,0,1,1,overedgecoord,1,keywrite,...
                invh,normk,0);
            
            if env.config.numcase==333
                plotandwrite(0,0,pressure,0,0,0,0,0,overedgecoord);
            end

            %Mesage for the user:
            disp('------------------------------------------------');
            disp('>> Global hydraulic head extrema values:');
            max_conval = max(pressure)
            min_conval = min(pressure)

        else
            %% ===============================================================
            % transient-state problem
            hydraulic(wells,parmRichardEq);
        end

    case 5 % contaminant transport with hydarulic head
        tempo=0;
        %Define elements associated to INJECTOR and PRODUCER wells.
        [injecelem,producelem,satinbound,Con,wellsc] = wellsparameter(wellsc,...
            Con,klb);

        %Define flags and known concentration or saturation on the vertices and edges.
        [satonvertices,satonedges,flagknownvert,flagknownedge] = ...
            getsatandflag(satinbound,injecelem,Con,nflagnoc,nflagfacec,0);

        %"precsaturation" - Preprocessor of the concentration or saturation equation
        [wvector,wmap,constraint,massweigmap,othervertexmap,lsw,swsequence,...
            ntriang,areatriang,prodwellbedg,prodwellinedg,mwmaprodelem,...
            vtxmaprodelem,coordmaprodelem,amountofneigvec,rtmd_storepos,...
            rtmd_storeleft,rtmd_storeright,isonbound] = ...
            preSaturation(flagknownedge,injecelem,producelem);

        %"IMHEC" function. There, HYDRAULIC HEAD and CONCENTRATION are solved.
        IMHEC(wells,keywrite,invh,env,parmconcentra,env,preTPFA,preMPFAD,...
            preMPFAQL,preNLTPFA, preMPFAH,presturation,parmRichardEq);
    case 6
        % steady-state problem

        %% ===============================================================
        % transient-state problem
        hydraulic_RE(env,preTPFA,preMPFAD,preMPFAH,parmRichardEq,source_wells);
end  %End of SWITCH
end


