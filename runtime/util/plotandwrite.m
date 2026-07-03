%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter
%Programer: Márcio Souza
%Modified: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals:
%This function writes and plots the data calculated.

%--------------------------------------------------------------------------
%Additional Comments: It is called by the function "IMPES.m"


%--------------------------------------------------------------------------

function [analsw, presanal]=plotandwrite(producelem,Sw,pressure,...
    satonvertices,Dmedio,velmedio,gamma,time,overedgecoord,hanalitica,hauxiliar)
%Define global parameters:
global filepath satlimit resfolder bcflag bcflagc numcase totaltime ...
    elemarea modflowcompared;

%"getgeodata" reads each row of *.geo and store it as a string array.
%Get the data according to "producelem" size
%There is one producer well
%Create the file name
prfilename = [resfolder '_' 'ProdutionReport.dat'];
%Open the file
readfile = fopen([filepath '\' prfilename]);

if numcase==245 || numcase==247

    %"getgeodata" reads each row of *.geo and store it as a string array.
    %Get the data according to "producelem" size
    %There is one producer well
    if length(producelem) == 1
        getdata = textscan(readfile,'%f%f');
        %Stabilishes number of columns
        colnum = 2;
        %There are more than one producer wells (but we just write two)
    else
        getdata = textscan(readfile,'%f%f');
        %Stabilishes number of columns
        colnum = 2;
    end  %End of IF

    %Attribute the data to "getgeodata"
    getmatxdata = cell2mat(getdata);
    %Initialize prodution parameters:
    contime = getmatxdata(:,1);
    oilflowratevec = getmatxdata(:,2);
    %    oilaccumvec = getmatxdata(:,3);
    %Attribute watercut column
    %    watercutvec = getmatxdata(:,colnum);
elseif 30<numcase && numcase<100
    %"getgeodata" reads each row of *.geo and store it as a string array.
    %Get the data according to "producelem" size
    %There is one producer well
    if length(producelem) == 1
        getdata = textscan(readfile,'%f%f%f%f');
        %Stabilishes number of columns
        colnum = 4;
        %There are more than one producer wells (but we just write two)
    else
        getdata = textscan(readfile,'%f%f%f%f%f');
        %Stabilishes number of columns
        colnum = 4:5;
    end  %End of IF
    %Attribute the data to "getgeodata"
    getmatxdata = cell2mat(getdata);
    %Initialize prodution parameters:
    contime = getmatxdata(:,1);
    oilflowratevec = getmatxdata(:,2);
    oilaccumvec = getmatxdata(:,3);
    %Attribute watercut column
    watercutvec = getmatxdata(:,colnum);
end
%--------------------------------------------------------------------------
%PLOT RESULTS

%Buckley-Leverett Problem (Case 31)
if numcase >= 31 && numcase <= 31.9
    %Get Buckley-Leverett parameters
    %Define the flowrate value (to semi-analitical solution)
    Qvalue = bcflag(logical(bcflag(:,1) > 200 & bcflag(:,2) ~= 0),2);
    %Get saturation and position to plot

    [posit,satfield,elemonline] = getlineresult(Sw,satonvertices,producelem);
    %Calculate semi-analitical saturation field.

    [x,Swanal] = solanalBL(Qvalue,totaltime(2),elemonline,Sw);

    %------------------------------------------------------
    %Define the polynimial of interpolation. Bastian (2002)
    %     xanal = flipud(x(2:length(x)));
    %     Sat = fliplr(Swanal(2:length(Swanal)));
    %     blcurve = fit(xanal,Sat','cubicinterp');
    %------------------------------------------------------

    %Evaluate the errors and Convergence Rate:
    %     buckey_levalidation(x(1),elemonline,Sw,blcurve);

    %Store the results
    storeresult(posit,satfield);
    hold on
    %Plot the results (Analitical Solution)
    plot (x,Swanal,'LineWidth',2)
    %Plot the results (Actual Numerical Solution)
    hold on;
    plot(posit,satfield,'-bv');
    %View "analytical curve"
    %     plot(xanal,blcurve(xanal),'-r','LineWidth',2)
    hold on;

    grid on;
    title ('Buckley-Leverett Solution')
    xlabel ('Reservoir Length');
    ylabel ('Water Saturation (Sw)');

    %Define the limits according to "numcase" number
    %Bastian (2002) problem
    if numcase == 31.1
        xlim ([0 300]);
        %Define a limit for "y" on the graph
        ylim ([satlimit(1) (1 - satlimit(2))]);
        %Serna (2009)
    elseif numcase == 31.8
        xlim ([-1 1]);
        %Define a limit for "y" on the graph
        ylim ([-0.1 1.1]);
        %Yu-Shu Wu (2017), Livro "Multiphase Fluid Flow in Porous and Fract..."
    elseif numcase == 31.9
        xlim ([0 120]);
        %Define a limit for "y" on the graph
        ylim ([satlimit(1) 1]);
        %Any other problem
    else
        xlim ([0 1]);
        %Define a limit for "y" on the graph
        ylim ([satlimit(1) (1 - satlimit(2))]);
    end  %End of IF

    %Kozdon's problems (case 45 and their sub-cases)
elseif numcase >= 45 && numcase < 46
    %Plot Water Cut
    plot(contime,watercutvec(:,1),'--g');
    hold on;
    plot(contime,watercutvec(:,2),'-g');
    hold off;

    grid on;
    xlabel('Time, [VPI]');
    ylabel('Water Cut, [adm]');
    title('Water Cut in Producer Well');
    %    legend('Upwind (Sat.) and TPS (Pres.)');
    xlim([0 0.06]);
    ylim([0 0.31]);

    %Eymard's problems (case 46 and their sub-cases)
elseif numcase >= 46 && numcase < 47
    %Get saturation and position to plot
    [satfieldh,satfield,posith,positd] = getlineresultEymardProblem(Sw);
    %Get pressure and position to plot
    [pressfieldh,pressfield,] = getlineresultEymardProblem(pressure);

    %Plot Water Saturation
    figure(1)
    %Water Saturation Profile (horizontal)
    plot(posith,satfieldh,'-k');
    hold on;
    %Water Saturation Profile (diagonal)
    plot(positd,satfield,'--b');
    hold off;

    grid on;
    xlabel('Domain');
    ylabel('Water Saturation');
    title('Water Saturarion (horizontal and diagonal profiles)');
    legend('Water Saturation (profile along horizontal line)',...
        'Water Saturation (profile along diagonal line)');
    xlim([-0.8 0.8]);
    ylim([0 1.301]);

    %Plot Pressure
    figure(2)
    %Pressure Profile (horizontal)
    plot(posith,pressfieldh,'-k');
    hold on;
    %Pressure Profile (diagonal)
    plot(positd,pressfield,'--b');
    hold off;

    grid on;
    xlabel('Domain');
    ylabel('Pressure');
    title('Pressure (horizontal and diagonal profiles)');
    legend('Pressure (profile along horizontal line)',...
        'Pressure (profile along diagonal line)');
    xlim([-0.8 0.8]);
    ylim([0 45]);

    %Khoozan's problems (case 48 and their sub-cases)
elseif numcase > 48 && numcase < 49
    %Plot Oil FLOW RATE, CUMULATIVE OIL and Water Cut
    figure(1);
    plot(contime,oilflowratevec,'-k');

    grid on;
    xlabel('Time, [VPI]');
    ylabel('Oil Recovery, [m^3/s]');
    title('Oil Recovery');
    legend('Upwind (Sat.) and TPS (Pres.)');
    xlim([0 1])
    %     ylim([0 1.21])

    %Plot the Cumulative Oil
    figure(2);
    plot(contime,oilaccumvec,'-b');
    grid on;
    xlabel('Time, [VPI]');
    ylabel('Cumulative Oil, [m^3]');
    title('Cumulative Oil');
    legend('Upwind (Sat.) and TPS (Pres.)');
    xlim([0 1])
    %     ylim([0 1])

    %Plot Water Cut
    figure(3);
    %Well 1:
    plot(contime,watercutvec(:,1),'--b');
    hold on;
    %Well 2:
    plot(contime,watercutvec(:,2),'-b');
    hold off;

    grid on;
    xlabel('Time, [VPI]');
    ylabel('Water Cut, [adm]');
    title('Water Cut in Producer Well');
    %    legend('Upwind (Sat.) and TPS (Pres.)');
    %     xlim([0 0.06]);
    %     ylim([0 0.31]);

    %Another two-phase flow cases (cases 32 on)
elseif numcase > 200
    switch numcase
        case 231
            [posit,confield,elemonline] = getlineresult(Sw,satonvertices);
            [analsolution]=ferncodes_analyticalSolution(satonvertices,...
                Dmedio,velmedio,posit,gamma);
            % readequacao
            analsolutionaux=length(analsolution)+1;
            confieldaux= length(analsolution)+1;
            positaux= length(posit)+1;
            analsolutionaux(1)=bcflagc(2,2); % concentratio na fronteira;
            confieldaux(1)=bcflagc(2,2);
            positaux(1)=0;

            analsolutionaux(2:length(analsolution)+1,1)=analsolution;
            confieldaux(2:length(confield)+1,1)=confield;
            positaux(2:length(posit)+1,1)=posit;

            hold on

            %Plot the results (Analitical Solution)
            plot (positaux,analsolutionaux,'k','LineWidth',2)
            %Plot the results (Actual Numerical Solution)
            %hold on;

            plot(positaux,confieldaux,'-b','LineWidth',1.5);

            hold on;

            grid on;
            title ('Concentration Solution')
            xlabel ('Reservoir Length');
            ylabel ('Concentration (C)');
            % Calculo dos erros

            %Calculate the relative error of "el2" ("relerrorL2")
            relerrorL2 = (abs(confield - analsolution).^2).*elemarea(elemonline,:);

            relerrorL1=  (abs(confield - analsolution).^1).*elemarea(elemonline,:);
            abserrorMAX = abs(analsolution - confield);

            %Calculate "emax"
            emaxc = max(abserrorMAX)
            %calculate "el1"
            el1c = sum(relerrorL1)/sum(elemarea(elemonline,:))
            %Calculate "el2"
            el2c = sqrt(sum(relerrorL2)/sum(elemarea(elemonline,:)))


        case 232
            [posit,confield,elemonline] = getlineresult(Sw,satonvertices);
            [analsolution]=ferncodes_analyticalSolution(satonvertices,Dmedio,velmedio,posit,gamma);
            % readequacao
            analsolutionaux=length(analsolution)+1;
            confieldaux= length(analsolution)+1;
            positaux= length(posit)+1;
            analsolutionaux(1)=bcflagc(2,2); % concentratio na fronteira;
            confieldaux(1)=bcflagc(2,2);
            positaux(1)=0;

            analsolutionaux(2:length(analsolution)+1,1)=analsolution;
            confieldaux(2:length(confield)+1,1)=confield;
            positaux(2:length(posit)+1,1)=posit;

            hold on

            %Plot the results (Analitical Solution)
            plot (positaux,analsolutionaux,'k','LineWidth',2)
            %Plot the results (Actual Numerical Solution)
            %hold on;
            plot(positaux,confieldaux,'-r','LineWidth',1.5);
            %plot(positaux,confieldaux,'Color',[0,0.7,0.9],'LineWidth',1.5);
            %plot(positaux,confieldaux,'Color', [0.9290 0.6940 0.1250],'LineWidth',1.5);

            %View "analytical curve"
            %     plot(xanal,blcurve(xanal),'-r','LineWidth',2)
            hold on;

            grid on;
            title ('Concentration Solution')
            xlabel ('Reservoir Length');
            ylabel ('Concentration (C)');
            % Calculo dos erros

            %Calculate the relative error of "el2" ("relerrorL2")
            relerrorL2 = (abs(confield - analsolution).^2).*elemarea(elemonline,:);

            relerrorL1=  (abs(confield - analsolution).^1).*elemarea(elemonline,:);
            abserrorMAX = abs(analsolution - confield);

            %Calculate "emax"
            emaxc = max(abserrorMAX)
            %calculate "el1"
            el1c = sum(relerrorL1)/sum(elemarea(elemonline,:))
            %Calculate "el2"
            el2c = sqrt(sum(relerrorL2)/sum(elemarea(elemonline,:)))


        case 233
            [posit,confield,elemonline] = getlineresult(Sw,satonvertices);
            [analsolution]=ferncodes_analyticalSolution(satonvertices,...
                Dmedio,velmedio,posit,gamma);
            hold on
            %Plot the results (Analitical Solution)
            plot (posit,analsolution,'LineWidth',2)
            %Plot the results (Actual Numerical Solution)
            hold on;
            plot(posit,confield,'-mv');
            %View "analytical curve"
            %     plot(xanal,blcurve(xanal),'-r','LineWidth',2)
            hold on;

            grid on;
            title ('Concentration Solution')
            xlabel ('Reservoir Length');
            ylabel ('Concentration (C)');
        case 234
            [posit,confield,elemonline] = getlineresult(Sw,satonvertices);
            [analsolution]=ferncodes_analyticalSolution(satonvertices,...
                Dmedio,velmedio,posit,gamma);
            hold on
            %Plot the results (Analitical Solution)
            plot (posit,analsolution,'LineWidth',2)
            %Plot the results (Actual Numerical Solution)
            hold on;
            plot(posit,confield,'-bv');
            %View "analytical curve"
            %     plot(xanal,blcurve(xanal),'-r','LineWidth',2)
            hold on;

            grid on;
            title ('Concentration Solution')
            xlabel ('Reservoir Length');
            ylabel ('Concentration (C)');
        case 235
            [posit,confield,elemonline] = getlineresult(Sw,satonvertices);
            [analsolution]=ferncodes_analyticalSolution(satonvertices,...
                Dmedio,velmedio,posit,gamma);
            hold on
            %Plot the results (Analitical Solution)
            plot (posit,analsolution,'LineWidth',2)
            %Plot the results (Actual Numerical Solution)
            hold on;
            plot(posit,confield,'-gv');
            %View "analytical curve"
            %     plot(xanal,blcurve(xanal),'-r','LineWidth',2)
            hold on;

            grid on;
            title ('Concentration Solution')
            xlabel ('Reservoir Length');
            ylabel ('Concentration (C)');
        case 236
            [posit,confield,elemonline] = getlineresult(Sw,satonvertices);
            [analsolution]=ferncodes_analyticalSolution(satonvertices,Dmedio,velmedio,posit,gamma);
            hold on
            %Plot the results (Analitical Solution)
            plot (posit,analsolution,'LineWidth',2)
            %Plot the results (Actual Numerical Solution)
            hold on;
            plot(posit,confield,'-gv');
            %View "analytical curve"
            %     plot(xanal,blcurve(xanal),'-r','LineWidth',2)
            hold on;

            grid on;
            title ('Concentration Solution')
            xlabel ('Reservoir Length');
            ylabel ('Concentration (C)');
        case 237
            [posit,confield,elemonline] = getlineresult(Sw,satonvertices);
            [analsolution]=ferncodes_analyticalSolution(satonvertices,Dmedio,velmedio,posit,gamma);
            hold on
            %Plot the results (Analitical Solution)
            plot (posit,analsolution,'LineWidth',2)
            %Plot the results (Actual Numerical Solution)
            hold on;
            plot(posit,confield,'-gv');
            %View "analytical curve"
            %     plot(xanal,blcurve(xanal),'-r','LineWidth',2)
            hold on;

            grid on;
            title ('Concentration Solution')
            xlabel ('Reservoir Length');
            ylabel ('Concentration (C)');
        case 238
            [posit,confield,elemonline] = getlineresult(Sw,satonvertices);
            [analsolution]=ferncodes_analyticalSolution(satonvertices,Dmedio,velmedio,posit,gamma);
            hold on
            %Plot the results (Analitical Solution)
            plot (posit,analsolution,'LineWidth',2)
            %Plot the results (Actual Numerical Solution)
            hold on;
            plot(posit,confield,'-gs');
            %View "analytical curve"
            %     plot(xanal,blcurve(xanal),'-r','LineWidth',2)
            hold on;

            grid on;
            title ('Concentration Solution')
            xlabel ('Reservoir Length');
            ylabel ('Concentration (C)');
        case 241
            [posit,confield,elemonline] = getlineresult(Sw,satonvertices);
            [analsolution]=ferncodes_analyticalSolution(satonvertices,Dmedio,velmedio,posit,gamma);

            analsolutionaux=length(analsolution)+1;
            confieldaux= length(analsolution)+1;
            positaux= length(posit)+1;
            analsolutionaux(1)=1; % concentratio na fronteira;
            confieldaux(1)=1;
            positaux(1)=0;

            analsolutionaux(2:length(analsolution)+1,1)=analsolution;
            confieldaux(2:length(confield)+1,1)=confield;
            positaux(2:length(posit)+1,1)=posit;
            % hold on
            %Plot the results (Analitical Solution)
            plot (positaux,analsolutionaux,'k','LineWidth',2)
            %Plot the results (Actual Numerical Solution)
            hold on;
            plot(positaux,confieldaux,'b','LineWidth',1.5);
            %View "analytical curve"
            %     plot(xanal,blcurve(xanal),'-r','LineWidth',2)
            hold on;

            grid on;
            title ('Concentration Solution')
            xlabel ('Reservoir Length');
            ylabel ('Concentration (C)');
        case 242
            [posit,confield,elemonline] = getlineresult(Sw,satonvertices);
            [analsolution]=ferncodes_analyticalSolution(satonvertices,Dmedio,velmedio,posit,gamma);

            % readequacao ative para malha grossa dx=0.5

            confieldaux= length(analsolution)+1;
            positaux= length(posit)+1;

            confieldaux(1)=0;
            positaux(1)=0;


            confieldaux(2:length(confield)+1,1)=confield;
            positaux(2:length(posit)+1,1)=posit;

            confield=confieldaux;
            posit=positaux;
            % hold on
            %Plot the results (Analitical Solution)
            %  plot (posit,analsolution,'k','LineWidth',2)
            %Plot the results (Actual Numerical Solution)
            hold on;
            plot(posit,confield,'b','LineWidth',1.5);
            %View "analytical curve"
            %     plot(xanal,blcurve(xanal),'-r','LineWidth',2)
            hold on;

            grid on;
            title ('Concentration Solution')
            xlabel ('Reservoir Length');
            ylabel ('Concentration (C)');
        case 248

            [analsw, presanal]=ferncodes_analyticalSolution(satonvertices,Dmedio,velmedio,0,gamma,time);
            abserrorMAX = abs(analsw - Sw);
            %Calculate the relative error of "el2" ("relerrorL2")
            relerrorL2 = (abs(Sw - analsw).^2).*elemarea;
            relerrorL2P= (abs(pressure - presanal).^2).*elemarea;
            relerrorL1=  (abs(Sw - analsw).^1).*elemarea;
            %Calculate "emax"
            emax = max(abserrorMAX)
            %calculate "el1"
            el1 = sum(relerrorL1)/sum(elemarea)
            %Calculate "el2"
            el2c = sqrt(sum(relerrorL2)/sum(elemarea))
            el2p = sqrt(sum(relerrorL2P)/sum(elemarea))
            %el2 = sqrt(sum(relerrorL2))
        case 245
            %Plot Oil FLOW RATE, CUMULATIVE OIL and Water Cut
            figure(1);
            plot(contime,oilflowratevec,'--k');

            grid on;
            xlabel('Time, [VPI]');
            ylabel('Oil Recovery, [m^3/s]');
            title('Oil Recovery');
            legend('Upwind (Sat.) and TPS (Pres.)');
            %xlim([0 1])
            %ylim([0 1.21])
            hold on
            %Plot the Cumulative Oil
            figure(2);
            plot(contime,oilaccumvec,'--k');
            grid on;
            xlabel('Time, [VPI]');
            ylabel('Cumulative Oil, [m^3]');
            title('Cumulative Oil');
            legend('Upwind (Sat.) and TPS (Pres.)');
            %xlim([0 1])
            %ylim([0 1])
            hold on
            %Plot Water Cut
            figure(3);
            plot(contime,watercutvec(:,1)+watercutvec(:,2),'--b');
            hold on
            grid on;
            xlabel('Time, [VPI]');
            ylabel('Water Cut');
            title('Water Cut');
            legend('Upwind (Sat.) and TPS (Pres.)');
            %xlim([0 1])
        case 247
            figure(1);
            hold on
            plot(contime,oilflowratevec,'--b');

            grid on;
            xlabel('Time, [VPI]');
            ylabel('Concenration, [m^3/s]');
            title('Oil Recovery');
            legend('Upwind (Sat.) and TPS (Pres.)');
            %xlim([0 1])
            %ylim([0 1.21])
            hold on
        case 330
            % solucao numerica numa determinada regiao do dominio
            [posit,confield,] = getlineresult(pressure,satonvertices);

            if strcmp(modflowcompared,'n')
                confieldaux= zeros(length(confield)+2,1);
                positaux= zeros(length(posit)+2,1);

                confieldaux(1)=100;
                positaux(1)=0;
                confieldaux(end)=100;
                positaux(end)=1000;

                confieldaux(2:length(confield)+1,1)=confield;
                positaux(2:length(posit)+1,1)=posit;

                plot(positaux,confieldaux,'--b','LineWidth',1.5);
            else
                plot(posit,confield,'--b','LineWidth',1.5);
            end
            hold on
            grid
            xlabel('y(m)');
            ylabel('Hydraulic head (m)');

        case 331
            % solucao numerica numa determinada regiao do dominio
            [posit,confield,] = getlineresult(pressure,satonvertices);
            if strcmp(modflowcompared,'n')
                % neste caso as condicoes de contorno esta imposta na
                % fronteita da malha
                confieldaux= zeros(length(confield)+2,1);
                positaux= zeros(length(posit)+2,1);

                confieldaux(1)=90;
                positaux(1)=0;
                confieldaux(end)=90;
                positaux(end)=1000;

                confieldaux(2:length(confield)+1,1)=confield;
                positaux(2:length(posit)+1,1)=posit;

                plot(positaux,confieldaux,'--b','LineWidth',1.5);
            else
                plot(posit,confield,'--b','LineWidth',1.5);
            end
            hold on
            grid
            xlabel('y(m)');
            ylabel('Hydraulic head (m)');

        case 332
            % solucao numerica numa determinada regiao do dominio
            [posit,confield,] = getlineresult(pressure,satonvertices);
            if strcmp(modflowcompared,'n')
                % neste caso as condicoes de contorno esta imposta na
                % fronteita da malha
                confieldaux= zeros(length(confield)+2,1);
                positaux= zeros(length(posit)+2,1);

                confieldaux(1)=100;
                positaux(1)=0;
                confieldaux(end)=100;
                positaux(end)=1000;

                confieldaux(2:length(confield)+1,1)=confield;
                positaux(2:length(posit)+1,1)=posit;

                plot(positaux,confieldaux,'--g','LineWidth',1.5);
            else
                plot(posit,confield,'--g','LineWidth',1.5);
            end
            hold on
            grid
            xlabel('y(m)');
            ylabel('Hydraulic head (m)');

        case 333
            % solucao analitica
            [presanalit,] = benchmark(overedgecoord);
            % solucao numerica
            [posit,confield,elemonline] = getlineresult(pressure,satonvertices);
            if strcmp(modflowcompared,'n')
                % neste caso as condicoes de contorno esta imposta na
                % fronteita da malha
                % numerical solution
                confieldaux= zeros(length(confield)+2,1);
                positaux= zeros(length(posit)+2,1);

                confieldaux(1)=2;
                positaux(1)=0;
                confieldaux(end)=2;
                positaux(end)=40;

                confieldaux(2:length(confield)+1,1)=confield;
                positaux(2:length(posit)+1,1)=posit;

                % analytical solution
                confieldaux1= zeros(length(elemonline)+2,1);
                positaux1= zeros(length(posit)+2,1);

                confieldaux1(1)=2;
                positaux1(1)=0;
                confieldaux1(end)=2;
                positaux1(end)=40;

                confieldaux1(2:length(confield)+1,1)=presanalit(elemonline);
                positaux1(2:length(posit)+1,1)=posit;
                %---------------------------------------------
                plot(positaux1,confieldaux1,'k','LineWidth',1.)
                hold on
                %------------------------------------------------------
                plot(positaux,confieldaux,'-b','LineWidth',1.5);
                hold on
            else
                %---------------------------------------------
                plot(posit,presanalit(elemonline),'k','LineWidth',1.)
                hold on
                %------------------------------------------------------
                plot(posit,confield,'b','LineWidth',1.5);
                hold on
            end
            grid
            xlabel('y(m)');
            ylabel('Hydraulic head (m)');
        case 342
            % solucao numerica
            [posit,confield,elemonline] = getlineresult(pressure,satonvertices);

            %---------------------------------------------
            %plot(posit,presanalit(elemonline),'k','LineWidth',1.)
            %hold on
            %------------------------------------------------------
            figure(1)
            plot(posit,confield,'-<','LineWidth',1.5);
            hold on
            grid
            figure(2)
            plot(hauxiliar(:,1), hauxiliar(:,2),'k')
            hold on
            plot(hanalitica(:,1), hanalitica(:,2),'r')
            grid
            xlabel('y(m)');
            ylabel('Hydraulic head (m)');
        case 343
              % solucao analitica
            [presanalit,] = benchmark(overedgecoord,time);
            % solucao numerica
            [posit,confield,elemonline] = getlineresult(pressure,satonvertices);
            
            % neste caso as condicoes de contorno esta imposta na
            % fronteita da malha
            % numerical solution
                confieldaux= zeros(length(confield)+2,1);
                positaux= zeros(length(posit)+2,1);

                confieldaux(1)=10;
                positaux(1)=0;
                confieldaux(end)=2;
                positaux(end)=1;

                confieldaux(2:length(confield)+1,1)=confield;
                positaux(2:length(posit)+1,1)=posit;

                % analytical solution
                confieldaux1= zeros(length(elemonline)+2,1);
                positaux1= zeros(length(posit)+2,1);

                confieldaux1(1)=10;
                positaux1(1)=0;
                confieldaux1(end)=2;
                positaux1(end)=1;

                confieldaux1(2:length(confield)+1,1)=presanalit(elemonline);
                positaux1(2:length(posit)+1,1)=posit;


            figure(1)
            plot(positaux1,confieldaux,'k-<','LineWidth',1.5);
            hold on
            plot(positaux1, confieldaux1,'k','LineWidth',1.5)
            grid
            xlabel('x(m)');
            ylabel('Hydraulic head (m)');

    end

else
    %Plot Oil FLOW RATE, CUMULATIVE OIL and Water Cut
    figure(1);
    plot(contime,oilflowratevec,'--k');

    grid on;
    xlabel('Time, [VPI]');
    ylabel('Oil Recovery, [m^3/s]');
    title('Oil Recovery');
    legend('Upwind (Sat.) and TPS (Pres.)');
    xlim([0 1])
    ylim([0 1.21])
    hold on
    %Plot the Cumulative Oil
    figure(2);
    plot(contime,oilaccumvec,'--k');
    grid on;
    xlabel('Time, [VPI]');
    ylabel('Cumulative Oil, [m^3]');
    title('Cumulative Oil');
    legend('Upwind (Sat.) and TPS (Pres.)');
    xlim([0 1])
    ylim([0 1])
    hold on
    %Plot Water Cut
    figure(3);
    plot(contime,watercutvec(:,1),'--k');
    hold on
    grid on;
    xlabel('Time, [VPI]');
    ylabel('Water Cut');
    title('Water Cut');
    legend('Upwind (Sat.) and TPS (Pres.)');
    xlim([0 1])

    %     ylim([0 1.21])
end  %End of IF

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTIONS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "getlineresult"
%--------------------------------------------------------------------------

function [posit,satfield,elemonline] = getlineresult(Sw,satonvertices,producelem)
%Define global parameters:
global centelem numcase coord satlimit;



%Initialize "pos" and "satfield"
%getxvalue(1) = 0;
getsatfield(1) = max(satonvertices);
getelemonline(1) = 0;
%Initialize "j"
j = 1;
jj=2;
%Swept "centelem"
for i = 1:size(centelem,1)
    if numcase==242 || numcase==241
        % ative para malha grossa
        %ymax= 2.7;
        %ymin= 2.5;
        % native para malha fina
        ymax= 2.4;
        ymin= 2.2;
        %ymax=1.72;
        %ymin=1.58;
    elseif numcase==231 || numcase==232
        % malha utilizado na tese de Nilson e Uewerton.
        % malha11
        ymax=max(coord(:,2));
        %ymax=10;
        %ymin=ymax/2;
        ymin=0;
    else
        ymin=0;
        ymax=0;
    end
    if 200<numcase && ((ymin<=centelem(i,2)) && (centelem(i,2) <= ymax))
        %if (2.2<=centelem(i,2)) && (centelem(i,2) <= 2.4)
        %Attribute to "pos" the value of "centelem" which match with
        %"y_value"
        getxvalue(j) = centelem(i,1);
        %Attribute to "satfield" the value of "Sw".
        getsatfield(j) = Sw(i);

        %Attribute the number of element to "getelemonline"
        getelemonline(j) = i;
        %Increment "j"
        j = j + 1;
    elseif 30<numcase && numcase<100
        %Define "xcolumn" and "ycolumn" according to "numcase" value
        boolean = (numcase ~= 31.5);
        xcolumn = boolean + (1 - boolean)*2;
        ycolumn = 2*boolean + (1 - boolean);
        %Choose a arbitrary position in "centelem"
        pointer = ceil(0.5*length(producelem));

        %Catch the "y" value in "centelem"
        y_value = centelem(producelem(pointer),ycolumn);
        %Initialize "pos" and "satfield"
        getxvalue(1) = 0;
        getsatfield(1) = 1 - satlimit(2);
        getelemonline(1) = 0;
        %Initialize "j"
        if (centelem(i,ycolumn) >= 0.98*y_value) && ...
                (centelem(i,ycolumn) <= 1.02*y_value)
            %Attribute to "pos" the value of "centelem" which match with
            %"y_value"
            getxvalue(jj) = centelem(i,xcolumn);
            %Attribute to "satfield" the value of "Sw".
            getsatfield(jj) = Sw(i);
            %Attribute the number of element to "getelemonline"
            getelemonline(jj) = i;
            %Increment "jj"
            jj = jj + 1;
        end  %End of IF
    elseif 300<numcase
        if numcase==330 || numcase==331
            %if (498<centelem(i,1) && centelem(i,1)<501) && centelem(i,2)<1000
            if (abs(centelem(i,1)-500)/500)<1e-2 && centelem(i,2)<1000
                %"y_value"
                getxvalue(j) = centelem(i,2);
                %Attribute to "satfield" the value of "Sw".
                getsatfield(j) = Sw(i);
                %Attribute the number of element to "getelemonline"
                getelemonline(j) = i;
                %Increment "jj"
                j = j + 1;
            end

        elseif numcase==332
            %if (abs(centelem(i,1)-250)/250)<1e-2 && centelem(i,2)<1000
            const=25; % para malha 40x40
            if ((const*29<centelem(i,1)&& centelem(i,1)<const*30) && centelem(i,2)<1000)

                %"y_value"
                getxvalue(j) = centelem(i,2);
                %Attribute to "satfield" the value of "Sw".
                getsatfield(j) = Sw(i);
                %Attribute the number of element to "getelemonline"
                getelemonline(j) = i;
                %Increment "jj"
                j = j + 1;
            end
        elseif numcase==333

            if (1.25*3<centelem(i,2) && centelem(i,2)<1.25*4) && centelem(i,1)<40
                getxvalue(j) = centelem(i,1);
                %Attribute to "satfield" the value of "Sw".
                getsatfield(j) = Sw(i);
                %Attribute the number of element to "getelemonline"
                getelemonline(j) = i;
                %Increment "jj"
                j = j + 1;
            end
        elseif numcase==342
            if (110<centelem(i,2) && centelem(i,2)<111 && centelem(i,1)<1800)
                %Attribute to "pos" the value of "centelem" which match with
                %"y_value"
                getxvalue(j) = centelem(i,1);
                %Attribute to "satfield" the value of "Sw".
                getsatfield(j) = Sw(i);
                %Attribute the number of element to "getelemonline"
                getelemonline(j) = i;
                %Increment "jj"
                j = j + 1;
            end  %End of IF
        elseif numcase==343
            if (0.4<centelem(i,2) && centelem(i,2)<0.56)
                getxvalue(j) = centelem(i,1);
                %Attribute to "satfield" the value of "Sw".
                getsatfield(j) = Sw(i);
                %Attribute the number of element to "getelemonline"
                getelemonline(j) = i;
                %Increment "jj"
                j = j + 1;

            end

        end
    end  %End of IF
end  %End of FOR


%Fix "getxvalue" and "getsatfield"
posit = sort(getxvalue);
satfield = zeros(length(getxvalue),1);
elemonline = satfield;
%Reposition the saturation field
for i = 1:length(getxvalue)
    satpointer = logical(getxvalue == posit(i));
    satfield(i) = getsatfield(satpointer);

    elemonline(i) = getelemonline(satpointer);
end  %End of FOR

%Update "elemonline" without the first positio (it is on face)
if 30<numcase && numcase<100
    elemonline = elemonline(2:length(elemonline));
else
    elemonline = elemonline(1:length(elemonline));
    posit=posit';
end

%--------------------------------------------------------------------------
%Function "getlineresultEymardProblem"
%--------------------------------------------------------------------------

%Get saturation and position to plot
function [satfieldh,satfield,posith,positd] = ...
    getlineresultEymardProblem(Sw)
%Define global parameters:
global centelem;
%Initialize "posit" and "satfield" vectors
getxvaluehorz = zeros(sqrt(size(centelem,1)),1);
getxvaluediag = getxvaluehorz;
getsatfieldhorz = getxvaluehorz;
getsatfieldiag = getxvaluehorz;
satfieldh = getxvaluehorz;
satfield = getxvaluehorz;
%Initialize tol
tol = 1/(3*length(getxvaluehorz));

%Initialize "h" and "d"
h = 1;
d = 1;
%Swept "centelem"
for i = 1:size(centelem,1)
    %Verify the horizontal path
    if abs(centelem(i,2)) <= tol
        %Attribute to "pos" the value of "centelem" which match with
        %"y_value"
        getxvaluehorz(h) = centelem(i,1);
        %Attribute to "satfield" the value of "Sw".
        getsatfieldhorz(h) = Sw(i);
        %Increment "h"
        h = h + 1;
    end  %End of IF

    %Verify the diagonal path
    if abs(centelem(i,2) - centelem(i,1)) <= tol
        %Attribute to "pos" the value of "centelem" which match with
        %"y_value"
        getxvaluediag(d) = centelem(i,1);
        %Attribute to "satfield" the value of "Sw".
        getsatfieldiag(d) = Sw(i);
        %Increment "d"
        d = d + 1;
    end  %End of IF
end  %End of FOR

%Fix "getxvalue" and "getsatfield"
posith = sort(getxvaluehorz);
positd = sort(getxvaluediag);
%Reposition the saturation field
for i = 1:length(getxvaluehorz)
    %Reorder to horizontal path
    satpointerh = logical(getxvaluehorz == posith(i));
    satfieldh(i) = getsatfieldhorz(satpointerh);
    %Reorder to diagonal path
    satpointerd = logical(getxvaluediag == positd(i));
    satfield(i) = getsatfieldiag(satpointerd);
end  %End of FOR

%Project the "x" coordinate on the diagonal axe
positd = sort(getxvaluediag)./cosd(45);

%--------------------------------------------------------------------------
%Function "storeresult"
%--------------------------------------------------------------------------

function storeresult(pos,satfield)
%Define global parameters:
global filepath resfolder;

%Create the file name
foldername = resfolder(9:length(resfolder));
filename = [foldername '.dat'];
%Create the file if the scheme is of first order
%Write the file
resfile = fopen([filepath '\' filename],'w');
%Print "pos" and "satfield" values
fprintf(resfile,'%26.16E %26.16E\r\n',[pos' satfield]');
%Close the file
fclose(resfile);

