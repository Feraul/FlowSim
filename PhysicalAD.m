function  [Sleft,Sright,mLLF]=PhysicalAD(Sw,taylorterms,limiterflag,flagknownvert,satonvertices,...
    constraint,mlplimiter,leftelem,rightelem,vertices,...
    verticescoord,Sleftlim,Sleftnonlim,Srightlim,Srightnonlim,dotvn,tol,rijL,rijR)
global centelem numcase
mLLF=1;
%%
if numcase==103
    intermin=-1;
    intermax=1;
elseif numcase>200
    intermin=10;
    intermax=0;
else
    intermin=0;
    intermax=1;
end
    
if strcmp(limiterflag{9},'on') % Estrategia MOOD classico
    Samon=Sleft*(dotvn>0)+Sright*(1-(dotvn>0));
    Sajus=Sleft*(dotvn<0)+Sright*(1-(dotvn<0));
    SmonUP=Sw(leftelem)*(dotvn>0)+Sw(rightelem)*(1-(dotvn>0));
    SjusUP=Sw(leftelem)*(dotvn<0)+Sw(rightelem)*(1-(dotvn<0));
    
    if SjusUP<= Sajus && Sajus<=Samon && Samon<= SmonUP
        
        %Sleft=Sleftlimit;
        %Sright=Srightlimit;
        %%
        %==================================================================
        %Define the elements that share the edge evaluated
        elemeval = [leftelem rightelem];
        %Get the saturation value recovered on each quadrature point ("on_q")
        %Get the statment regarding to mlp limiter:
        boolmlp = (length(mlplimiter) == 1);
        mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
        limiterflagaux=limiterflag;
        limiterflagaux{1}='on';
        limiterflagaux{4}='wf';
        %limiterflagaux{7}=constaltaordem; % ja tentamos com -1 mas ele oscila quando tratamos o problema de buckley-leverett
        % sobre malha quadrilateral ortogonal.
        faceorderauxlef=2;
        Sleftlimit = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,...
            limiterflagaux,faceorderauxlef,constraint,flagknownvert,satonvertices,...
            mlpbyelem,centelem(elemeval,1:2));
        %======================================================================
        % Define the elements that share the edge evaluated
        elemeval = [rightelem leftelem];
        % Get the saturation value recovered on each quadrature point ("on_q")
        %     Get the statment regarding to mlp limiter:
        boolmlp = (length(mlplimiter) == 1);
        %     Get the statment regarding to mlp limiter:
        mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
        limiterflagaux=limiterflag;
        limiterflagaux{1}='on';
        limiterflagaux{4}='wf';
        %limiterflagaux{7}=constaltaordem;% ja tentamos com -1 mas ele oscila quando tratamos o problema de buckley-leverett
        % sobre malha quadrilateral ortogonal.
        faceorderauxright=2;
        Srightlimit = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,...
            limiterflagaux,faceorderauxright,constraint,flagknownvert,satonvertices,...
            mlpbyelem,centelem(elemeval,1:2));
        %==================================================================
    else
        
        Sleft=0.5*(Sleft+Sw(leftelem));
        Sright=0.5*(Sright+Sw(rightelem));
    end
end

%==================================================================
if strcmp(limiterflag{11},'on') % MUSCL Modificado
    theta=acos(dot(rijL,rijR)/(norm(rijL)*norm(rijR)));   
    SmonUPLIM= Sleftlim*(dotvn>0)+Srightlim*(1-(dotvn>0));
    SjusUPLIM= Sleftlim*(dotvn<0)+Srightlim*(1-(dotvn<0));
    
    if intermin<Sleftnonlim && Sleftnonlim<intermax   % problema reservatorio
    %if -1<Sleftnonlim && Sleftnonlim<1 % problema senoidal
        %==================================================================
        if min([SmonUPLIM,SjusUPLIM])> Sleftnonlim | Sleftnonlim>max([SmonUPLIM,SjusUPLIM])  %viola principio de maximo local
            %Sleft=Sleftlimit;
            Sleft=0.5*(SmonUPLIM+SjusUPLIM);
           
        else
           
            [Sat_max, Sat_min]=ferncodes_Saturation_max_min(leftelem,Sw);
            phi_left=max(0,(Sat_max-Sleftnonlim)/(Sat_max-Sat_min)); % pior resultado
            
            A=min(1,phi_left);
            fiGleft=(1-cos(A*pi))/2;
            %if theta<1e-10
            %    fiGleft=fiGleft;
            %else
            %    fiGleft=min(fiGleft,theta);
            %end
            Sleft=fiGleft*Sleftnonlim+(1-fiGleft)*Sleftlim;             
            
        end
    else
        Sleft=0.5*(SmonUPLIM+SjusUPLIM);
    end
    if intermin<Srightnonlim && Srightnonlim<intermax   % problema reservatorio
    %if -1<Srightnonlim && Srightnonlim<1 % problema senoidal
        
        if min([SmonUPLIM,SjusUPLIM])> Srightnonlim | Srightnonlim>max([SmonUPLIM,SjusUPLIM]) % viola principio de maximo local
            %Sright=Srightlimit;
            Sright=0.5*(SmonUPLIM+SjusUPLIM);
        else
            
            [Sat_max, Sat_min]=ferncodes_Saturation_max_min(rightelem,Sw);
            phi_right = max(0,(Sat_max-Srightnonlim)/(Sat_max-Sat_min)); % pior resultado
            
            B=min(1,phi_right);
            fiGright=(1-cos(B*pi))/2;
            %if theta<1e-10
            %   fiGright= fiGright;
            %else
            %   fiGright=min(fiGright,theta);
            %end
            Sright=fiGright*Srightnonlim+(1-fiGright)*Srightlim;
            
        end
    else
        Sright=0.5*(SmonUPLIM+SjusUPLIM);
    end
end
%==================================================================
if strcmp(limiterflag{12},'on') % Estrategia MOOD-Japones
    
    
    %Samon=Sleft*(dotvn>0)+Sright*(1-(dotvn>0));
    %Sajus=Sleft*(dotvn<0)+Sright*(1-(dotvn<0));
    
    SamonR=Sleft*(dotvn>0)+Sright*(1-(dotvn>0));
    SajusR=Sleft*(dotvn<0)+Sright*(1-(dotvn<0));
    
    Smon=Sw(leftelem)*(dotvn>0)+Sw(rightelem)*(1-(dotvn>0));
    Sjus=Sw(leftelem)*(dotvn<0)+Sw(rightelem)*(1-(dotvn<0));
    
    Smid = (SamonR+SajusR)/2;
    %Discrete:
    %"fw" has three values: fw(Sleft) is fw(1), fw(Sright) is fw(3)
    [fw,~,gama,] = twophasevar([SamonR Smid SajusR],numcase);
    % fw(1) ---> fluxo fracional no elemento fw(Sleft)
    % fw(2) ---> fluxo fracional no elemento fw(Smid)
    % fw(3) ---> fluxo fracional no elemento fw(Sright)
    %Calculate the Rankine-Hugoniot ratio:
    % veja o tese Darlan pag. 165
    [dfwdS_rh,dgamadS_rh] = ...
        calcdfunctiondS([fw(1) fw(3)],[gama(1) gama(3)],[SamonR SajusR],0);
    % VL: velocidade montante
    % VR: velocidade jusante
    % VF: velocidade da frente
    Vmon=(fw(2) - fw(1))/(Smid - SamonR); % velocidade montante
    Vjusan=(fw(2) - fw(3))/(Smid - SajusR); % Velocidade jusante
    VF=dfwdS_rh;
    
    if  (Vjusan<=VF && VF<=Vmon) && (Smid~=0) && abs(dotvn)>tol
        %if  Smon< SamonR || SamonR<SajusR || SajusR<Sjus
        if SamonR<SajusR % Não violo Entropia, mas violo DMP local
            % Portanto, utilizamos um argumento limitado
            %Sleft=Sleftlimit;
            %Sright=Srightlimit;
            %==================================================
            
            % Portanto, utilizamos um argumento limitado
            [Sat_max, Sat_min]=Saturation_max_min(leftelem,Sw);
            
            %phi_left=max(0,(Sat_max-(max([Smon,SamonR,Sjus])/min([Smon,SamonR,Sjus])))/(Sat_max-Sat_min));
            %terceiro melhor resultado
            phi_left=max(0,(Sat_max-(Sleft/(Sright+1e-16)))/(Sat_max-Sat_min));
            %segundo melhor resultado
            %phi_left=max(0,(Sat_max-Sleft)/(Sat_max-Sat_min)); % pior resultado
            %phi_left=max(0,(1-(Sleft/Sright))); % primeiro lugar em resultado
            A=min(1,phi_left);
            fiGleft=(1-cos(A*pi))/2;
            %fiGleft=min(fiGleft,fiFace);
            SleftRR=fiGleft*Sleft+(1-fiGleft)*Sleftlimit;
            
            [Sat_max, Sat_min]=Saturation_max_min(rightelem,Sw);
            
            %phi_right =max(0,(Sat_max-(max([Smon,SajusR,Sjus])/min([Smon,SajusR,Sjus])))/(Sat_max-Sat_min));%
            %terceiro melhor resultado
            phi_right =max(0,(Sat_max-(Sright/(Sleft+1e-16)))/(Sat_max-Sat_min));
            %segundo melhor resultado
            %phi_right = max(0,(Sat_max-Sright)/(Sat_max-Sat_min)); % pior resultado
            %phi_right = max(0,1-(Sright/Sleft)); % primeiro lugar em resultado
            B=min(1,phi_right);
            fiGright=(1-cos(B*pi))/2;
            
            %fiGright=min(fiGright,fiFace);
            
            SrightRR=fiGright*Sright+(1-fiGright)*Srightlimit;
            
            SamonR=SleftRR*(dotvn>0)+SrightRR*(1-(dotvn>0));
            SajusR=SleftRR*(dotvn<0)+SrightRR*(1-(dotvn<0));
            
            if Smon< SamonR || SamonR<SajusR || SajusR<Sjus
                
                SleftRR=Sleftlimit;
                SrightRR=Srightlimit;
                % moficacion 03_05_2019
                SamonR=SleftRR*(dotvn>0)+SrightRR*(1-(dotvn>0));
                SajusR=SleftRR*(dotvn<0)+SrightRR*(1-(dotvn<0));
                
                if Smon< SamonR || SamonR<SajusR || SajusR<Sjus
                    SleftRR=Sw(leftelem);
                    SrightRR=Sw(rightelem);
                    
                end
            end
            Sleft=SleftRR;
            Sright=SrightRR;
            %=====================================================
        else
            mLLF=0;
        end
        %==========================================================
        
    elseif (Vmon<VF || VF<Vjusan) && (Smid~=0) && abs(dotvn)>tol
        % quando viola a entropia utilizamos um argumento limitado.
        %if  Smon< SamonR || SamonR<SajusR || SajusR<Sjus
        if SamonR<SajusR
            %Sleft=Sleftlimit;
            %Sright=Srightlimit;
            %==================================================
            
            % Portanto, utilizamos um argumento limitado
            [Sat_max, Sat_min]=Saturation_max_min(leftelem,Sw);
            
            %phi_left=max(0,(Sat_max-(max([Smon,SamonR,Sjus])/min([Smon,SamonR,Sjus])))/(Sat_max-Sat_min));
            %terceiro melhor resultado
            phi_left=max(0,(Sat_max-(Sleft/(Sright+1e-16)))/(Sat_max-Sat_min));
            %segundo melhor resultado
            %phi_left=max(0,(Sat_max-Sleft)/(Sat_max-Sat_min)); % pior resultado
            %phi_left=max(0,(1-(Sleft/Sright))); % primeiro lugar em resultado
            A=min(1,phi_left);
            fiGleft=(1-cos(A*pi))/2;
            %fiGleft=min(fiGleft,fiFace);
            SleftRR=fiGleft*Sleft+(1-fiGleft)*Sleftlimit;
            
            [Sat_max, Sat_min]=Saturation_max_min(rightelem,Sw);
            
            %phi_right =max(0,(Sat_max-(max([Smon,SajusR,Sjus])/min([Smon,SajusR,Sjus])))/(Sat_max-Sat_min));%
            %terceiro melhor resultado
            phi_right =max(0,(Sat_max-(Sright/(Sleft+1e-16)))/(Sat_max-Sat_min));
            %segundo melhor resultado
            %phi_right = max(0,(Sat_max-Sright)/(Sat_max-Sat_min)); % pior resultado
            %phi_right = max(0,1-(Sright/Sleft)); % primeiro lugar em resultado
            B=min(1,phi_right);
            fiGright=(1-cos(B*pi))/2;
            
            %fiGright=min(fiGright,fiFace);
            
            SrightRR=fiGright*Sright+(1-fiGright)*Srightlimit;
            
            SamonR=SleftRR*(dotvn>0)+SrightRR*(1-(dotvn>0));
            SajusR=SleftRR*(dotvn<0)+SrightRR*(1-(dotvn<0));
            
            if Smon< SamonR || SamonR<SajusR || SajusR<Sjus
                
                SleftRR=Sleftlimit;
                SrightRR=Srightlimit;
                % moficacion 03_05_2019
                SamonR=SleftRR*(dotvn>0)+SrightRR*(1-(dotvn>0));
                SajusR=SleftRR*(dotvn<0)+SrightRR*(1-(dotvn<0));
                
                if Smon< SamonR || SamonR<SajusR || SajusR<Sjus
                    SleftRR=Sw(leftelem);
                    SrightRR=Sw(rightelem);
                    
                end
            end
            Sleft=SleftRR;
            Sright=SrightRR;
            %==================================================
            
        else
            Sleft=Sleftlimit;
            Sright=Srightlimit;
        end
    end
end
% PAD

end