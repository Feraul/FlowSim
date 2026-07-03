function [M,I,elembedge]=ferncodes_assemblematrixMPFAH(parameter,nflagface,...
    weightDMP,SS,dt,h,MM,gravrate,viscosity)
global inedge coord bedge bcflag elem numcase methodhydro ...
    keygravity dens elemarea modflowcompared Nmod normals;
%---------------------------------------------------------
% viscosity=viscosity when contaminant transport
% viscosity=mobility when two-phase flow
%---------------------------------------------------------

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
elembedge=0;
%Initialize "M" (global matrix) and "I" (known vector)
M = sparse(size(elem,1),size(elem,1)); %Prealocação de M.
I = zeros(size(elem,1),1);
jj=1;
%viscosidade ou mobilidade na face
% Montagem da matriz global
visonface = 1;
for ifacont=1:bedgesize
    % when 200<numcase<300, is groundwater transport model
    if 200<numcase && numcase<300
        % equacao de concentracao
        if numcase == 246 || numcase == 245 || numcase==247 || ...
                numcase==248 || numcase==249 ||numcase==251
            % vicosity on the boundary edge
            visonface = viscosity(ifacont,:);
        end  %End of IF
        % when numcase<200, is two-phase model
    elseif 30<numcase && numcase<200
        % equacao de saturacao "viscosity=mobility"
        visonface=sum(viscosity(ifacont,:));
    end
    % element flag
    lef=bedge(ifacont,3);
    % length of face
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    
    if bedge(ifacont,5)>200
        
        if numcase==341
           
            aaa=0.5*(coord(bedge(ifacont,1),:) + coord(bedge(ifacont,2),:));
            auxkmap = ferncodes_K(aaa(1,1),aaa(1,2));
            %----------------------------------------------------------
            %auxkmap=kmap(lef, 2);
            I(lef) = I(lef)+ normals(ifacont,2)*auxkmap(1)*nflagface(ifacont,2);
        else
            x=bcflag(:,1)==bedge(ifacont,5);
            r=find(x==1);
            % I(lef)=I(lef)- normcont*bcflag(r,2);% testes feitos em todos os
            %problemas monofásico
            I(lef)=I(lef)+ normcont*bcflag(r,2); % problema de buckley leverett Bastian
        end
    else
        % amazena os elementos vizinhos no contorno de Dirichlet e o valor
        % da condicao de contorno nela
        if strcmp(modflowcompared,'y')
            elembedge(jj,1)=bedge(ifacont,3);
            elembedge(jj,2)=nflagface(ifacont,2);
            jj=jj+1;
        end
        %-------------------------------------------------------------%
        % somando 1
        ifacelef1=parameter(1,3,ifacont);
        auxparameter1=parameter(1,1,ifacont);
        [M,I]=ferncodes_tratmentcontourlfvHP(ifacelef1,parameter,nflagface,...
            normcont,auxparameter1,visonface,lef,weightDMP,M,I);
        %-------------------------------------------------------------%
        % somando 2
        ifacelef2=parameter(1,4,ifacont);
        auxparameter2=parameter(1,2,ifacont);
        [M,I]=ferncodes_tratmentcontourlfvHP(ifacelef2,parameter,nflagface,...
            normcont,auxparameter2,visonface,lef,weightDMP,M,I);
        
        M(lef,lef)=M(lef,lef)+ visonface*normcont*(auxparameter1 + auxparameter2);
        
    end
end


for iface=1:inedgesize
    % when 200<numcase<300, is groundwater transport model
    if 200<numcase && numcase<300
        % concentration equation
        if numcase == 246 || numcase == 245 || numcase==247 ||...
                numcase==248 || numcase==249 || numcase==251
            % vicosity on the boundary edge
            visonface = viscosity(bedgesize + iface,:);
            %It is a Two-phase flow
        end  %End of IF
    elseif 30<numcase && numcase<200
        % saturation equation "viscosity=mobility"
        visonface=sum(viscosity(bedgesize + iface,:));
    end
    % element flags
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    
    %Determinação dos centróides dos elementos à direita e à esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    % calculo da norma do vetor "vd1"
    norma=norm(vd1);
    % Calculo dos fluxos parciais
    ifactual=iface+size(bedge,1);
    
    %-----------------------------------------------------------------%
    % faces de que contem os pontos de interpolação correspondente a
    % elemento a esquerda
    
    ifacelef1=parameter(1,3,ifactual);
    ifacelef2=parameter(1,4,ifactual);
    % faces de que contem os pontos de interpolação correspondente a
    % elemento a direita
    ifacerel1= parameter(2,3,ifactual);
    ifacerel2= parameter(2,4,ifactual);
    
    % Calculo das contribuições do elemento a esquerda
    mulef=weightDMP(ifactual-size(bedge,1),1);
    murel=weightDMP(ifactual-size(bedge,1),2);
    
    % calculo da contribuição, Eq. 2.12 (resp. Eq. 21) do artigo Gao and Wu 2015 (resp. Gao and Wu 2014)
    ALL=norma*murel*(parameter(1,1,ifactual)+parameter(1,2,ifactual));
    ARR=norma*mulef*(parameter(2,1,ifactual)+parameter(2,2,ifactual));
    % implementação da matriz global
    % contribuição da transmisibilidade no elemento esquerda
    M(lef,lef)=M(lef,lef)+ visonface*ALL;
    M(lef,rel)=M(lef,rel)- visonface*ARR;
    % contribuição da transmisibilidade no elemento direita
    M(rel,rel)=M(rel,rel)+ visonface*ARR;
    M(rel,lef)=M(rel,lef)- visonface*ALL;
    
    
    % contribuições do elemento a esquerda
    %------------------------ somando 1 ----------------------------------%
    termo0=visonface*norma*murel*parameter(1,1,ifactual);
    if ifacelef1<size(bedge,1) || ifacelef1==size(bedge,1)
        if nflagface(ifacelef1,1)<200
            % neste caso automaticamente introduzimos o valor dada no
            % contorno de Dirichlet
            % neste caso interpolamos pelo lpew os vertices da face
            pressureface1=nflagface(ifacelef1,2);
            I(lef)=I(lef)+  termo0*pressureface1;
            I(rel)=I(rel)-  termo0*pressureface1;
        else% quando a face é contorno de Neumann
            % agora a face em questão é "ifacelef1"
            % a ideia principal nesta rutina é isolar a pressão na face "ifacelef1"
            % e escrever eem funcção das pressões de face internas ou
            % Dirichlet ou simplesmente em funcação da pressão da VC.
            
            normcont=norm(coord(bedge(ifacelef1,1),:)-coord(bedge(ifacelef1,2),:));
            x=bcflag(:,1)==bedge(ifacelef1,5);
            r=find(x==1);
            
            % calcula o fluxo na face "ifacelef1"
            fluxoN=normcont*bcflag(r,2);
            
            % retorna as faces que formam os eixos auxiliares
            auxifacelef1=parameter(1,3,ifacelef1);
            auxifacelef2=parameter(1,4,ifacelef1);
            
            % retorna as faces os coeficientes correpondente ao face "ifacelef1"
            ksi1=parameter(1,1,ifacelef1);
            ksi2=parameter(1,2,ifacelef1);
            % identifica as faces "atual" que neste caso é "ifacelef1" e a
            % outra face oposto
            if auxifacelef1==ifacelef1
                faceoposto=auxifacelef2;
                atualksi=ksi1;
                opostoksi=ksi2;
            else
                faceoposto=auxifacelef1;
                atualksi=ksi2;
                opostoksi=ksi1;
            end
            
            
            % verifica se a face oposto a "ifacelef1" pertence ao contorno
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                if nflagface(faceoposto,1)<200
                    % caso I. Quando a face oposto pertence ao contorno da
                    % malha
                    % atribui a pressão da face de Dirichlet
                    pressure= nflagface(faceoposto,2);
                    % o coeficiente "atualksi" pentence ao face "ifacelef1"
                    termo1=((atualksi+opostoksi)/atualksi);
                    
                    termo2=(fluxoN/(normcont*atualksi))+(opostoksi/atualksi)*pressure;
                    
                    % contribuição no elemeto a esquerda
                    M(lef,lef)=M(lef,lef)- termo0*termo1;
                    I(lef)=I(lef)-termo0*termo2;
                    % contribuição no elemeto a esquerda
                    M(rel,lef)=M(rel,lef)+ termo0*termo1;
                    I(rel)=I(rel)+termo0*termo2;
                    
                    
                else
                    %Define "mobonface" (for "bedge")
                    %It is a One-phase flow. In this case, "mobility" is "1"
                    
                    normcontopost=norm(coord(bedge(faceoposto,1),:)-coord(bedge(faceoposto,2),:));
                    
                    x=bcflag(:,1)==bedge(faceoposto,5);
                    r=find(x==1);
                    fluxOpost=normcontopost*bcflag(r,2);
                    % caso II. Quando a face oposto pertence ao contorno de Neumann
                    
                    if auxifacelef1==faceoposto
                        atualksO=parameter(1,1,faceoposto);
                        opostksO=parameter(1,2,faceoposto);
                    else
                        atualksO=parameter(1,2,faceoposto);
                        opostksO=parameter(1,1,faceoposto);
                    end
                    
                    sumaksiatual=  opostksO*(atualksi + opostoksi);
                    
                    sumaksiopost= opostoksi*(atualksO + opostksO);
                    
                    
                    sumatotalnumer=  sumaksiatual - sumaksiopost;
                    
                    sumatotaldenom= atualksi*opostksO - atualksO*opostoksi;
                    
                    termo1=sumatotalnumer/sumatotaldenom;
                    
                    if abs(sumatotaldenom)<1e-5
                        termo2=0;
                    else
                        termo2=(opostoksi*(fluxOpost/(normcontopost)) - opostksO*(fluxoN/(normcont)))/sumatotaldenom;
                    end
                    M(lef,lef)=M(lef,lef)- termo0*termo1;
                    I(lef)=I(lef)- termo0*termo2;
                    
                    M(rel,lef)=M(rel,lef)+ termo0*termo1;
                    I(rel)=I(rel)+ termo0*termo2;
                    
                end
                
            else
                % caso III. Quando a face oposto pertence ao interior da
                % malha
                termo1=(ksi1+ksi2)/(atualksi);
                
                termo2=(fluxoN/(normcont*atualksi));
                
                % Calculo das contribuições do elemento a esquerda
                auxweightlef1=weightDMP(faceoposto-size(bedge,1),1);
                auxweightlef2=weightDMP(faceoposto-size(bedge,1),2);
                
                auxlef=weightDMP(faceoposto-size(bedge,1),3);
                auxrel=weightDMP(faceoposto-size(bedge,1),4);
                if auxlef==lef
                    auxelematual=auxlef;
                    auxelemopost=auxrel;
                    pesatual=auxweightlef1;
                    pesopost=auxweightlef2;
                else
                    auxelematual=auxrel;
                    pesatual=auxweightlef2;
                    pesopost=auxweightlef1;
                    auxelemopost=auxlef;
                end
                
                termo3= (opostoksi/atualksi);
                
                
                % contribuição no elemeto a direita
                M(lef,auxelematual)=M(lef,auxelematual)+ termo0*(termo3*pesatual-termo1);
                
                M(lef,auxelemopost)=M(lef,auxelemopost)+ termo0*termo3*pesopost;
                I(lef)=I(lef)-termo0*termo2;
                
                % contribuição no elemeto a esquerda
                
                M(rel,auxelematual)=M(rel,auxelematual)- termo0*(termo3*pesatual-termo1);
                
                M(rel,auxelemopost)=M(rel,auxelemopost)- termo0*termo3*pesopost;
                I(rel)=I(rel)+termo0*termo2;
                
            end
            
        end
    else
        % Calculo das contribuições do elemento a esquerda
        auxweightlef1=weightDMP(ifacelef1-size(bedge,1),1);
        auxweightlef2=weightDMP(ifacelef1-size(bedge,1),2);
        
        auxlef=weightDMP(ifacelef1-size(bedge,1),3);
        auxrel=weightDMP(ifacelef1-size(bedge,1),4);
        
        % contribuição do elemento a esquerda
        M(lef,[auxlef,auxrel])=  M(lef,[auxlef,auxrel])- termo0*[auxweightlef1,auxweightlef2];
        
        % contribuição do elemento a direita
        M(rel,[auxlef,auxrel])=  M(rel,[auxlef,auxrel])+ termo0*[auxweightlef1,auxweightlef2];
        
    end
    
    %--------------------------- somando 2 -------------------------------%
    termo0=visonface*norma*murel*parameter(1,2,ifactual);
    if ifacelef2<size(bedge,1) || ifacelef2==size(bedge,1)
        if nflagface(ifacelef2,1)<200
            % neste caso automaticamente introduzimos o valor dada no
            % contorno de Dirichlet
            % neste caso interpolamos pelo lpew os vertices da face
            pressureface2=nflagface(ifacelef2,2);
            I(lef)=I(lef)+ termo0*pressureface2;
            I(rel)=I(rel)- termo0*pressureface2;
        else% quando a face é contorno de Neumann
            % agora a face em questão é "ifacelef2"
            % a ideia principal nesta rutina é isolar a pressão na face "ifacelef2"
            % e escrever eem funcção das pressões de face internas ou
            % Dirichlet ou simplesmente em funcação da pressão da VC.
            
            normcont=norm(coord(bedge(ifacelef2,1),:)-coord(bedge(ifacelef2,2),:));
            x=bcflag(:,1)==bedge(ifacelef2,5);
            r=find(x==1);
            
            % calcula o fluxo na face "ifacelef2"
            fluxoN=normcont*bcflag(r,2);
            
            % retorna as faces que formam os eixos auxiliares
            auxifacelef1=parameter(1,3,ifacelef2);
            auxifacelef2=parameter(1,4,ifacelef2);
            
            % retorna as faces os coeficientes correpondente ao face "ifacelef2"
            ksi1=parameter(1,1,ifacelef2);
            ksi2=parameter(1,2,ifacelef2);
            % identifica as faces "atual" que neste caso é "ifacelef2" e a
            % outra face oposto
            if auxifacelef1==ifacelef2
                faceoposto=auxifacelef2;
                atualksi=ksi1;
                opostoksi=ksi2;
            else
                faceoposto=auxifacelef1;
                atualksi=ksi2;
                opostoksi=ksi1;
            end
            
            % verifica se a face oposto a "ifacelef1" pertence ao contorno
            
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                if nflagface(faceoposto,1)<200
                    % caso I. Quando a face oposto pertence ao contorno da
                    % malha
                    % atribui a pressão da face de Dirichlet
                    pressure= nflagface(faceoposto,2);
                    % o coeficiente "atualksi" pentence ao face "ifacelef1"
                    termo1=((atualksi+opostoksi)/atualksi);
                    
                    termo2=(fluxoN/(normcont*atualksi))+(opostoksi/atualksi)*pressure;
                    
                    % contribuição no elemeto a esquerda
                    M(lef,lef)=M(lef,lef)- termo0*termo1;
                    I(lef)=I(lef)-termo0*termo2;
                    % conntribuição no elemento a direita
                    M(rel,lef)=M(rel,lef)+ termo0*termo1;
                    I(rel)=I(rel)+termo0*termo2;
                    
                else
                    
                    normcontopost=norm(coord(bedge(faceoposto,1),:)-coord(bedge(faceoposto,2),:));
                    
                    x=bcflag(:,1)==bedge(faceoposto,5);
                    r=find(x==1);
                    fluxOpost=normcontopost*bcflag(r,2);
                    % caso II. Quando a face oposto pertence ao contorno de Neumann
                    
                    if auxifacelef1==faceoposto
                        atualksO=parameter(1,1,faceoposto);
                        opostksO=parameter(1,2,faceoposto);
                    else
                        atualksO=parameter(1,2,faceoposto);
                        opostksO=parameter(1,1,faceoposto);
                    end
                    
                    sumaksiatual=  opostksO*(atualksi + opostoksi);
                    
                    sumaksiopost= opostoksi*(atualksO + opostksO);
                    
                    
                    sumatotalnumer=  sumaksiatual - sumaksiopost;
                    
                    sumatotaldenom= atualksi*opostksO - atualksO*opostoksi;
                    
                    termo1=sumatotalnumer/sumatotaldenom;
                    
                    if abs(sumatotaldenom)<1e-5
                        termo2=0;
                    else
                        termo2=(opostoksi*(fluxOpost/(normcontopost)) - opostksO*(fluxoN/(normcont)))/sumatotaldenom;
                    end
                    M(lef,lef)=M(lef,lef)- termo0*termo1;
                    I(lef)=I(lef)+ termo0*termo2;
                    
                    M(rel,lef)=M(rel,lef)+ termo0*termo1;
                    I(lef)=I(lef)- termo0*termo2;
                    
                end
                
            else
                % caso III. Quando a face oposto pertence ao interior da
                % malha
                termo1=(ksi1+ksi2)/(atualksi);
                
                termo2=(fluxoN/(normcont*atualksi));
                
                % Calculo das contribuições do elemento a esquerda
                auxweightlef1=weightDMP(faceoposto-size(bedge,1),1);
                auxweightlef2=weightDMP(faceoposto-size(bedge,1),2);
                
                auxlef=weightDMP(faceoposto-size(bedge,1),3);
                auxrel=weightDMP(faceoposto-size(bedge,1),4);
                if auxlef==lef
                    auxelematual=auxlef;
                    auxelemopost=auxrel;
                    pesatual=auxweightlef1;
                    pesopost=auxweightlef2;
                else
                    auxelematual=auxrel;
                    pesatual=auxweightlef2;
                    pesopost=auxweightlef1;
                    auxelemopost=auxlef;
                end
                
                termo3= (opostoksi/atualksi);
                
                
                % contribuição no elemeto a direita
                M(lef,auxelematual)=M(lef,auxelematual)+ termo0*(termo3*pesatual-termo1);
                
                M(lef,auxelemopost)=M(lef,auxelemopost)+ termo0*termo3*pesopost;
                I(lef)=I(lef)-termo0*termo2;
                
                % contribuição no elemeto a esquerda
                
                M(rel,auxelematual)=M(rel,auxelematual)- termo0*(termo3*pesatual-termo1);
                
                M(rel,auxelemopost)=M(rel,auxelemopost)- termo0*termo3*pesopost;
                I(rel)=I(rel)+termo0*termo2;
                
            end
            
        end
    else
        % Calculo das contribuições do elemento a esquerda
        auxweightlef1=weightDMP(ifacelef2-size(bedge,1),1);
        auxweightlef2=weightDMP(ifacelef2-size(bedge,1),2);
        
        auxlef=weightDMP(ifacelef2-size(bedge,1),3);
        auxrel=weightDMP(ifacelef2-size(bedge,1),4);
        
        % contribuição do elemento a esquerda
        M(lef,[auxlef, auxrel])=  M(lef,[auxlef,auxrel])- termo0*[auxweightlef1,auxweightlef2];
        
        % contribuição do elemento a direita
        M(rel,[auxlef,auxrel])=  M(rel,[auxlef,auxrel])+ termo0*[auxweightlef1,auxweightlef2];
        
    end
    %%
    % contribuições do elemento a direita
    % somando 1
    termo0=visonface*norma*mulef*parameter(2,1,ifactual);
    if ifacerel1<size(bedge,1) || ifacerel1==size(bedge,1)
        if nflagface(ifacerel1,1)<200
            % neste caso automaticamente introduzimos o valor dada no
            % contorno de Dirichlet
            pressurefacerel1=nflagface(ifacerel1,2);
            I(rel)=I(rel) + termo0*pressurefacerel1;
            I(lef)=I(lef) - termo0*pressurefacerel1;
        else % quando a face é contorno de Neumann
            % agora a face em questão é "ifacerel1"
            % a ideia principal nesta rutina é isolar a pressão na face "ifacerel1"
            % e escrever eem funcção das pressões de face internas ou
            % Dirichlet ou simplesmente em funcação da pressão da VC.
            normcont=norm(coord(bedge(ifacerel1,1),:)-coord(bedge(ifacerel1,2),:));
            x=bcflag(:,1)==bedge(ifacerel1,5);
            r=find(x==1);
            fluxoN=normcont*bcflag(r,2);
            % lembre-se que a face "ifacerel1" pertence ao contorno, então
            % as faces que definem os eixos auxiliares são
            auxifacerel1=parameter(1,3,ifacerel1);
            auxifacerel2=parameter(1,4,ifacerel1);
            % um deles é mesmo "ifacerel1" e o outro e oposto.
            % as coeficientes que representam a cada face
            ksi1=parameter(1,1,ifacerel1);
            ksi2=parameter(1,2,ifacerel1);
            % aqui calculamos quem é face atual "ifacerel1" e o face oposto
            % e também os coeficientes adequados
            if auxifacerel1==ifacerel1
                faceoposto=auxifacerel2;
                atualksi=ksi1;
                opostoksi=ksi2;
            else
                faceoposto=auxifacerel1;
                atualksi=ksi2;
                opostoksi=ksi1;
            end
            % vejamos a qual face pertence a face oposto
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                if nflagface(faceoposto,1)<200
                    % caso I. Quando a face oposto pertence ao contorno da
                    % malha atribui a pressão da face de Dirichlet
                    pressure= nflagface(faceoposto,2);
                    
                    termo1=(atualksi+opostoksi)/atualksi;
                    
                    termo2=(fluxoN/(normcont*atualksi))+(opostoksi/atualksi)*pressure;
                    
                    % contribuição no elemeto a esquerda
                    M(rel,rel)=M(rel,rel)- termo0*termo1;
                    I(rel)=I(rel)-termo0*termo2;
                    
                    M(lef,rel)=M(lef,rel)+ termo0*termo1;
                    I(lef)=I(lef)+termo0*termo2;
                else
                    
                    normcontopost=norm(coord(bedge(faceoposto,1),:)-coord(bedge(faceoposto,2),:));
                    
                    x=bcflag(:,1)==bedge(faceoposto,5);
                    r=find(x==1);
                    fluxOpost=normcontopost*bcflag(r,2);
                    % caso II. Quando a face oposto pertence ao contorno de Neumann
                    ksi1a=parameter(1,1,faceoposto);
                    ksi2b=parameter(1,2,faceoposto);
                    
                    if auxifacelef1==faceoposto
                        atualksO=ksi1a;
                        opostksO=ksi2b;
                    else
                        atualksO=ksi2b;
                        opostksO=ksi1a;
                    end
                    
                    sumaksiatual=  opostksO*(atualksi + opostoksi);
                    sumaksiopost= opostoksi*(atualksO + opostksO);
                    sumatotalnumer=  sumaksiatual - sumaksiopost;
                    sumatotaldenom= atualksi*opostksO - atualksO*opostoksi;
                    
                    termo1=sumatotalnumer/sumatotaldenom;
                    
                    if abs(sumatotaldenom)<1e-5
                        termo2=0;
                    else
                        termo2=(opostoksi*(fluxOpost/(normcontopost)) - opostksO*(fluxoN/(normcont)))/sumatotaldenom;
                        
                    end
                    M(rel,rel)=M(rel,rel)- termo0*termo1;
                    I(rel)=I(rel)+ termo0*termo2;
                    
                    M(lef,rel)=M(lef,rel)+ termo0*termo1;
                    I(lef)=I(lef)- termo0*termo2;
                end
            else
                % caso III. Quando a face oposto pertence à interior da
                % malha
                termo1=(opostoksi+atualksi)/(atualksi);
                
                termo2=(fluxoN/(normcont*atualksi));
                
                termo3= (opostoksi/atualksi);
                
                % lembre-se que a face "ifacerel1" pertence ao contorno.
                % a face oposto a "ifacerel1" é chamado "faceoposto" e
                % pertence ao face interior, então a pressão nessa deve
                % ser interpolado para isso calculamos os pesos e
                % elementos que facem parte da face "oposto"
                
                auxweight1=weightDMP(faceoposto-size(bedge,1),1);
                auxweight2=weightDMP(faceoposto-size(bedge,1),2);
                
                auxlef=weightDMP(faceoposto-size(bedge,1),3);
                auxrel=weightDMP(faceoposto-size(bedge,1),4);
                
                % lembrando que um "auxlef" ou "auxrel" é "rel"
                if auxlef==rel
                    auxelematual=auxlef;
                    pesatual=auxweight1;
                    
                    auxelemopost=auxrel;
                    pesopost=auxweight2;
                else
                    auxelematual=auxrel;
                    pesatual=auxweight2;
                    
                    pesopost=auxweight1;
                    auxelemopost=auxlef;
                end
                
                % contribuição no elemeto a direita
                M(rel,auxelematual)=M(rel,auxelematual)+ termo0*(termo3*pesatual-termo1 );
                
                M(rel,auxelemopost)=M(rel,auxelemopost)+ termo0*termo3*pesopost;
                
                I(rel)=I(rel)-termo0*termo2;
                
                % contribuição no elemeto a esquerda
                
                M(lef,auxelematual)=M(lef,auxelematual)- termo0*(termo3*pesatual-termo1 );
                
                M(lef,auxelemopost)=M(lef,auxelemopost)- termo0*termo3*pesopost;
                I(lef)=I(lef)+termo0*termo2;
                
                
            end
            
        end
    else
        % Calculo das contribuições do elemento a esquerda
        auxweightrel1=weightDMP(ifacerel1-size(bedge,1),1);
        auxweightrel2=weightDMP(ifacerel1-size(bedge,1),2);
        
        auxlef=weightDMP(ifacerel1-size(bedge,1),3);
        auxrel=weightDMP(ifacerel1-size(bedge,1),4);
        % contribuição do elemento a esquerda
        M(lef,[auxlef,auxrel])=  M(lef,[auxlef,auxrel])+ termo0*[auxweightrel1,auxweightrel2];
        % contribuição do elemento a direita
        M(rel,[auxlef,auxrel])=  M(rel,[auxlef,auxrel])- termo0*[auxweightrel1,auxweightrel2];
    end
    % somando 2
    termo0=visonface*norma*mulef*parameter(2,2,ifactual);
    if ifacerel2<size(bedge,1) || ifacerel2==size(bedge,1)
        if nflagface(ifacerel2,1)<200
            % neste caso automaticamente introduzimos o valor dada no
            % contorno de Dirichlet
            % neste caso interpolamos pelo lpew os vertices da face
            
            pressurefacerel2=nflagface(ifacerel2,2);
            I(rel)=I(rel) + termo0*pressurefacerel2;
            I(lef)=I(lef) - termo0*pressurefacerel2;
        else% quando a face é contorno de Neumann
            % agora a face em questão é "ifacerel1"
            % a ideia principal nesta rutina é isolar a pressão na face "ifacerel2"
            % e escrever eem funcção das pressões de face internas ou
            % Dirichlet ou simplesmente em funcação da pressão da VC.
            
            
            normcont=norm(coord(bedge(ifacerel2,1),:)-coord(bedge(ifacerel2,2),:));
            x=bcflag(:,1)==bedge(ifacerel2,5);
            r=find(x==1);
            fluxoN=normcont*bcflag(r,2);
            % lembre-se que a face "ifacerel1" pertence ao contorno, então
            % as faces que definem os eixos auxiliares são
            auxifacerel1=parameter(1,3,ifacerel2);
            auxifacerel2=parameter(1,4,ifacerel2);
            % um deles é mesmo "ifacerel2" e o outro e oposto.
            % as coeficientes que representam a cada face
            ksi1=parameter(1,1,ifacerel2);
            ksi2=parameter(1,2,ifacerel2);
            % aqui calculamos quem é face atual "ifacerel2" e o face oposto
            % e também os coeficientes adequados
            if auxifacerel1==ifacerel2
                faceoposto=auxifacerel2;
                atualksi=ksi1;
                opostoksi=ksi2;
            else
                faceoposto=auxifacerel1;
                atualksi=ksi2;
                opostoksi=ksi1;
            end
            
            
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                if nflagface(faceoposto,1)<200
                    % caso I. Quando a face oposto pertence ao contorno da
                    % malha atribui a pressão da face de Dirichlet
                    pressure= nflagface(faceoposto,2);
                    
                    termo1=(atualksi+opostoksi)/atualksi;
                    
                    termo2=(fluxoN/(normcont*atualksi))+(opostoksi/atualksi)*pressure;
                    
                    % contribuição no elemeto a esquerda
                    M(rel,rel)=M(rel,rel)- termo0*termo1;
                    I(rel)=I(rel)-termo0*termo2;
                    
                    M(lef,rel)=M(lef,rel)+ termo0*termo1;
                    I(lef)=I(lef)+termo0*termo2;
                    
                else
                    
                    normcontopost=norm(coord(bedge(faceoposto,1),:)-coord(bedge(faceoposto,2),:));
                    
                    x=bcflag(:,1)==bedge(faceoposto,5);
                    r=find(x==1);
                    fluxOpost=normcontopost*bcflag(r,2);
                    % caso II. Quando a face oposto pertence ao contorno de Neumann
                    ksi1a=parameter(1,1,faceoposto);
                    ksi2b=parameter(1,2,faceoposto);
                    
                    if auxifacelef1==faceoposto
                        atualksO=ksi1a;
                        opostksO=ksi2b;
                    else
                        atualksO=ksi2b;
                        opostksO=ksi1a;
                    end
                    
                    sumaksiatual=  opostksO*(atualksi + opostoksi);
                    
                    sumaksiopost= opostoksi*(atualksO + opostksO);
                    
                    
                    sumatotalnumer=  sumaksiatual - sumaksiopost;
                    
                    sumatotaldenom= atualksi*opostksO - atualksO*opostoksi;
                    
                    termo1=sumatotalnumer/sumatotaldenom;
                    
                    if abs(sumatotaldenom)<1e-5
                        termo2=0;
                    else
                        termo2=(opostoksi*(fluxOpost/(normcontopost)) - opostksO*(fluxoN/(normcont)))/sumatotaldenom;
                    end
                    M(rel,rel)=M(rel,rel)- termo0*termo1;
                    I(rel)=I(rel)+ termo0*termo2;
                    
                    M(lef,rel)=M(lef,rel)+ termo0*termo1;
                    I(lef)=I(lef)- termo0*termo2;
                    
                end
                
            else
                % caso III. Quando a face oposto pertence à interior da
                % malha
                termo1=(opostoksi+atualksi)/(atualksi);
                
                termo2=(fluxoN/(normcont*atualksi));
                
                termo3= (opostoksi/atualksi);
                
                % lembre-se que a face "ifacerel2" pertence ao contorno.
                % a face oposto a "ifacerel2" é chamado "faceoposto" e
                % pertence ao face interior, então a pressão nessa deve
                % ser interpolado para isso calculamos os pesos e
                % elementos que facem parte da face "oposto"
                
                auxweight1=weightDMP(faceoposto-size(bedge,1),1);
                auxweight2=weightDMP(faceoposto-size(bedge,1),2);
                
                auxlef=weightDMP(faceoposto-size(bedge,1),3);
                auxrel=weightDMP(faceoposto-size(bedge,1),4);
                
                % lembrando que um "auxlef" ou "auxrel" é "rel"
                if auxlef==rel
                    auxelematual=auxlef;
                    pesatual=auxweight1;
                    
                    auxelemopost=auxrel;
                    pesopost=auxweight2;
                else
                    auxelematual=auxrel;
                    pesatual=auxweight2;
                    
                    pesopost=auxweight1;
                    auxelemopost=auxlef;
                end
                
                % contribuição no elemeto a direita
                M(rel,auxelematual)=M(rel,auxelematual)+ termo0*(termo3*pesatual-termo1 );
                
                M(rel,auxelemopost)=M(rel,auxelemopost)+ termo0*termo3*pesopost;
                
                I(rel)=I(rel)-termo0*termo2;
                
                % contribuição no elemeto a esquerda
                
                M(lef,auxelematual)=M(lef,auxelematual)- termo0*(termo3*pesatual-termo1 );
                
                M(lef,auxelemopost)=M(lef,auxelemopost)- termo0*termo3*pesopost;
                I(lef)=I(lef)+termo0*termo2;
                
            end
            
        end
        
    else
        % Calculo das contribuições do elemento a esquerda
        auxweightrel1=weightDMP(ifacerel2-size(bedge,1),1);
        auxweightrel2=weightDMP(ifacerel2-size(bedge,1),2);
        
        auxlef=weightDMP(ifacerel2-size(bedge,1),3);
        auxrel=weightDMP(ifacerel2-size(bedge,1),4);
        
        % contribuição do elemento a direita
        M(lef,[auxlef,auxrel])=  M(lef,[auxlef,auxrel])+ termo0*[auxweightrel1,auxweightrel2];
        
        % contribuição do elemento a esquerda
        M(rel,[auxlef,auxrel])=  M(rel,[auxlef,auxrel])- termo0*[auxweightrel1,auxweightrel2];
        
    end
    % contribuicao do termo gravitacional
    if strcmp(keygravity,'y')
        if numcase<200
            % escoamento bifasico oleo-agua
            averagedensity=(viscosity(bedgesize + iface,:)*dens')/visonface;
            m=averagedensity*gravrate(bedgesize + iface,1);
        else
            % concentracao soluto-solvente
            m=dens(1,1)*gravrate(bedgesize + iface,1)/visonface;
        end
        I(inedge(iface,3))=I(inedge(iface,3))+visonface*m;
        I(inedge(iface,4))=I(inedge(iface,4))-visonface*m;
    end
end
% calcula um problema transiente
[M,I]=ferncodes_implicitandcranknicolson(M,I,SS,dt,MM,h);
%% utilize somente quando o teste vai ser comparado com resultados do
% modflow
if strcmp(modflowcompared,'y')
    for iw = 1:size(elembedge,1)
        M(elembedge(iw,1),:)=0*M(elembedge(iw,1),:);
        M(elembedge(iw,1),elembedge(iw,1))=1;
        I(elembedge(iw,1))=elembedge(iw,2);
    end
end
end