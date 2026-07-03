function [pressurinterp,cinterp]=ferncodes_pressureinterpHP(p,nflagface,...
    parameter,weightDMP,weightDMPc,Con,nflagfacec,dparameter)
global inedge coord bedge bcflag bcflagc numcase


pressurinterp=zeros(size(inedge,1)+size(bedge,1),1);
cinterp=zeros(size(inedge,1)+size(bedge,1),1);
% interpolação das pressões no pontos armonicos internos
pressurinterp(size(bedge,1)+1:size(bedge,1)+size(inedge,1),1)=weightDMP(:,1).*p(weightDMP(:,3))+ weightDMP(:,2).*p(weightDMP(:,4));

if 200<numcase && numcase<300
    cinterp(size(bedge,1)+1:size(bedge,1)+size(inedge,1),1)=weightDMPc(:,1).*Con(weightDMPc(:,3))+ weightDMPc(:,2).*Con(weightDMPc(:,4));
end
% interpolação das pressões no pontos medios do contorno
for ifacont=1:size(bedge,1)
    
    
    lef=bedge(ifacont,3);
    
    % Quando o ponto harmonico oposto pertece a face interior da malha ou
    % contorno de Dirichlet
    
    if bedge(ifacont,5)>200
        auxfacelef1=parameter(1,3,ifacont);
        auxfacelef2=parameter(1,4,ifacont);
        if auxfacelef1==ifacont
            faceoposto=auxfacelef2;
            atualksi=parameter(1,1,ifacont);
            opostoksi=parameter(1,2,ifacont);
        else
            faceoposto=auxfacelef1;
            atualksi=parameter(1,2,ifacont);
            opostoksi=parameter(1,1,ifacont);
        end
        normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        
        % calcula o fluxo na face "ifacelef2"
        fluxoN=normcont*bcflag(r,2);
        
        % "faceoposto" é aquele face oposto ao "ifacont"
        if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
            % faceoposto pode ser contorno de Dirichlet
            if nflagface(faceoposto,1)<200
                %  caso I. Quando a face oposto pertence ao contorno de Dirichlet
                pressure=nflagface(faceoposto,2);
                
                termo1=(opostoksi+atualksi)/atualksi;
                termo2= fluxoN/(normcont*atualksi);
                termo3= opostoksi/atualksi;
                
                pressurinterp(ifacont,1)= termo1*p(lef)-termo2-termo3*pressure;
            else
                %  caso II. Quando a face oposto pertence ao contorno de Neumann
                
                
                normcontopost=norm(coord(bedge(faceoposto,1),:)-coord(bedge(faceoposto,2),:));
                
                x=bcflag(:,1)==bedge(faceoposto,5);
                r=find(x==1);
                fluxOpost=normcontopost*bcflag(r,2);
                
                if auxfacelef1==faceoposto
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
                pressurinterp(ifacont,1)= termo1*p(lef)+termo2;
                
            end
        else
            % caso III. Quando a face oposto pertence ao interior da
            % malha
            termo1=(opostoksi+atualksi)/atualksi;
            
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
            
            pfaceopost= pesatual*p(auxelematual)+pesopost*p(auxelemopost);
            
            pressurinterp(ifacont,1)= termo1*p(auxelematual)-termo2-termo3*pfaceopost;
            
        end
        
        
    else
        pressurinterp(ifacont,1)= nflagface(ifacont,2);
    end
    
    %% ======================================================================
    
    % Quando o ponto harmonico oposto pertece a face interior da malha ou
    % contorno de Dirichlet
    if 200<numcase && numcase<300
        if bedge(ifacont,7)>200
            auxfacelef1=dparameter(1,3,ifacont);
            auxfacelef2=dparameter(1,4,ifacont);
            if auxfacelef1==ifacont
                faceoposto=auxfacelef2;
                atualksi=dparameter(1,1,ifacont);
                opostoksi=dparameter(1,2,ifacont);
            else
                faceoposto=auxfacelef1;
                atualksi=dparameter(1,2,ifacont);
                opostoksi=dparameter(1,1,ifacont);
            end
            normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
            x=bcflagc(:,1)==bedge(ifacont,7);
            r=find(x==1);
            
            % calcula o fluxo na face "ifacelef2"
            fluxoN=normcont*bcflagc(r,2);
            
            % "faceoposto" é aquele face oposto ao "ifacont"
            if faceoposto<size(bedge,1) || faceoposto==size(bedge,1)
                % faceoposto pode ser contorno de Dirichlet
                if nflagfacec(faceoposto,1)<200
                    %  caso I. Quando a face oposto pertence ao contorno de Dirichlet
                    concentr=nflagfacec(faceoposto,2);
                    
                    termo1=(opostoksi+atualksi)/atualksi;
                    termo2= fluxoN/(normcont*atualksi);
                    termo3= opostoksi/atualksi;
                    
                    cinterp(ifacont,1)= termo1*p(lef)-termo2-termo3*concentr;
                else
                    %  caso II. Quando a face oposto pertence ao contorno de Neumann
                    
                    
                    normcontopost=norm(coord(bedge(faceoposto,1),:)-coord(bedge(faceoposto,2),:));
                    
                    x=bcflagc(:,1)==bedge(faceoposto,7);
                    r=find(x==1);
                    fluxOpost=normcontopost*bcflagc(r,2);
                    
                    if auxfacelef1==faceoposto
                        atualksO=dparameter(1,1,faceoposto);
                        opostksO=dparameter(1,2,faceoposto);
                    else
                        atualksO=dparameter(1,2,faceoposto);
                        opostksO=dparameter(1,1,faceoposto);
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
                    cinterp(ifacont,1)= termo1*Con(lef)+termo2;
                    
                end
            else
                % caso III. Quando a face oposto pertence ao interior da
                % malha
                termo1=(opostoksi+atualksi)/atualksi;
                
                termo2=(fluxoN/(normcont*atualksi));
                
                % Calculo das contribuições do elemento a esquerda
                auxweightlef1=weightDMPc(faceoposto-size(bedge,1),1);
                auxweightlef2=weightDMPc(faceoposto-size(bedge,1),2);
                
                auxlef=weightDMPc(faceoposto-size(bedge,1),3);
                auxrel=weightDMPc(faceoposto-size(bedge,1),4);
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
                
                pfaceopost= pesatual*Con(auxelematual)+pesopost*Con(auxelemopost);
                
                cinterp(ifacont,1)= termo1*Con(auxelematual)-termo2-termo3*pfaceopost;
                
            end
            
            
        else
            cinterp(ifacont,1)= nflagfacec(ifacont,2);
        end
    end
end
end