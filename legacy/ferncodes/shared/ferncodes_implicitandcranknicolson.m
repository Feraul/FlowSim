function [M,I]=ferncodes_implicitandcranknicolson(M,I,...
                          env,dt)
global  methodhydro 
    auxnumcase=env.config.numcase;
    auxelem=env.geometry.elem;
    auxelemarea=env.geometry.elemarea;
    
    %
    if 300<auxnumcase && auxnumcase<379
        if auxnumcase==333 || auxnumcase==331 %|| numcase==347
            %para aquifero nao confinado
            coeficiente=dt^-1*SS.*auxelemarea(:);
        else
            % para quifero confinado
            coeficiente=dt^-1*MM*SS.*auxelemarea(:);
        end
        % Euler backward method
        if strcmp(methodhydro,'backward')
            % equacao 30 Qian et al 2023
            M=coeficiente.*eye(size(auxelem,1))+M;
       
            I=I+coeficiente.*eye(size(auxelem,1))*h_n;
            
        else
            % Crank-Nicolson method
            % equacao 33 Qian et al 2023
            
             I=I+(coeficiente.*eye(size(auxelem,1))-0.5*M)*h_n;
           
            M=  (coeficiente.*eye(size(auxelem,1))+0.5*M);
        end
    end

end