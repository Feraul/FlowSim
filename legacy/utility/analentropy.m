function [Sleft, Sright]=analentropy(Sright,Sleft,numcase,dotvn,leftelem,rightelem,Sw,Srightaux,Sleftaux)

%%
if sign(dotvn)<0 && abs(dotvn)>=1e-10
    if Sw(rightelem)<Sright || Sright< Sleft || Sleft<Sw(leftelem)
        
        Sleft= Sleftaux;
        if Sw(rightelem)<Sright ||Sright< Sleft
            Sleft= Sw(leftelem);
            if Sw(rightelem)<Sright ||Sright< Sleft
                Sright=Srightaux;
                if Sw(rightelem)<Sright ||Sright< Sleft
                    Sright= Sw(rightelem);
                end
            end
        end
        
    end
elseif abs(dotvn)<1e-10
    Sleft=Sw(leftelem);
    Sright=Sw(rightelem);
else
    if Sw(leftelem)<Sleft || Sleft<Sright || Sright<Sw(rightelem)
        Sright=Srightaux;
        if Sw(leftelem)<Sleft ||Sleft<Sright
            Sright=Sw(rightelem);
            if Sw(leftelem)<Sleft ||Sleft<Sright
                Sleft=Sleftaux;
                if Sw(leftelem)<Sleft ||Sleft<Sright
                    Sleft=Sw(leftelem);
                end
            end
        end
    end
    
end
if Sw(leftelem)<=Sleft && Sleft<=Sright && Sright<=Sw(rightelem)
elseif Sw(rightelem)<=Sright && Sright<=Sleft && Sleft<=Sw(leftelem)
else
    disp('erro não é monotono os arugmentos locais')
    pause
end
%%
if sign(dotvn)<0 && abs(dotvn)>1e-10
    % quando a vazão é menor que zero quer dizer que o
    % fluido vai de direita para esquerda
    A= Sright;
    B= Sleft;
    
    AB= 0.5*(A+B);
    %"fw" has three values: fw(Sleft) is fw(1), fw(Sright) is fw(3)
    [fwaux,~,gama,] = twophasevar([A AB B],numcase);
    %There saturation difference bigger than zero
    [dfwdS_rh,dgamadS_rh] = calcdfunctiondS([fwaux(1) fwaux(3)],[gama(1) gama(3)],[A B],0);
    if (AB-A) ~= 0
        dfdSA = (fwaux(2) - fwaux(1))/( AB-A);
        
        %The saturation difference is zeros
    else
        dfdSA = 0;
        
    end  %End of IF
    if (AB-B) ~= 0
        dfdSB= (fwaux(2) - fwaux(3))/(AB-B);
        
        %The saturation difference is zeros
    else
        dfdSB = 0;
        
    end  %End of IF
    if ( dfdSB<= dfwdS_rh  && dfwdS_rh<= dfdSA) %|| ( dfdSB>= dfwdS_rh  && dfwdS_rh>= dfdSA)
        mm=1;
        Srightauxentropy=A;
        Sleftauxentropy =B ;
    else
        mm=0;
    end
elseif sign(dotvn)==0 || abs(dotvn)<1e-10
    
    mm=1;
    Srightauxentropy=Sright;
    Sleftauxentropy =Sleft ;
    
else
    
    A= Sleft;
    B= Sright;
    
    AB= 0.5*(A+B);
    %"fw" has three values: fw(Sleft) is fw(1), fw(Sright) is fw(3)
    [fwaux,~,gama,] = twophasevar([A AB B],numcase);
    [dfwdS_rh,dgamadS_rh] = calcdfunctiondS([fwaux(1) fwaux(3)],[gama(1) gama(3)],[A B],0);
    %There saturation difference bigger than zero
    if (AB-A) ~= 0
        dfdSA = (fwaux(2) - fwaux(1))/( AB-A);
        
        %The saturation difference is zeros
    else
        dfdSA = 0;
        
    end  %End of IF
    if (AB-B) ~= 0
        dfdSB= (fwaux(2) - fwaux(3))/(AB-B);
        
        %The saturation difference is zeros
    else
        dfdSB = 0;
        
    end  %End of IF
    if ( dfdSB<= dfwdS_rh && dfwdS_rh<= dfdSA) %|| ( dfdSB>= dfwdS_rh && dfwdS_rh>= dfdSA)
        mm=1;
        Srightauxentropy=B;
        Sleftauxentropy =A ;
    else
        mm=0;
    end
end

if sign(dotvn)<0 && abs(dotvn)>1e-10 && mm==0
    % quando a vazão é menor que zero quer dizer que o
    % fluido vai de direita para esquerda
    A= Sw(rightelem);
    B= Sw(leftelem);
    AB= 0.5*(A+B);
    %"fw" has three values: fw(Sleft) is fw(1), fw(Sright) is fw(3)
    [fwaux,~,gama,] = twophasevar([A AB B],numcase);
    %There saturation difference bigger than zero
    [dfwdS_rh,dgamadS_rh] = ...
        calcdfunctiondS([fwaux(1) fwaux(3)],[gama(1) gama(3)],[A B],0);
    if (AB-A) ~= 0
        dfdSA = (fwaux(2) - fwaux(1))/( AB-A);
        
        %The saturation difference is zeros
    else
        dfdSA = 0;
        
    end  %End of IF
    if (AB-B) ~= 0
        dfdSB= (fwaux(2) - fwaux(3))/(AB-B);
        
        %The saturation difference is zeros
    else
        dfdSB = 0;
        
    end  %End of IF
    
    if ( dfdSB<= dfwdS_rh  && dfwdS_rh<= dfdSA)||( dfdSB>= dfwdS_rh  && dfwdS_rh>= dfdSA)
        Srightauxentropy=A;
        Sleftauxentropy= B;
    end
elseif (sign(dotvn)==0 || abs(dotvn)<1e-10) && mm==0
    dfdSA=0;
    dfdSB=0;
    dfwdS_rh=0;
    
    Srightauxentropy=Sright;
    Sleftauxentropy= Sleft;
elseif mm==0
    % quando a vazão é menor que zero quer dizer que o
    % fluido vai de esquerda para direita
    A= Sw(leftelem);
    B= Sw(rightelem);
    
    AB= 0.5*(A+B);
    %"fw" has three values: fw(Sleft) is fw(1), fw(Sright) is fw(3)
    [fwaux,~,gama,] = twophasevar([A AB B],numcase);
    [dfwdS_rh,dgamadS_rh] = ...
        calcdfunctiondS([fwaux(1) fwaux(3)],[gama(1) gama(3)],[A B],0);
    %There saturation difference bigger than zero
    if (AB-A) ~= 0
        dfdSA = (fwaux(2) - fwaux(1))/( AB-A);
        
        %The saturation difference is zeros
    else
        dfdSA = 0;
        
    end  %End of IF
    if (AB-B) ~= 0
        dfdSB= (fwaux(2) - fwaux(3))/(AB-B);
        
        %The saturation difference is zeros
    else
        dfdSB = 0;
        
    end  %End of IF
    
    if ( dfdSB<= dfwdS_rh  && dfwdS_rh<= dfdSA) || ( dfdSB>= dfwdS_rh  && dfwdS_rh>= dfdSA)
        
        Sleftauxentropy= A;
        Srightauxentropy=B;
    end
end


Sright=Srightauxentropy;
Sleft= Sleftauxentropy;
end