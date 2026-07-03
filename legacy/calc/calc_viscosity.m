function viscosidade=calc_viscosity(nflagfacec,Sleft,Sright,timelevel,earlysw,Con,timelevelold)
global bedge inedge visc elemarea numcase benchkey

% if timelevel~=1
%
%     viscosidade=1./(0.5-0.2.*earlysw);
%
% else
for i=1:size(bedge,1)+size(inedge,1)
    if i<= size(bedge,1)
        if bedge(i,7)<200
            
            concen=nflagfacec(i,2);
            
        else
            if timelevel==1
                
                concen=0;
            else
                concen=Sleft(i);
                %concen=earlysw(i);
            end
            
        end
        if numcase==251
            expo=((1-concen)*visc(2)^-0.25+(visc(1)^(-0.25))*concen)^4;
        else
            expo=1/(0.5-0.2*concen);
        end
        viscosidade(i,1)=expo;
    else
        lef=inedge(i-size(bedge,1),3);
        rel=inedge(i-size(bedge,1),4);
        if timelevel==1
            concenlef=0;
            concenrel=0;
            %concen=0;
        else
            concenlef=Sleft(i);
            concenrel=Sright(i-size(bedge,1));
            %concen=earlysw(i);
        end
        if numcase==251
            expolef=((1-concenlef)*visc(2)^-0.25+(visc(1)^(-0.25))*concenlef)^4;
            exporel=((1-concenrel)*visc(2)^-0.25+(visc(1)^(-0.25))*concenrel)^4;
        else
            exporel=1/(0.5-0.2*concenlef);
            expolef=1/(0.5-0.2*concenrel);
            %expo=1/(0.5-0.2*concen);
        end
        
        viscosidade(i,1)=(elemarea(lef)*expolef+elemarea(rel)*exporel)/(elemarea(lef)+elemarea(rel));
        %viscosidade(i,1)=expo;
    end
end
%end

end