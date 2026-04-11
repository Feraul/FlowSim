function [Sat_max, Sat_min]=ferncodes_Saturation_max_min(ielem,Sw)
global elem esurn1 esurn2 
if elem(ielem,4)~=0
    m=4;
else
    m=3;
end

    r=1;
    for j=elem(ielem,1:m)
        list=esurn1(esurn2(j)+1:esurn2(j+1));
        
        for k=1:size(list,1)
            if list(k)~=ielem %não inclui o elemento em questão
                a(r,1)=Sw(list(k));
                r=r+1;
            end
        end
        
    end
    Sat_max=max(max([a]));
    Sat_min=min(min([a]));
end
