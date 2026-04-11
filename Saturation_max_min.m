function [Sat_max, Sat_min]=Saturation_max_min(ielem,Sw)
global elem esurn1 esurn2
    n1=elem(ielem,1);
    n2=elem(ielem,2);
    n3=elem(ielem,3);
    r=1;
    for j=[n1,n2,n3]
        list=esurn1(esurn2(j)+1:esurn2(j+1));
        
        for k=1:size(list,1)
            %if list(k)~=ielem
                a(k,r)=Sw(list(k));
            %end
        end
        r=r+1;
    end
    Sat_max=max(max([a]));
    Sat_min=min(min([a]));
end
