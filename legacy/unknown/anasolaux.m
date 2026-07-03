
function [analsol]=anasolaux(velmedio,Dmedio,t)

global elem centelem
if t<=1
    for i=1:size(elem,1)
        
        x=centelem(i,1);
        
        if x>2.5
            A1= 0.5*(erfc((x-velmedio*t)/(2*sqrt(Dmedio*t)))+exp((velmedio*x)/(Dmedio))*erfc((x+velmedio*t)/(2*sqrt(Dmedio*t))));
            
            analsol(i,1)= A1;
        else
            A1= 0.5*(erfc((x-velmedio*t)/(2*sqrt(Dmedio*t)))+ exp((velmedio*x)/(Dmedio))*erfc((x+velmedio*t)/(2*sqrt(Dmedio*t))));
            
            analsol(i,1)= A1;
        end
        
        
    end
else
    for i=1:size(elem,1)
        
        x=centelem(i,1);

        if x>2.7
            A1= 0.5*(erfc((x-velmedio*t)/(2*sqrt(Dmedio*t)))+        (exp((velmedio*x)/(Dmedio)))*erfc((x+velmedio*t)/(2*sqrt(Dmedio*t))));
            A2= 0.5*(erfc((x-velmedio*(t-1))/(2*sqrt(Dmedio*(t-1))))+(exp((velmedio*x)/(Dmedio)))*erfc((x+velmedio*(t-1))/(2*sqrt(Dmedio*(t-1)))));
            analsol(i,1)= A1-A2;
        else
            A1= 0.5*(erfc((x-velmedio*t)/(2*sqrt(Dmedio*t)))+        (exp((velmedio*x)/(Dmedio)))*erfc((x+velmedio*t)/(2*sqrt(Dmedio*t))));
            A2= 0.5*(erfc((x-velmedio*(t-1))/(2*sqrt(Dmedio*(t-1))))+(exp((velmedio*x)/(Dmedio)))*erfc((x+velmedio*(t-1))/(2*sqrt(Dmedio*(t-1)))));
            analsol(i,1)= A1-A2;
        end
        
    end
    
end

end