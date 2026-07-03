function [numflux, Saproxima]=riemannsolver(Sright,Sleft,method,bedgesize, inedg,dotvn,dotvg,charvel_rh,dotdif)

if strcmp(method,'upwd')
   
    
        
%         ve_mais = max([charvel_rh; 0],[],1);
%         ve_menos= min([charvel_rh; 0],[],1);
%     
%         numflux= ve_mais*Sleft + ve_menos*Sright + dotdif;
    
           if charvel_rh > 0 || charvel_rh==0
                 %Calculate the numerical flux through interface
                 numflux = Sleft*dotvn + dotdif ;
                %Fill "earlysw"
                Saproxima=Sleft;
                
                %It uses the saturation on the right
            else
                
                 %Calculate the numerical flux through interface
                 numflux = Sright*dotvn + dotdif;
                %Fill "earlysw"
                Saproxima = Sright;
                
            end  %End of IF (Upwind flux)
   
end
end