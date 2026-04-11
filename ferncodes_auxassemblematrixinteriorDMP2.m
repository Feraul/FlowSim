function [M,I]=ferncodes_auxassemblematrixinteriorDMP2(auxiface,M,I,parameter1,parameter2,ielem1,ielem2,mu1,mu2,weight1,weight2,...
                                       pinterp,norma,beta,ifactual,gamma,weightDMP,auxmobility)
global bedge inedge 
if auxiface >size(bedge,1) && auxiface==ifactual
       
    Transmi1=weight2*gamma*beta*parameter1*norma;
    %% contribuição da transmisibilidade no elemento esquerda
    M(ielem1,ielem1)=M(ielem1,ielem1)+ Transmi1*auxmobility;
    M(ielem1,ielem2)=M(ielem1,ielem2)- Transmi1*auxmobility;
    
elseif auxiface >size(bedge,1) && auxiface~=ifactual
    
    auxlef=inedge(auxiface-size(bedge,1),3);
    auxrel=inedge(auxiface-size(bedge,1),4);
    
    %[weight1,weight2]=weightnlfvDMP(kmap,auxiface-size(bedge,1));
    % Calculo das contribuições do elemento a esquerda
     weight1=weightDMP(auxiface-size(bedge,1),1);
     weight2=weightDMP(auxiface-size(bedge,1),2);
    if auxlef==ielem1
        elemavaliarlef= auxrel;
        auxweight= weight2;
    elseif auxrel==ielem1
        elemavaliarlef= auxlef;
        auxweight= weight1;
    end
    
    auxTransmilef= auxweight*beta*parameter1*norma;
    
    M(ielem1,ielem1)=M(ielem1,ielem1)+ auxTransmilef*auxmobility;
    
    M(ielem1,elemavaliarlef)=  M(ielem1,elemavaliarlef)- auxTransmilef*auxmobility;  
    
else
    
    auxTransmilef= beta*parameter1*norma;
    
    M(ielem1,ielem1)=M(ielem1,ielem1)+ auxTransmilef*auxmobility;
    
    I(ielem1)=I(ielem1)+ auxTransmilef*pinterp(auxiface)*auxmobility;
end

end