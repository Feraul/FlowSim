function [M,I]=ferncodes_auxassemblematrixinteriorDMP1(auxiface,M,I,auxparameter,ielem,...
                                       pinterp,norma,beta,weightDMP,auxmobility)
global bedge inedge 
if auxiface >size(bedge,1) 
    auxlef=inedge(auxiface-size(bedge,1),3);
    auxrel=inedge(auxiface-size(bedge,1),4);
    
    %[weight1,weight2]=weightnlfvDMP(kmap,auxiface-size(bedge,1));
    % Calculo das contribuições do elemento a esquerda
     weight1=weightDMP(auxiface-size(bedge,1),1);
     weight2=weightDMP(auxiface-size(bedge,1),2);
    if auxlef==ielem
        elemavaliarlef= auxrel;
        auxweight= weight2;
    elseif auxrel==ielem
        elemavaliarlef= auxlef;
        auxweight= weight1;
    end
    
    auxTransmilef= auxweight*beta*auxparameter*norma;
    
    M(ielem,ielem)=M(ielem,ielem)+ auxTransmilef*auxmobility;
    
    M(ielem,elemavaliarlef)=  M(ielem,elemavaliarlef)- auxTransmilef*auxmobility;
    
else 
    
    auxTransmilef= beta*auxparameter*norma;
    
    M(ielem,ielem)=M(ielem,ielem)+ auxTransmilef*auxmobility;
    
    I(ielem)=I(ielem)+ auxTransmilef*pinterp(auxiface)*auxmobility;
end

end