function [M,I]=ferncodes_auxassemblematrixcontourDMP(auxiface,M,I,auxparameter,...
                                             ielem,pinterp,normcont,weightDMP,auxmobility)
global bedge inedge 
if auxiface >size(bedge,1)
    auxlef=inedge(auxiface-size(bedge,1),3);
    auxrel=inedge(auxiface-size(bedge,1),4);
    %[weight1, weight2]=weightnlfvDMP(kmap,auxiface-size(bedge,1));
    weight1=weightDMP(auxiface-size(bedge,1),1);
    weight2=weightDMP(auxiface-size(bedge,1),2);
    if auxlef==ielem
        elemavaliar=auxrel;
        auxweight= weight2;
    elseif auxrel==ielem
        elemavaliar=auxlef;
        auxweight=weight1;
    end
    
    Transmicont= auxweight*auxparameter*normcont;
    
    M(ielem,ielem)=M(ielem,ielem)+Transmicont*auxmobility;
    
    M(ielem,elemavaliar)=M(ielem,elemavaliar)-Transmicont*auxmobility;
    
else
    Transmicont= auxparameter*normcont;
    
    M(ielem,ielem)=M(ielem,ielem)+ Transmicont*auxmobility;
    
    I(ielem)=I(ielem)+ Transmicont*pinterp(auxiface)*auxmobility;
end

end