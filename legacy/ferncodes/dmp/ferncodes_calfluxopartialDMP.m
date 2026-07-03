function [partialflux,auxparameter]=ferncodes_calfluxopartialDMP(iface1,iface2,...
              auxparameter1, auxparameter2,gamma,ielem,pinterp,ifactual,p,norma,auxmobility)

if iface1==ifactual
    % Fluxo parcial do elemento a esquerda
    
    partialflux=auxmobility*norma*(gamma*auxparameter1*(p(ielem)-pinterp(iface1))+ auxparameter2*(p(ielem)-pinterp(iface2)));
    auxparameter=auxparameter1;
elseif iface2==ifactual
    % Fluxo parcial do elemento a esquerda
    partialflux=auxmobility*norma*(gamma*auxparameter2*(p(ielem)-pinterp(iface2))+ auxparameter1*(p(ielem)-pinterp(iface1)));
    auxparameter=auxparameter2;
else
    % Fluxo parcial do elemento a esquerda
    partialflux=auxmobility*norma*(auxparameter1*(p(ielem)-pinterp(iface1))+ auxparameter2*(p(ielem)-pinterp(iface2)));
    auxparameter=0;
end
if abs(partialflux)<1e-20
    partialflux=0;
end
end
