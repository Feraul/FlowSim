


function [env,parmRichardEq] = calcnormk(env, parmRichardEq,time)

%Define the norm of permeability tensor
%Obtain "kmap" for each case
[env,parmRichardEq] = PLUG_kfunction(env,parmRichardEq,time);
%Swept all elements
pointer = env.geometry.elem(:,5);
if env.config.numcase < 400 || isempty(parmRichardEq)
    permcompon = env.config.perm(:,2:5);
else
    permcompon = parmRichardEq.auxperm(:,2:5);
end

% Norma Frobenius vetorizada
tmp = sqrt(sum(permcompon.^2, 2));   % 1×1×N

% Converter para vetor coluna N×1
normk = reshape(tmp, [], 1);

% Garantir que normk tenha o mesmo número de linhas que env.config.perm
% (caso pointer não seja 1:N)
fullnorm = zeros(size(env.config.perm,1),1);
fullnorm(pointer) = normk;

normk = fullnorm;
if env.config.numcase <400
    env.geometry.normperm=normk;
else
    parmRichardEq.normperm=normk;
end

end%End of FOR