function [Capw] = SWcapacity(h,parmRichardEq, env)

numcase=env.config.numcase;
theta_s=parmRichardEq.theta_s;
theta_r=parmRichardEq.theta_r;
alpha=parmRichardEq.alpha;
pp=parmRichardEq.pp;
q=parmRichardEq.q;
nvg=parmRichardEq.nvg;
% calculo da capacidade da agua

Capw = zeros(size(h));      % inicializa
idx = (h <= 0);             % condição usada em todos os casos

if numcase == 431
    
    Capw(idx) = -q*pp*(alpha^pp) * ((theta_s-theta_r) ./ ...
        ((1 + (alpha*abs(h(idx))).^pp).^(q+1))) .* ...
        (abs(h(idx)).^pp) .* (h(idx).^-1);

elseif numcase == 436
    
    Capw(idx) = alpha*(nvg-1) .* ((-alpha*h(idx)).^(nvg-1)) .* ...
        (((-alpha*h(idx)).^nvg + 1).^((1/nvg) - 2));

else
    
    Capw(idx) = (theta_s-theta_r)*pp*q*alpha .* ((-alpha*h(idx)).^(pp-1)) .* ...
        (1 + (-alpha*h(idx)).^pp).^((1/pp) - 2);

end