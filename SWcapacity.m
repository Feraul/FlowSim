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

elseif numcase==437
    
    Capw(idx)=alpha.*(theta_s - theta_r).*exp(alpha.*h(idx));

else
    
   A = theta_s - theta_r;

    % Termos auxiliares vetorizados
    term1 = (-alpha .* h(idx)).^(pp-1);      % N×1
    term2 = (1 + (-alpha .* h(idx)).^pp).^(-q-1);  % N×1

    % Derivada final
    Capw(idx) = A .* q .* pp .* alpha .* term1 .* term2;

end