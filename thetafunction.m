
% calculo do parametro de van Genutchen: conteudo de agua
function  theta=thetafunction(h,parmRichardEq,env)


auxnumcase=env.config.numcase;
alpha=parmRichardEq.alpha;
nvg=parmRichardEq.nvg;
theta_s= parmRichardEq.theta_s;
theta_r=parmRichardEq.theta_r;
pp=parmRichardEq.pp;
q=parmRichardEq.q;

theta = zeros(size(h));  % inicializa o vetor

if auxnumcase == 436
    idx_neg = h < 0;          % índices onde h < 0
    idx_pos = ~idx_neg;       % índices onde h >= 0

    % h < 0
    theta(idx_neg) = (1 + (-alpha*h(idx_neg)).^nvg).^(-(nvg-1)/nvg);

    % h >= 0
    theta(idx_pos) = 1;
elseif auxnumcase==437
    idx_neg = h < 0;          % índices onde h < 0
    idx_pos = ~idx_neg;       % índices onde h >= 0

    theta(idx_neg)=theta_r + (theta_s - theta_r).*exp(alpha*h(idx_neg));
else
    idx_neg = h < 0;
    idx_pos = ~idx_neg;

    % h < 0
    theta(idx_neg) = theta_r + ((theta_s - theta_r) ./ ((1 + (-alpha*h(idx_neg)).^pp).^q));

    % h >= 0
    theta(idx_pos) = theta_s;
end
end