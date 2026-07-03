function [q,q_numerico] = fluxo_Richards(x,z,t, parmRichardEq,env,flowrate)
% RICHARDS_FLOW_SOLUTION Calcula pressão e fluxos usando solução analítica
%   Entradas:
%       parmrichardEq - estrutura com campos: alpha, theta_s, theta_r, h_init, Ks
%       centelem       - matriz Nx2 com coordenadas (x, z) dos pontos
%       t              - tempo (escalar)
%   Saídas:
%       psi            - carga de pressão (mesmo tamanho que x, z)
%       qx, qz         - componentes do fluxo de Darcy

    % --- Parâmetros físicos e geométricos ---
    Ks = 0.1;
    alpha = parmRichardEq.alpha;
    theta_s = parmRichardEq.theta_s;
    theta_r = parmRichardEq.theta_r;
    L=parmRichardEq.valuecontor;
    l=L;
    
    %l = 15.24;   % comprimento na direção x (período)
    %L = 15.24;   % comprimento na direção z (domínio)
    % --- Pré-cálculo de xi = exp(alpha * h) ---
    % Se h for escalar, xi será escalar (broadcast automático)
    xi = exp(-alpha *L); % condicao inicial
    
    % --- Parâmetros da solução ---
    beta1 = sqrt(alpha^2/4 + (pi/l)^2);
    d = alpha * (theta_s - theta_r) / Ks;   % difusividade
    
    % --- Número de termos da série (convergência) ---
    N_terms = 200;
    j = (1:N_terms)';
    lambda = (pi / L) .* j;
    gamma = (1/d).* (beta1^2 + lambda.^2);
    
    % Coeficientes da série temporal
    
    coef_series = (2/(L*d)) * ((-1).^j) .* (lambda ./ gamma) .* exp(-gamma * t);
    % coef_series: vetor coluna N_terms x 1
    
    % --- Solução estacionária (parte transiente excluída) ---
    sinh_beta1_z = sinh(beta1 * z);
    sinh_beta1_L = sinh(beta1 * L);
    F_stat = sinh_beta1_z ./ sinh_beta1_L;
    dFdz_stat = beta1 * cosh(beta1 * z) ./ sinh_beta1_L;
    
    % --- Contribuição da série (parte dependente do tempo) ---
    % Matrizes de senos e cossenos: N_terms x N_points
    sin_lz = sin(lambda .* z');
    cos_lz = cos(lambda .* z');
    
    % Soma sobre j: resultado é vetor coluna (N_points x 1)
    F_series = sin_lz' * coef_series;
    dFdz_series = (lambda .* cos_lz)' * coef_series;
    
    % --- Funções F e dF/dz completas ---
    F = F_stat + F_series;
    dFdz = dFdz_stat + dFdz_series;
    
    % --- Funções trigonométricas em x ---
    sin_pix_l = sin(pi * x ./ l);
    cos_pix_l = cos(pi * x ./ l);
    exp_factor = exp((alpha/2).* (L - z));
    
    % --- Carga de pressão bar_psi e derivadas ---
    bar_psi = (1 - xi) .* sin_pix_l .* exp_factor .* F;
    
    dbar_psi_dx = (1 - xi) .* (pi/l) .* cos_pix_l .* exp_factor .* F;
    dbar_psi_dz = (1 - xi) .* sin_pix_l .* exp_factor .* (dFdz - alpha/2 * F);
    
    % --- Carga de pressão total psi (transformada inversa) ---
    psi = (1/alpha) * log(xi + bar_psi);
    
    % --- Condutividade hidráulica K(psi) = Ks * exp(alpha*psi) ---
    Kpsi = Ks * exp(alpha * psi);   % equivalente a Ks * (xi + bar_psi)
    
    % --- Gradiente de (psi + z) ---
    grad_psi_x = (1/alpha) * dbar_psi_dx ./ (xi + bar_psi);
    grad_psi_z = (1/alpha) * dbar_psi_dz ./ (xi + bar_psi) + 1;
    
    % --- Fluxos de Darcy ---
    qx = -Kpsi .* grad_psi_x;
    qz = -Kpsi .* grad_psi_z;

    q=sum([qx,qz].*env.geometry.unitnormals(:,1:2),2);
    q_numerico=flowrate./vecnorm(env.geometry.normals(:,1:2),2,2);

end