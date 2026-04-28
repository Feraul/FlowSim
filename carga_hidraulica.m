function [h_exact]=carga_hidraulica(x,z,t,parmRichardEq)

%------------------------------------------------------------------
% x : vetor (nx)
% z : vetor (nz)
% t : escalar
% N : número de termos do somatório

% x, z já são vetores coluna de mesmo tamanho


%L = 15.24;
%l = 15.24;
L=parmRichardEq.valuecontor;
l=L;

theta_s = parmRichardEq.theta_s;
theta_r = parmRichardEq.theta_r;
Ks      = 0.1;
alpha   = parmRichardEq.alpha;


xi = exp(-alpha * L);
N  = 200;

% Parâmetros auxiliares
beta1 = sqrt((alpha^2)/4 + (pi/l)^2);
d     = alpha*(theta_s - theta_r)/Ks;

% --- Termo base (N×1) ---
base = (1 - xi) .* sin(pi*x/l) .* exp((alpha/2) .* (L - z));

% --- Termo 1 (N×1) ---
termo1 = sinh(beta1 * z) ./ sinh(beta1 * L);

% --- Termo 2 (somatório ponto a ponto) ---
j = (1:N)';                 % N×1
lambda_j = (pi/L) * j;
gamma_j  = (beta1^2 + lambda_j.^2)/d;

% sin(lambda_j * zᵢ) → (N × nPoints)
sin_lz = sin(lambda_j .* z');   % broadcasting

% exp(-gamma_j * t) → (N × 1)
exp_gt = exp(-gamma_j * t);

% Combinação → (N × nPoints)
S = sin_lz .* exp_gt;

% Somatório → (1 × nPoints)
termo2 = (2/(L*d)) * sum( (-1).^j .* (lambda_j ./ gamma_j) .* S , 1);

% Transpor para vetor coluna (nPoints × 1)
termo2 = termo2';

% --- ψ̄1(x,z,t) ---
psi_bar = base .* (termo1 + termo2);

% --- h_exact ---
h_exact = (1/alpha) .* log(xi + psi_bar);
% estado estacionario
h_ss =(1/alpha).*log(xi+ (1-xi) .* sin(pi*x/l) .* exp((alpha/2) .* (L - z)).*termo1);

end