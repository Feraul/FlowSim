function u0 = IC(x,x_c,Np,ICcase)
% Create vector u0 with an Initial Condition. 8 cases are available.
%**************************************************************************
%
% 4 cases available: 
% {1} Gaussian wave for Advection problems
% {2} Gaussian wave for Diffusion problems
% {3} Lifted Sinusoidal wave
% {4} Centered Sinusoidal wave
% {5} HyperTangent for Burgers' Equation problem
% {6} Riemann IC
% {7} HyperTangent
% {8} Square Jump
%**************************************************************************
% References: 
% [1] High-Order Methods for Diffusion Equation with Energy Stable Flux
% Reconstruction Scheme. K. Ou, P. Vincent, and A. Jameson. AIAA 2011-46.
%
%**************************************************************************

% Create the selected IC
switch ICcase
    case 1 % Gaussian wave for [-1,1] 
           % for advection test with periodic BCs. See ref [1].
        xmid = (x(end) + x(1))/2;
        u0 = exp(-20*(x-xmid).^2);
        
    case 2 % Gaussian wave for [-3,3] 
           % for diffusion test with Dirichlet BCs. See ref [1].
        mu = 0.01;
        xmid = (x(end) + x(1))/2;
        u0 = exp(-(x-xmid).^2/(4*mu));
               
    case 3 % Centered sinusoidal wave for [-1,1] for advection.
        u0 = -sin(pi*x);
        
	case 4 % Lifted sinusoidal wave for [0,2*pi] for advection.
        u0 = 0.5 + sin(x);
        
    case 5 % Hyperbolic Tangent in [-,2,2] 
           % for Viscous Burgers' equation with Dirichlet BCs. See ref [1].
        mu = 0.02;
        u0 = 0.5*(1-tanh(x/(4*mu)));
        
    case 6 % Riemann problem
        % u = 1 for x <  x_mid
        % u = 2 for x >= x_mid
        u0 = ones(size(x));
        xmid = (x(end)+x(1))/2;
        rhs = find(x<=xmid);
        u0(rhs) = 2;
    
    case 7 % Hyperbolic Tangent
        % u = 1 for x <  x_mid
        % u = 0 for x >= x_mid
        a = x(1); b = x(end);
        xi = (4-(-4))/(b-a)*(x - a) - 4;
        u0 = 1/2*(tanh(-4*xi)+1);
        
    case 8 % Square Jump
        xmid = 0; % center of square
        u0 = heaviside(x-xmid+1/4) - heaviside(x+xmid-1/4) + 1;
    
    case 9 % step [-1 1]
        range1 = x<0;
        range2 = x>=0;
        u0 = -0.5.*range1 + 1.*range2;
        
    case 10
        range1 = x>=0;
        u0 = (1e-6).*range1; %u0(1,1)=1;
        
    case 11 % step [-1 1]
        range1 = x>=-0.5 & x<=1e-6;
        range2 = x<-0.5 & x>1e-6;
        u0 = 1.*range1 + 1e-6.*range2;
        
        
    case 12 % step [0 1]
        range1 = x>0.25 & x<0.75;
        range2 = x<=0.25 & x>=0.75;
        u0 = range1 + 0*range2;
        
    case 13 % bl_[-1 1] Sw<0
        
        range1 = x>=-0.5 & x<=1e-6;
        range2 = x<-0.5; 
        range3 = x>1e-6;
        u0 = 1.*range1 + -0.5*range2 + 0*range3;
        
    case 14
        sLB= 0.9999; sRB= 0.0001;
        %% Initial solution saturation
        u0 = sRB + (sLB-sRB)*(x_c<=0.0); u0 = kron(u0, ones(Np,1));
        
    case 15
        
         u0 = (1/(2*pi))*sin(2*pi*x);
    otherwise
        error('case not in the list')
end
