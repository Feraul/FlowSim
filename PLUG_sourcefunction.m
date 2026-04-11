function [sourcevector] = PLUG_sourcefunction(P,env,time,parmRichardEq)


elem=env.geometry.elem;
centelem=env.geometry.centelem;
elemarea=env.geometry.elemarea;
numcase=env.config.numcase;
if numcase>400
    alpha=parmRichardEq.alpha;
    nvg= parmRichardEq.nvg;
end
n = size(elem,1);

x = centelem(:,1);
y = centelem(:,2);

sourcevector = zeros(n,1);

switch numcase

    case 434
        mask = (y > 95) & (y < 105);
        sourcevector(mask) = -355 .* elemarea(mask);
    case 436

        % Parâmetros
        alpha = parmRichardEq.alpha;
        n   = parmRichardEq.nvg;
        Kabs = 0.03;
        mu    = 1;

        % Coordenadas dos centros dos elementos
        x = centelem(:,1);
        y = centelem(:,2);
        t = time;

        % ========================================================================

% -------------------- Solução p(x,y,t) ----------------------------------
p = -3*t .* x.*(1-x).*y.*(1-y) - 1;

% Derivadas de p
pt1  = -3 * x.*(1-x).*y.*(1-y); % ok

px  = -3*t * (1-2*x).*y.*(1-y);
py  = -3*t * x.*(1-x).*(1-2*y);

grad2 = px.^2 + py.^2;

lap_p = 6*t * ( x.*(1-x) + y.*(1-y) ); %ok

% -------------------- Theta(p) ------------------------------------------
A1 = (-alpha*p).^n;

theta = (1 + A1).^(-(n-1)/n);

% -------------------- dtheta/dp -----------------------------------------
dtheta_dp1 = (n-1)*alpha * (-alpha*p).^(n-1) .* (1 + A1).^(-(2*n-1)/n); %ok

% -------------------- kappa(theta) --------------------------------------
m = n/(n-1);
beta = (n-1)/n;

theta_m = theta.^m;

% regularização numérica
base = 1 - theta_m;

B = base.^beta;
T = 1 - B;

kappa = (Kabs/mu) * sqrt(theta) .* T.^2;

% -------------------- dkappa/dtheta -------------------------------------
dT_dtheta = m*beta * base.^(beta-1) .* theta.^(m-1);

dkappa_dtheta = (Kabs/mu) * (0.5 * theta.^(-1/2) .* T.^2 + ...
                              2 * sqrt(theta) .* T .* dT_dtheta );

% -------------------- Termo fonte f -------------------------------------
f = dtheta_dp1 .* pt1  - kappa .* lap_p - dkappa_dtheta .* dtheta_dp1 .* grad2;


        % ------------------ vetor fonte ----------------------------
        sourcevector = f .* elemarea;

    case 11
        sourcevector = elemarea .* ...
            (-2 .* exp(x.*y) .* (1 + x.^2 + y.^2 + x.*y));

    case {12.1,12.2,12.3,12.4}
        % Define alpha conforme o caso
        if numcase == 12.1, alpha = 1; end
        if numcase == 12.2, alpha = 10; end
        if numcase == 12.3, alpha = 100; end
        if numcase == 12.4, alpha = 1000; end

        mask1 = x < 0;
        mask2 = ~mask1;

        sourcevector(mask1) = elemarea(mask1) .* ...
            (((2*sin(y(mask1)) + cos(y(mask1))) .* ...
            (alpha .* x(mask1))) + sin(y(mask1)));

        sourcevector(mask2) = -elemarea(mask2) .* ...
            ((2*alpha) .* exp(x(mask2)) .* cos(y(mask2)));

    case 13
        sourcevector = elemarea .* ...
            (11*pi^2 .* cos(pi*x) .* cos(pi*y));

    case 14.1
        sourcevector = elemarea .* ...
            (-48*x.^2 + x.*(80 - 64*y) - ...
            48*(y - 1.43426).*(y - 0.232408));

    case 15.3
        epsilon = 1e-3;

        x1 = x + epsilon;
        y1 = y + epsilon;

        sourcevector = -elemarea .* ...
            (sin(pi*x).*sin(pi*y).*((1+epsilon)*pi^2.*(x1.^2 + y1.^2)) + ...
            cos(pi*x).*sin(pi*y).*((1-3*epsilon)*pi.*x1) + ...
            sin(pi*x).*cos(pi*y).*((1-3*epsilon)*pi.*y1) + ...
            cos(pi*x).*cos(pi*y).*(2*pi^2*(1-epsilon).*x1.*y1));

    case 15.4
        mask = (x > 0.125 & x < 0.375) & ...
            (y > 0.125 & y < 0.375);

        sourcevector(mask) = elemarea(mask);

    case 16
        mask1 = x <= 0.5;
        mask2 = x > 0.5;

        sourcevector(mask1) = elemarea(mask1) .* ...
            (2*pi^2 .* cos(pi*x(mask1)) .* sin(pi*y(mask1)));

        sourcevector(mask2) = elemarea(mask2) .* ...
            (9.87059 .* cos(pi*x(mask2)) .* sin(pi*y(mask2)));

    case 20.1
        idx = ismember((1:n)', [121 122]);
        sourcevector(idx) = 1;

    case 20.2
        sourcevector(121) = elemarea(121);

    case {21,21.1}
        x = centelem(:,1);
        y = centelem(:,2);

        cx = cos(pi*x);
        cy = cos(pi*y);
        sx = sin(pi*x);
        sy = sin(pi*y);

        % Término común
        C = 2.46711*(x.^2 - y.^2).*cx.*cy;

        % t1
        t1 = C ...
            - pi^2*(1 + 1.0669*x.^2 + 1.93289*y.^2).*sx.*sy ...
            + 1.57061*x.*sx.*cy ...
            + 6.70353*x.*cx.*sy;

        % t2
        t2 = C ...
            - pi^2*(1 + 1.93289*x.^2 + 1.0669*y.^2).*sx.*sy ...
            + 6.70353*y.*sx.*cy ...
            - 1.57061*y.*cx.*sy;

        % Resultado final
        sourcevector = -(t1 + t2) .* elemarea;

    case 22
        mask = (x >= 7/18 & x <= 11/18) & ...
            (y >= 7/18 & y <= 11/18);

        sourcevector(mask) = 1 ./ vecnorm(centelem(mask,:),2,2);

    case 29
        mask = (x >= 0.125 & x <= 0.375) & ...
            (y >= 0.125 & y <= 0.375);

        sourcevector(mask) = 10 .* elemarea(mask);

    case 336
        mask = x <= 0.5;
        sourcevector(mask) = -20 .* elemarea(mask);
        sourcevector(~mask) = -8 .* elemarea(~mask);

    case {333,335,337,338}
        sourcevector = P .* elemarea;

    case {341,341.1}
        sourcevector = P(:) .* elemarea;

    case 347
        idx = ~ismember((1:n)', wells(:,1));
        sourcevector(idx) = P .* elemarea(idx);

    case 342
        if time ~= 0
            termo1 = (0.00880285 .* x .* exp(-0.0000270492*(x.^2)/time)) ./ (time^(3/2));
            termo2 = (9.5244e-7 .* x .* exp(-0.0000270492*(x.^2)/time)) ./ (time^(3/2));

            sourcevector = elemarea .* ((3.28e-3)*3 .* termo1 - 91.5 .* termo2);
        end

end