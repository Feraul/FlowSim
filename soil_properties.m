function [M,I]=soil_properties(M,I,parms,flowresultZ,env)


Area=env.geometry.elemarea;

h_n=parms.h_init;
h_m=parms.h_old;
dt=parms.dt;
invdt = 1/dt;

%-------------------------------------------------------------------------
theta_n  = env.benchmark.calcularTheta(h_n, parms);
theta_m  = env.benchmark.calcularTheta(h_m, parms);
dthetadh = env.benchmark.calcularCapacidade(h_m, parms);

Dtheta=theta_m - theta_n;

coef=invdt.*dthetadh.*Area(:);
coefI=dthetadh.*Area(:);

M = M + spdiags(coef, 0, size(M,1), size(M,2));

I = I + invdt .* (coefI .* h_m - Dtheta.*Area(:)) - flowresultZ;
end