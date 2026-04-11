function [M,I]=soil_properties(M,I,parmRichardEq,preMPFAD,env)


auxelemarea=env.geometry.elemarea;
auxflowresultZ=preMPFAD.flowresultZ;

h_n=parmRichardEq.h_init;
h_kickoff=parmRichardEq.h_old;
dt=parmRichardEq.dt;

%-------------------------------------------------------------------------
theta_n=thetafunction(h_n,parmRichardEq,env);

theta_n1=thetafunction(h_kickoff,parmRichardEq,env);

soilwatercapac=SWcapacity(h_kickoff,parmRichardEq,env);

diff_theta=theta_n1 - theta_n;

M=M+(dt^-1)*diag(soilwatercapac.*auxelemarea(:));
I=I+(dt^-1)*diag(soilwatercapac.*auxelemarea(:))*h_kickoff-...
    (dt^-1)*diff_theta.*auxelemarea(:)+auxflowresultZ;

end