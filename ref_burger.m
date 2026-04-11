function ue = ref_burger(u0,x,t)
u_guess = u0(:); ue = zeros(length(x(:)),1);       
for i=1: length(x(:))
    %f1 = @(v) v - 0.5 - sin(x(i)-t*v);
    f1 = @(v) v -(1/(2*pi))*sin(2*pi*(x(i)-t*v));
    ue(i) = fzero(f1, u_guess(i));  
end
return