% a = [16 32 64 128 256];
% %Erro 1a ordem
% b = [5.42e-2 3.29e-2 2.78e-2 1.23e-2 9.33e-3];
% 
% %Erro 2a Ordem MUSCL-va, extrap. Sw, Euler
% c = [2.79e-2 2.04e-2 9.50e-3 6.11e-3 2.47e-3];
% %Erro 2a Ordem MUSCL-va, extrap. fw, Euler
% d = [3.11e-2 1.81e-2 1.10e-2 4.81e-3 2.55e-3];
% 
% %Erro 2a Ordem MUSCL-MLP, extrap. Sw, Euler
% e = [4.25e-2 3.20e-2 1.04e-2 1.15e-2 4.15e-3];
% %Erro 2a Ordem MUSCL-MLP, extrap. fw, Euler
% f = [2.96e-2 2.27e-2 9.40e-3 8.69e-3 3.06e-3];
% 
% %Erro 2a Ordem MUSCL-MLPvk, extrap. Sw, Euler
% g = [2.61e-2 2.00e-2 8.36e-3 6.19e-3 2.24e-3];
% %Erro 2a Ordem MUSCL-MLPvk, extrap. fw, Euler
% h = [2.43e-2 1.46e-2 9.24e-3 4.29e-3 2.27e-3];
% 
% %Erro 2a Ordem MUSCL-va, extrap. Sw, RK2
% i = [2.81e-2 1.75e-2 1.17e-2 4.92e-3 3.20e-3];
% %Erro 2a Ordem MUSCL-MLPvk, extrap. fw, RK2
% j = [3.56e-2 1.92e-2 1.45e-2 5.41e-3 3.93e-3];

%--------------------------------------------------------------

% a = [16 32 64 128 256];
% %Erro 1a ordem
% b = [5.42e-2 3.29e-2 2.78e-2 1.23e-2 9.33e-3];
% 
% %Erro 2a Ordem MUSCL-va, extrap. Sw, Euler
% c = [2.79e-2 2.04e-2 9.50e-3 6.11e-3 2.47e-3];
% 
% %Erro 3a Ordem MUSCL-MLP, extrap. Sw, Euler
% d = [4.43e-2 3.35e-2 1.07e-2 1.14e-2 4.06e-3];
% %Erro 3a Ordem MUSCL-MLP, extrap. fw, Euler
% e = [2.69e-2 1.87e-2 9.30e-3 6.34e-3 2.57e-3];
% 
% %Erro 3a Ordem MUSCL-MLPvk, extrap. Sw, Euler
% f = [2.92e-2 2.20e-2 8.93e-3 6.55e-3 2.35e-3];
% %Erro 3a Ordem MUSCL-MLPvk, extrap. fw, Euler
% g = [3.10e-2 1.68e-2 1.19e-2 4.36e-3 2.79e-3];
% 
% %Erro 3a Ordem MUSCL-MLP, extrap. fw, RK2
% h = [2.72e-2 1.68e-2 9.95e-3 4.84e-3 2.51e-3];

%---------------------------------------------------------------


% a = [16 32 64 128 256];
% %Erro 1a ordem
% b = [5.42e-2 3.29e-2 2.78e-2 1.23e-2 9.33e-3];
% 
% %Erro 2a Ordem MUSCL-va, extrap. Sw, Euler
% c = [2.79e-2 2.04e-2 9.50e-3 6.11e-3 2.47e-3];
% 
% %Erro 4a Ordem, EV, MUSCL-MLP, extrap. Sw, Euler
% d = [6.13e-2 5.29e-2 2.95e-2 2.81e-2 4.59e-3];
% %Erro 4a Ordem, EV, MUSCL-MLP, extrap. fw, Euler
% e = [2.42e-2 1.61e-2 9.20e-3 4.68e-3 2.12e-3];
% 
% %Erro 4a Ordem, EC, MUSCL-MLP, extrap. Sw, Euler
% f = [4.71e-2 2.97e-2 2.79e-2 7.74e-3 3.46e-3];
% %Erro 4a Ordem, EC, MUSCL-MLP, extrap. fw, Euler
% g = [2.48e-2 1.42e-2 1.14e-2 4.53e-3 3.34e-3];
% 
% %Erro 4a Ordem, EV, MUSCL-MLPvk, extrap. Sw, Euler
% h = [5.38e-2 4.85e-2 2.76e-2 2.62e-2 3.62e-3];
% %Erro 4a Ordem, EV, MUSCL-MLPvk, extrap. fw, Euler
% i = [2.66e-2 1.50e-2 1.06e-2 4.30e-3 2.75e-3];
% 
% %Erro 4a Ordem, EC, MUSCL-MLPvk, extrap. Sw, Euler
% j = [3.64e-2 2.33e-2 1.17e-2 6.46e-3 3.42e-3];
% %Erro 4a Ordem, EC, MUSCL-MLPvk, extrap. fw, Euler
% k = [2.47e-2 1.51e-2 1.41e-2 5.24e-3 3.79e-3];
% 
% %Erro 4a Ordem, EV, MUSCL-MLP, extrap. fw, RK2
% l = [2.69e-2 1.54e-2 1.05e-2 4.10e-3 2.77e-3];


a = [16 32 64 128 256];
%Erro 1a ordem
b = [5.42e-2 3.29e-2 2.78e-2 1.23e-2 9.33e-3];

%Erro 2a Ordem MUSCL-MLPvk, extrap. fw, Euler
c = [2.43e-2 1.46e-2 9.24e-3 4.29e-3 2.27e-3];
%Erro 3a Ordem MUSCL-MLP, extrap. fw, RK2
d = [2.72e-2 1.68e-2 9.95e-3 4.84e-3 2.51e-3];
%Erro 4a Ordem, EV, MUSCL-MLP, extrap. fw, Euler
e = [2.42e-2 1.61e-2 9.20e-3 4.68e-3 2.12e-3];




% plot(log2(a),log2(b),'-oc','Linewidth',2)
% hold on
% plot(log2(a),log2(c),'-sb')
% plot(log2(a),log2(d),'--db')
% plot(log2(a),log2(e),'-vr')
% plot(log2(a),log2(f),'--xr')
% plot(log2(a),log2(g),'->k')
% plot(log2(a),log2(h),'--pk')
% plot(log2(a),log2(i),'-+g')
% plot(log2(a),log2(j),'--xg')
% hold off

% plot(log2(a),log2(b),'-oc','Linewidth',2)
% hold on
% plot(log2(a),log2(c),'--sc','Linewidth',2)
% plot(log2(a),log2(d),'-dk')
% plot(log2(a),log2(e),'--vk')
% plot(log2(a),log2(f),'-xr')
% plot(log2(a),log2(g),'-->r')
% plot(log2(a),log2(h),'-pb')
% hold off

% plot(log2(a),log2(b),'-oc','Linewidth',2)
% hold on
% plot(log2(a),log2(c),'--sc','Linewidth',2)
% plot(log2(a),log2(d),'-dk')
% plot(log2(a),log2(e),'--vk')
% plot(log2(a),log2(f),'-xr')
% plot(log2(a),log2(g),'-->r')
% plot(log2(a),log2(h),'-.b')
% plot(log2(a),log2(i),'--+b')
% plot(log2(a),log2(j),'-sg')
% plot(log2(a),log2(k),'--pg')
% plot(log2(a),log2(l),'-py')
% hold off


plot(log2(a),log2(b),'-oc','Linewidth',2)
hold on
plot(log2(a),log2(c),'-sk')
plot(log2(a),log2(d),'--pb')
plot(log2(a),log2(e),'->r')
hold off


grid on;

%Corte de Água
title ('Taxas de Convergência para norma L_1')
xlabel ('Log2(1/h)');
ylabel ('Erro (L_1)');
xlim ([4 8]);
ylim ([-9 -4]);

% legend ('1^a Ordem',...
%     '2^a Ordem (extrap. em Sw, vanAlbada), Euler',...
%     '2^a Ordem (extrap. em fw, vanAlbada), Euler',...
%     '2^a Ordem (extrap. em Sw, MLP), Euler',...
%     '2^a Ordem (extrap. em fw, MLP), Euler',...
%     '2^a Ordem (extrap. em Sw, MLP-vk), Euler',...
%     '2^a Ordem (extrap. em fw, MLP-vk), Euler',...
%     '2^a Ordem (extrap. em Sw, vanAlbada), RK2',...
%     '2^a Ordem (extrap. em fw, MLP-vk), RK2');

% legend ('1^a Ordem',...
%     '2^a Ordem (extrap. em Sw, vanAlbada), Euler',...
%     '3^a Ordem (extrap. em Sw, MLP), Euler',...
%     '3^a Ordem (extrap. em fw, MLP), Euler',...
%     '3^a Ordem (extrap. em Sw, MLP-vk), Euler',...
%     '3^a Ordem (extrap. em fw, MLP-vk), Euler',...
%     '3^a Ordem (extrap. em fw, MLP), RK2');


% legend ('1^a Ordem',...
%     '2^a Ordem (extrap. em Sw, vanAlbada), Euler',...
%     '4^a Ordem, EV (extrap. em Sw, MLP), Euler',...
%     '4^a Ordem, EV (extrap. em fw, MLP), Euler',...
%     '4^a Ordem, EC (extrap. em Sw, MLP), Euler',...
%     '4^a Ordem, EC (extrap. em fw, MLP), Euler',...
%     '4^a Ordem, EV (extrap. em Sw, MLP-vk), Euler',...
%     '4^a Ordem, EV (extrap. em fw, MLP-vk), Euler',...
%     '4^a Ordem, EC (extrap. em Sw, MLP-vk), Euler',...
%     '4^a Ordem, EC (extrap. em fw, MLP-vk), Euler',...
%     '4^a Ordem (extrap. em fw, MLP), RK2');


legend ('1^a Ordem',...
    '2^a Ordem (extrap. em fw, MLP-vk), Euler',...
    '3^a Ordem, (extrap. em fw, MLP), RK2',...
    '4^a Ordem, EV (extrap. em fw, MLP), Euler');
