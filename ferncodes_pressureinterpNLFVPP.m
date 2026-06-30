function [pressurinterp,cinterp]=ferncodes_pressureinterpNLFVPP(p,premethod,env)


coord=env.geometry.coord;
esurn1=env.geometry.esurn1;
esurn2=env.geometry.esurn2;
numcase=env.config.numcase;
w=premethod.MPFAD.weight;
s=premethod.MPFAD.s;
nflag=env.config.nflag;

nno = size(coord,1);
pressurinterp = zeros(nno,1);
cinterp = 0;

%% --- construir índices nodo-elemento ---
counts = diff(esurn2);                 % a diferenca de elementos consicutivos
node_ids = repelem((1:nno)',counts);   % nodo asociado a cada entrada CSR
elem_ids = esurn1;                     % elementos

%% --- interpolación presión ---
vals = w' .* p(elem_ids);
sum_press = accumarray(node_ids,vals,[nno 1]);

pressurinterp = sum_press;

% añadir contribución s si flag==202
mask202 = (nflag(:,1)==202);
pressurinterp(mask202) = pressurinterp(mask202) + s(mask202);

% nodos con condición Dirichlet
mask_dir = (nflag(:,1)<=200);
pressurinterp(mask_dir) = nflag(mask_dir,2);

%% --- interpolación concentración ---
if (200<numcase && numcase<300) || (379<numcase && numcase<400)

    valsc = wightc .* Con(elem_ids);
    sum_con = accumarray(node_ids,valsc,[nno 1]);

    cinterp = sum_con;

    mask202c = (nflagc(:,1)==202);
    cinterp(mask202c) = cinterp(mask202c) + sc(mask202c);

    mask_dirc = (nflagc(:,1)<=200);
    cinterp(mask_dirc) = nflagc(mask_dirc,2);

    if numcase<200
        cinterp = max(0,min(1,cinterp));
    end

end

end