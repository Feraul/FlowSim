function [MBE, MBE_rel, Rloc, L1, L2] = ferncodes_MBE(theta_old, theta_new, flowrate, flowresult, ...
                                               sourcevector, elemarea, dt, bedgesize, s)
%FERNCODES_MBE  Erro de balanco de massa (MBE) por passo de tempo.
%
%   Rloc(i) = elemarea(i)*(theta_new(i)-theta_old(i)) - dt*(sourcevector(i) - s*flowresult(i))
%
%   Para MPFA-D neste projeto (soil_properties.m + ferncodes_globalmatrix_MPFAD.m),
%   a convencao de sinal correta e s = +1: a matriz M inclui o termo temporal via
%   soil_properties, e flowresult(i) == (M_difusiva*h - I_difusiva)_i (mesma
%   convencao de Kde/Ded/Kn/Kt usada nas duas rotinas). No ponto de convergencia
%   do nao-linear (Picard/L-scheme):
%
%       elemarea(i)*(theta_new(i)-theta_old(i)) = dt*(sourcevector(i) - flowresult(i))
%
%   Entradas:
%     theta_old, theta_new : theta(h) nos dois niveis de tempo deste passo
%                             (theta_old = theta(h_n) usado em soil_properties,
%                              theta_new = theta(h_new) numerico convergido)
%     flowrate              : vetor de fluxo por face (bedge + inedge), saida de
%                              ferncodes_flowrate / calcularFlowrate
%     flowresult             : fluxo liquido por elemento (mesma saida)
%     sourcevector           : termo-fonte * elemarea, MESMO tempo usado no solve
%                              deste passo (PLUG_sourcefunction)
%     elemarea, dt, bedgesize: geometria/passo de tempo
%     s                      : sinal de acoplamento (+1 para este projeto)
%
%   Saidas:
%     MBE      : soma global do residuo (deve ser ~0 na tolerancia do solver)
%     MBE_rel  : MBE normalizado (adimensional)
%     Rloc     : residuo por elemento — use max(abs(Rloc)) para detectar erro
%                localizado (ex.: consistencia MPFA-D em malha distorcida)
%     L1       : norma L1 intensiva (media ponderada por area de |Rloc|/area),
%                analoga a "MPFA-D L1" nos graficos de convergencia log-log
%     L2       : norma L2 intensiva (RMS de Rloc/area, ponderado por area),
%                analoga a "MPFA-D L2" nos mesmos graficos
%
%   Uso tipico p/ grafico log-log de convergencia (estilo log(E) x log(1/h)):
%     para cada malha, guarde max(L1) e max(L2) ao longo dos passos de tempo
%     (ou o valor no passo final), e plote log10(1./h) x log10(L1), log10(L2)
%     junto com uma reta de referencia de inclinacao -2 (2a ordem).

    dS   = elemarea .* (theta_new - theta_old);
    Rloc = dS - dt*( sourcevector - s*flowresult );      % residuo por elemento

    MBE  = sum(Rloc);

    Fb        = sum(flowrate(1:bedgesize));                    % soma so das faces de contorno
    MBE_check = sum(dS) - dt*sum(sourcevector) + s*dt*Fb;       % deve bater com MBE (telescoping de accumarray)

    denom   = max( sum(abs(dS)), dt*(sum(abs(sourcevector)) + abs(Fb)) );
    MBE_rel = abs(MBE) / denom;

    if abs(MBE - MBE_check) > 1e-8*max(1,abs(MBE))
        warning('ferncodes_MBE:inconsistente', ...
            ['MBE (%.3e) e MBE_check (%.3e) nao batem -- verifique indices/sinais ' ...
             'na montagem de flowresult (accumarray).'], MBE, MBE_check);
    end

    % ── normas intensivas (residuo por unidade de area) ──────────────────
    % evita que a queda observada seja so efeito de escala (elemarea~h^2)
    Rdens = Rloc ./ elemarea;                       % residuo "por area", cada elemento
    Atot  = sum(elemarea);

    L1 = sum(abs(Rdens) .* elemarea) / Atot;        % media ponderada por area de |Rdens|
    L2 = sqrt( sum((Rdens.^2) .* elemarea) / Atot );% RMS ponderado por area
end
