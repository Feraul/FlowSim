%% plot_convergencia_MBE.m
% Reproduz o estilo do grafico log(E) x log(1/refino) com reta de
% referencia de ordem teorica, no mesmo padrao usado para o erro L1/L2
% de pressao (MPFA-D vs TPFA). Aqui "E" e a norma da MBE (L1_MBE,
% L2_MBE) -- preencha os vetores com os valores impressos por
% Caso436.finalizar ("L1_MBE maximo" / "L2_MBE maximo") para cada rodada.
%
% Dois modos:
%   modo = 'dx' -> estudo de convergencia ESPACIAL (dt fixo/pequeno ou
%                  casado dt=C*dx^2). Eixo x = log(1/dx). Reta de
%                  referencia de 2a ordem (MPFA-D).
%   modo = 'dt' -> estudo de convergencia TEMPORAL (dx fixo/pequeno ou
%                  casado dx=dt/2, mantendo dx^2/dt pequeno). Eixo x =
%                  log(1/dt). Reta de referencia de 1a ordem (Euler
%                  implicito).

modo = 'dx';   % 'dx' ou 'dt'

switch modo
    case 'dx'
        refino   = [1/8, 1/16, 1/32, 1/64, 1/128];   % dx de cada rodada
        ordemRef = 2;
        xlab     = 'log(1/dx)';
        tit      = 'Convergencia espacial da MBE (dt casado, dt=0.5\,dx^2)';
    case 'dt'
        refino   = [1/16, 1/32, 1/64];                % dt de cada rodada
        ordemRef = 1;
        xlab     = 'log(1/dt)';
        tit      = 'Convergencia temporal da MBE (dx casado, dx=dt/2)';
    otherwise
        error('modo deve ser ''dx'' ou ''dt''');
end

% ── preencha com os valores reais (max(extras.L1_MBE), max(extras.L2_MBE)) ──
L1_MPFAD = nan(size(refino));
L2_MPFAD = nan(size(refino));

% Se tiver rodado TPFA tambem para comparacao (mesmo estilo da figura):
L1_TPFA = nan(size(refino));
L2_TPFA = nan(size(refino));

logx = log10(1./refino);

figure('Color','w'); hold on

plot(logx, log10(L1_MPFAD), '-o', 'Color','k', 'MarkerFaceColor','k', ...
    'LineWidth', 1.3, 'DisplayName', 'MPFA-D L_1');
plot(logx, log10(L2_MPFAD), '-s', 'Color','k', 'MarkerFaceColor','w', ...
    'LineWidth', 1.3, 'DisplayName', 'MPFA-D L_2');

if any(~isnan(L1_TPFA))
    plot(logx, log10(L1_TPFA), '-o', 'Color',[0.6 0.6 0.6], 'MarkerFaceColor',[0.6 0.6 0.6], ...
        'LineWidth', 1.1, 'DisplayName', 'TPFA L_1');
    plot(logx, log10(L2_TPFA), '-s', 'Color',[0.6 0.6 0.6], 'MarkerFaceColor','w', ...
        'LineWidth', 1.1, 'DisplayName', 'TPFA L_2');
end

% ── reta de referencia (ordem 2 p/ 'dx', ordem 1 p/ 'dt'), ancorada no
%    primeiro ponto valido ──────────────────────────────────────────
valid = ~isnan(L1_MPFAD);
i0 = find(valid, 1, 'first');
if ~isempty(i0)
    refLine = log10(L1_MPFAD(i0)) - ordemRef*(logx - logx(i0));
    plot(logx, refLine, '--k', 'LineWidth', 1.0, 'DisplayName', ...
        sprintf('%dth order', ordemRef));

    % ordem observada via ajuste de minimos quadrados (log-log)
    p1 = polyfit(logx(valid), log10(L1_MPFAD(valid)), 1);
    p2 = polyfit(logx(valid), log10(L2_MPFAD(valid)), 1);
    fprintf('[%s] Ordem observada (ajuste log-log) -- L1: %.3f | L2: %.3f (teorica: %d)\n', ...
        modo, p1(1), p2(1), ordemRef);
end

xlabel(xlab);
ylabel('log(E)');
legend('Location','southwest');
grid on
title(tit);
