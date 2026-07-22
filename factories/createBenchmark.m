%--------------------------------------------------------------------------
% createBenchmark — Factory de benchmarks (casos de simulacao)
%
% Instancia e retorna o objeto benchmark correspondente ao numcase
% lido do Start.dat. E o UNICO lugar do simulador que conhece a
% relacao entre numcase e a classe concreta correspondente.
%
% USO NO MAIN:
%   env.benchmark = createBenchmark(env.config.numcase);
%
% COMO ADICIONAR UM NOVO CASO:
%   1. Crie o arquivo benchmarks/CasoNNN.m
%      (herdando de SimulacaoBase e implementando todos os metodos Abstract)
%   2. Adicione UMA linha aqui: case NNN, bench = CasoNNN();
%   3. Nao e necessario modificar nenhum outro arquivo do motor
%
% AGRUPAMENTOS:
%   Varios numcase podem compartilhar o mesmo benchmark quando a fisica
%   e identica e apenas os parametros (dt, h_init, theta_s...) diferem.
%   Nesse caso, use a sintaxe: case {NNN, MMM}, bench = CasoNNN();
%   e trate as diferencas de parametros dentro de preprocessar().
%
% INTERVALOS DE NUMCASE (convencao do simulador):
%   1   – 100  → problemas de referencia (tensores, anisotropia, convergencia)
%   100 – 200  → escoamento monofasico
%   200 – 300  → transporte de contaminantes
%   300 – 400  → carga hidraulica (groundwater)
%   400 – 500  → equacao de Richards (zona saturada/nao-saturada)
%
% NOTA SOBRE numcase COM PONTO DECIMAL (ex: 341.1, 21.1, 34.6):
%   O MATLAB aceita numcase com ponto decimal no switch/case.
%   Use nomes de classe com 'p' no lugar do ponto: Caso341p1, Caso21p1
%   para evitar conflitos de nomenclatura de arquivos no Windows/Linux.
%--------------------------------------------------------------------------
function bench = createBenchmark(numcase)

switch numcase
    % caso monofasico
    case 1, bench=Caso1();

    case 2, bench=Caso2();

    %% ── Carga hidraulica — Groundwater (300–350) ─────────────────

    % Caso 331 — aquifero transiente com termo de armazenamento
    case 331,               bench = Caso331();

    % Caso 333 — aquifero confinado estacionario (solucao de referencia)
    case 333,               bench = Caso333();

    % Casos 330,332,334,335,337,338 — variantes do aquifero confinado
    % compartilham o mesmo benchmark; diferencas em preprocessar()
    case {330,332,334,335,337,338}, bench = Caso330();

    % Caso 336 — aquifero com permeabilidade variavel espacialmente
    case 336,               bench = Caso336();

    % Caso 341 — aquifero nao-confinado com BC mistas (Dirichlet + Neumann)
    % usa ferncodes_K() para permeabilidade variavel na fronteira Neumann
    case 341,               bench = Caso341();

    % Caso 341.1 — variante 1D do caso 341 (usa ferncodes_K_1D)
    case 341.1,             bench = Caso3411();

    % Caso 342 — aquifero transiente com comparacao analitica (erfc)
    case 342,               bench = Caso342();

    % Casos 343, 347, 248 — kmap nao e alterado (no-op em configurarPermeabilidade)
    % benchmark minimalista: apenas inicializa sem modificar kmap ou BC
    case {343, 347, 248},   bench = Caso343();

    %% ── Equacao de Richards (400–500) ────────────────────────────

    % Casos 431, 435 — Van Genuchten com infiltracao vertical
    % compartilham o mesmo modelo de Kr e theta; diferencas em preprocessar()
    %   431: malha z=100, dados de Caviedes (theta_s=0.363, alpha=0.01)
    %   435: variante com condicao de saturacao parcial inicial
    case {431, 435},        bench = Caso431();

    % Casos 432, 433, 434 — Van Genuchten, solos distintos
    % compartilham a mesma formula de Kr; diferencas em preprocessar():
    %   432: Beit Netofa clay, BC de Dirichlet variavel no tempo
    %   433: silt loam com fluxo de Neumann
    %   434: solo com coluna de drenagem (h_init = 65-z)
    case {432, 433, 434},   bench = Caso432();

    % Caso 436 — Van Genuchten com expoente nvg (variante de Moholt)
    % BC dependente do tempo: h = -3t*x*(1-x)*y*(1-y) - 1
    % precisaAtualizarFlags() retorna true → recalcula flags a cada dt
    case 436,               bench = Caso436();

    % Caso 437 — modelo de Gardner (permeabilidade exponencial)
    % K(h) = Ks * exp(alpha*h) — adequado para solos de textura fina
    % tem solucao analitica: calcula erros L2 e H1 em finalizar()
    case 437,               bench = Caso437();

    % Caso 438 — modelo cubico com tensor anisotropico rotacionado
    % K(h) = Ks*(2-h)^(-1/3) — dominio com dois materiais (mask z>0.5)
    case 438,               bench = Caso438();

    % Caso 439 — Brooks-Corey com recarga de aquifero 2D
    % Kr(h) = c/(c+|h|^5) — dominio 300x200m, 6 pontos de monitoramento
    % salva indices em .mat para reutilizacao entre simulacoes
    case 439,               bench = Caso439();

    %% ── Problemas de referencia com tensores (1–100) ─────────────

    % Caso 21.1 — tensor anisotropico rotacionado (theta=5*pi/12)
    % K variavel: K11=1+2x²+y², K22=1+x²+2y²
    case 21.1,              bench = Caso21p1();

    % Caso 34.6 — campo aleatorio de permeabilidade (log-normal)
    % usa getrandist() para gerar distribuicao espacial estocastica
    case 34.6,              bench = Caso346();

    % Caso 34.7 — permeabilidade obtida de mapa de cores (getchue)
    case 34.7,              bench = Caso347c();

    % Caso 35 — permeabilidade periodica: K=exp(sqrt(s)*cos(2pi*m*x)*cos(2pi*l*y))
    case 35,                bench = Caso35();

    % Caso 36 — permeabilidade em camadas (K=1,2,5,10 por faixas de y)
    case 36,                bench = Caso36();

    %% ── Transporte de contaminantes (200–300) ────────────────────

    % Caso 241 — tensor anisotropico fixo (k=[10,0;0,0.01] rotacionado pi/7.2)
    case 241,               bench = Caso241();

    % Caso 245 — permeabilidade periodica com s=4 (variante do caso 35)
    case 245,               bench = Caso245();

    % Casos 247, 249, 250 — campo de permeabilidade lido de arquivo externo
    % (Perm_Var0p1.mat) — kmap carregado e flipado em preprocessar()
    case {247, 249, 250},   bench = Caso247();

    %% ── Caso nao registrado ──────────────────────────────────────
    % Se chegar aqui, o numcase do Start.dat nao tem benchmark associado.
    % Verifique:
    %   1. O numcase no Start.dat esta correto?
    %   2. O arquivo CasoNNN.m existe na pasta benchmarks/?
    %   3. Foi adicionada a linha correspondente acima?
    otherwise
        error(['createBenchmark: numcase %g nao registrado.\n' ...
               'Para adicionar: crie benchmarks/CasoNNN.m e adicione ' ...
               'uma linha nesta funcao.'], numcase);
end
end