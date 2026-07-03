%--------------------------------------------------------------------------
% createSimulacao — Factory de tipos de simulacao
%
% Instancia e retorna o objeto sim correspondente ao valor de phasekey
% lido do Start.dat. E o UNICO lugar do simulador que conhece a relacao
% entre o numero phasekey e a classe concreta correspondente.
%
% USO NO MAIN:
%   sim = createSimulacao(env.config.phasekey);
%
% DIFERENCA ENTRE sim E env.benchmark:
%
%   sim (criado aqui) → define o TIPO de simulacao (phasekey)
%     Responsabilidades:
%       - preprocessar(): chama preRE, prehydraulic, etc.
%       - definirFontes(): chama defineWells
%       - configurarFlags(): monta nflag e nflagface (generico)
%     Exemplo: SimRichards sabe que precisa chamar preRE() e defineWells()
%              mas NAO sabe se e o caso 439 ou 437
%
%   env.benchmark (criado em createBenchmark) → define o CASO especifico
%     Responsabilidades:
%       - preprocessar(): preenche theta_s, h_init, dt, h_old para AQUELE caso
%       - configurarPermeabilidade(): K(h) especifico do modelo
%       - calcularTheta(), calcularCapacidade(): modelo hidrico especifico
%       - inicializar(), atualizarEstado(), finalizar(): logica do loop temporal
%     Exemplo: Caso439 sabe que usa Brooks-Corey, tem 6 pontos de monitoramento,
%              dominio 300x200m, etc.
%
% COMO ADICIONAR UM NOVO TIPO DE SIMULACAO:
%   1. Crie o arquivo simulacoes/SimXXX.m
%      (herdando de SimulacaoBase e implementando todos os metodos Abstract)
%   2. Defina um novo valor de phasekey no Start.dat (ex: phasekey = 7)
%   3. Adicione UMA linha aqui: case 7, sim = SimXXX();
%   4. Nao e necessario modificar nenhum outro arquivo do motor
%
% TIPOS DISPONIVEIS (phasekey no Start.dat):
%
%   [0] SimHiperbolica      — apenas equacao hiperbolica (transporte puro)
%                             sem equacao de pressao / carga hidraulica
%
%   [1] SimMonofasica       — escoamento monofasico (1 fluido, lei de Darcy)
%                             resolve apenas a equacao de pressao
%                             casos de referencia: 2, 3, 4.1, 5.1, 11, 12, 13...
%
%   [2] SimBifasica         — escoamento bifasico (2 fluidos imisc., IMPES)
%                             resolve pressao (implicito) + saturacao (explicito)
%                             casos tipicos: agua-oleo em reservatorios
%
%   [3] SimContaminante     — transporte de contaminante (1 fase + concentracao)
%                             resolve pressao + equacao de adveccao-difusao
%                             casos: 241, 245, 247, 248, 249, 250, 251
%
%   [4] SimGroundwater      — carga hidraulica (lei de Darcy com armazenamento)
%                             resolve equacao de groundwater (Darcy + armazenamento)
%                             casos: 331, 333, 334, 335, 336, 337, 338, 341, 342, 343...
%
%   [5] SimContamGroundwater — carga hidraulica + transporte de contaminante
%                              acoplamento entre equacao de groundwater e
%                              equacao de adveccao-difusao (IMHEC)
%
%   [6] SimRichards         — equacao de Richards (zona saturada/nao-saturada)
%                             nao-linear: K = K(h), theta = theta(h)
%                             resolve por iteracoes de Picard a cada dt
%                             casos: 431, 432, 433, 434, 435, 436, 437, 438, 439
%--------------------------------------------------------------------------
function sim = createSimulacao(tipoID)

switch tipoID

    % ── [0] Apenas equacao hiperbolica ───────────────────────────
    % Transporte puro sem equacao de pressao.
    % Ficheiro: simulacoes/SimHiperbolica.m
    case 0
        sim = SimHiperbolica();

    % ── [1] Escoamento monofasico (One-phase flow) ────────────────
    % Um unico fluido, lei de Darcy estacionaria ou quasi-estacionaria.
    % O solver resolve apenas a equacao de pressao p (ou h).
    % E o tipo mais comum para problemas de referencia e validacao.
    % Ficheiro: simulacoes/SimMonofasica.m
    case 1
        sim = SimMonofasica();

    % ── [2] Escoamento bifasico (Two-phase flow) ──────────────────
    % Dois fluidos imisciveis (ex: agua e oleo).
    % Esquema IMPES: pressao implicita + saturacao explicita.
    % Ficheiro: simulacoes/SimBifasica.m
    case 2
        sim = SimBifasica();

    % ── [3] Transporte de contaminante ───────────────────────────
    % Um fluido com soluto transportado (adveccao + difusao).
    % Resolve equacao de pressao acoplada a equacao de concentracao C.
    % Ficheiro: simulacoes/SimContaminante.m
    case 3
        sim = SimContaminante();

    % ── [4] Carga hidraulica (Groundwater) ────────────────────────
    % Lei de Darcy com armazenamento especifico (zona saturada).
    % Pode ser estacionario (numcase 333, 336...) ou transiente (342...).
    % Ficheiro: simulacoes/SimGroundwater.m
    case 4
        sim = SimGroundwater();

    % ── [5] Contaminante + carga hidraulica acoplados ─────────────
    % Resolve simultaneamente a carga hidraulica h e a concentracao C.
    % Usa o metodo IMHEC (Implicit Method for Hydraulic head
    % and Contaminant transport).
    % Ficheiro: simulacoes/SimContamGroundwater.m
    case 5
        sim = SimContamGroundwater();

    % ── [6] Equacao de Richards ───────────────────────────────────
    % Fluxo em meio poroso parcialmente saturado (zona vadosa).
    % Nao-linear: K = K(h) e theta = theta(h) dependem da carga h.
    % Loop temporal com iteracoes de Picard a cada passo de tempo.
    % A fisica especifica (modelo de K e theta) e definida pelo benchmark:
    %   Van Genuchten (431-435), Gardner (437), Brooks-Corey (439)...
    % Ficheiro: simulacoes/SimRichards.m
    case 6
        sim = SimRichards();

    % ── Tipo nao reconhecido ──────────────────────────────────────
    % Se chegar aqui, o valor de phasekey no Start.dat nao tem
    % uma classe correspondente registrada.
    % Verifique:
    %   1. O valor de phasekey no Start.dat esta correto? (0 a 6)
    %   2. O arquivo SimXXX.m existe na pasta simulacoes/?
    %   3. Foi adicionada a linha correspondente acima?
    otherwise
        error(['createSimulacao: phasekey %d nao reconhecido.\n' ...
               'Valores validos: 0 (hiperbolica), 1 (monofasica), ' ...
               '2 (bifasica), 3 (contaminante),\n' ...
               '4 (groundwater), 5 (contam+groundwater), 6 (Richards)\n' ...
               'Para adicionar: crie simulacoes/SimXXX.m e adicione ' ...
               'uma linha nesta funcao.'], tipoID);
end
end