%--------------------------------------------------------------------------
% createMetodo — Factory de metodos numericos
%
% Instancia e retorna o objeto metodo correspondente ao valor de pmethod
% lido do Start.dat. E o UNICO lugar do simulador que conhece a relacao
% entre a string pmethod e a classe concreta correspondente.
%
% USO NO MAIN:
%   env.metodo = createMetodo(env.config.pmethod);
%
% COMO ADICIONAR UM NOVO METODO NUMERICO:
%   1. Crie o arquivo metodos/MetodoXXX.m
%      (herdando de MetodoBase e implementando todos os metodos Abstract:
%       preprocessar, atualizarPremethod, montarSistema, resolver, calcularFlowrate)
%   2. Adicione UMA linha aqui: case 'xxx', metodo = MetodoXXX();
%   3. No Start.dat, use pmethod = xxx
%   4. Nao e necessario modificar nenhum outro arquivo do motor
%
% METODOS DISPONIVEIS:
%
%   'tpfa'   → MetodoTPFA
%              Two-Point Flux Approximation
%              O mais simples e rapido. Consistente apenas em malhas
%              ortogonais ou com tensores K alinhados com a malha.
%              Sem interpolacao nos vertices. Sem Ded, Kt, pesos LPEW2.
%              Recomendado para: testes rapidos, malhas ortogonais,
%              problemas com K isotropico.
%
%   'mpfad'  → MetodoMPFAD
%              Multi-Point Flux Approximation — Diamond patch (Contreras et al.)
%              Consistente para K heterogeneo e anisotropico em malhas
%              nao-estruturadas. Usa interpolacao LPEW2 nos vertices internos
%              e correcao de anisotropia Ded nas arestas internas.
%              Recomendado para: a maioria dos casos Richards, tensores
%              anisotropicos, malhas triangulares ou quadrilaterais distorcidas.
%
%   'mpfah'  → MetodoMPFAH
%              Multi-Point Flux Approximation — Harmonic points
%              Variante do MPFA que usa pontos harmonicos para interpolacao.
%              Alternativa ao MPFA-D para alguns tipos de malha.
%
%   'nlfvpp' → MetodoNLFVPP
%              Non-Linear Finite Volume — Positive Preserving (Contreras et al. 2021)
%              Metodo nao-linear que garante principio do maximo discreto
%              (sem oscilacoes espurias). Mais caro computacionalmente
%              pois requer iteracoes extras para os pesos nao-lineares.
%              Recomendado para: problemas com alta anisotropia onde
%              MPFA-D apresenta valores negativos indesejaveis.
%
%   'mpfaql' → MetodoMPFAQL
%              Multi-Point Flux Approximation — QL scheme (Contreras et al. 2019)
%              Variante do MPFA que usa pesos baseados em quadrilateros locais.
%
% NOTA TECNICA:
%   lower(char(pmethod)) garante que o valor lido do Start.dat funcione
%   independente de maiusculas/minusculas:
%     'TPFA', 'tpfa', 'Tpfa' → todos instanciam MetodoTPFA
%   Isso evita erros silenciosos por diferenca de caixa na configuracao.
%--------------------------------------------------------------------------
function metodo = createMetodo(pmethod)

% converte para minusculo e char para robustez
% (pmethod pode chegar como string ou char array do Start.dat)
switch lower(char(pmethod))

    % ── TPFA: Two-Point Flux Approximation ───────────────────────
    % Mais simples e eficiente. Use para malhas ortogonais ou K isotropico.
    % Ficheiro: metodos/MetodoTPFA.m
    case 'tpfa'
        metodo = MetodoTPFA();

    % ── MPFA-D: Multi-Point — Diamond patch ──────────────────────
    % Metodo padrao para a maioria dos casos do simulador.
    % Consistente para K anisotropico em malhas nao-estruturadas.
    % Ficheiro: metodos/MetodoMPFAD.m
    case 'mpfad'
        metodo = MetodoMPFAD();

    % ── MPFA-H: Multi-Point — Harmonic points ────────────────────
    % Alternativa ao MPFA-D com interpolacao em pontos harmonicos.
    % Ficheiro: metodos/MetodoMPFAH.m
    case 'mpfah'
        metodo = MetodoMPFAH();

    % ── NL-TPFA / NLFVPP: Non-Linear Positive Preserving ─────────
    % Garante principio do maximo discreto. Mais caro que MPFA-D.
    % Recomendado quando MPFA-D produz valores negativos indesejaveis.
    % Ficheiro: metodos/MetodoNLFVPP.m
    case 'nlfvpp'
        metodo = MetodoNLFVPP();

    % ── MPFA-QL: Multi-Point — QL scheme ─────────────────────────
    % Variante com pesos baseados em quadrilateros locais.
    % Ficheiro: metodos/MetodoMPFAQL.m
    case 'mpfaql'
        metodo = MetodoMPFAQL();

    % ── Metodo nao reconhecido ────────────────────────────────────
    % Se chegar aqui, o valor de pmethod no Start.dat nao tem
    % uma classe correspondente registrada.
    % Verifique:
    %   1. O valor de pmethod no Start.dat esta correto?
    %      Valores validos: tpfa, mpfad, mpfah, nlfvpp, mpfaql
    %   2. O arquivo MetodoXXX.m existe na pasta metodos/?
    %   3. Foi adicionada a linha correspondente acima?
    otherwise
        error(['createMetodo: pmethod "%s" nao reconhecido.\n' ...
               'Valores validos: tpfa, mpfad, mpfah, nlfvpp, mpfaql\n' ...
               'Para adicionar: crie metodos/MetodoXXX.m e adicione ' ...
               'uma linha nesta funcao.'], char(pmethod));
end
end