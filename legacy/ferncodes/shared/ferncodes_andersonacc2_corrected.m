function [x, erro, iter, tabletol] = ferncodes_andersonacc2_corrected(x, auxtol, ...
    R0, env, parmRichardEq, ...
    preMPFAD, dt, time, source_wells)
% [x, erro, iter, tabletol] = FERNCODES_ANDERSONACC2_CORRECTED(x, auxtol, R0, env, ...
%    parmRichardEq, preMPFAD, dt, time, source_wells)
%
% Método de aceleração de Anderson para resolver sistema não linear F(x)=0,
% onde F(x) = M(x)*x - RHS(x). A cada iteração resolve-se M(x)*g = RHS(x)
% e aplica-se aceleração de Anderson com janela deslizante.
%
% Entradas:
%   x             - vetor coluna inicial (estimativa)
%   auxtol        - tolerância absoluta e relativa (escalar)
%   R0            - norma do resíduo inicial (para erro relativo)
%   env           - estrutura com configurações (deve conter env.config.pmethod,
%                   env.config.nltol, env.config.aa_start opcional)
%   parmRichardEq - parâmetros da equação de Richards
%   preMPFAD      - estrutura para método MPFA
%   dt, time      - passo de tempo e tempo atual
%   source_wells  - termos fonte/pogos
%
% Saídas:
%   x      - solução aproximada
%   erro   - erro relativo final (norma do resíduo / R0, ou norma do resíduo se R0=0)
%   iter   - número de iterações realizadas
%   tabletol - histórico [iter, erro] (opcional)

    % Parâmetros do método
    pmethod = env.config.pmethod;
    nltol   = env.config.nltol;
    
    % Parâmetros de aceleração de Anderson (podem ser definidos em env.config)
    mMax     = min(9, size(x,1));            % tamanho máximo da janela
    itmax    = 10000;                        % número máximo de iterações
    atol     = auxtol;                       % tolerância absoluta
    rtol     = auxtol;                       % tolerância relativa
    beta     = 0.6;                          % parâmetro de relaxamento (damping)
    % Início da aceleração: usar env.config.aa_start se existir, senão 25
   
    AAstart = 7;
    
    % Inicialização
    iter = 0;
    erro = 1.0;
    
    % Matrizes para armazenamento do histórico (janela deslizante)
    G = [];      % armazena g (soluções)
    F = [];      % armazena resíduos f = g - x
    
    % Histórico de convergência
    tabletol = zeros(itmax+2, 2);
    
    % -------------------------------------------------------------
    % Primeira avaliação do sistema (iteração 0)
    % -------------------------------------------------------------
    % Nota: evaluate_system retorna também um x, mas ignoramos para não sobrescrever
    [M_new, RHS_new, ~] = evaluate_system(x, pmethod, env, parmRichardEq, ...
        preMPFAD, dt, time, source_wells);
    
    % Resíduo original: ||M*x - RHS||
    RR = norm(M_new * x - RHS_new);
    if R0 == 0.0
        erro = RR;                  % erro absoluto
    else
        erro = abs(RR / R0);        % erro relativo
    end
    tabletol(1, :) = [0, erro];
    
    % Se o erro inicial já é aceitável, retorna imediatamente
    if erro < nltol
        iter = 0;
        return;
    end
    M_old=M_new;
    RHS_old=RHS_new;
    
    % -------------------------------------------------------------
    % Loop principal de iterações
    % -------------------------------------------------------------
    for iter = 1:itmax

        % --- bloco de estabilização L-scheme ---
        n = size(M_old,1);
        M_L   = M_old + 1 * speye(n);     % M^{(k)} + L I
        RHS_L = RHS_old + 1 * parmRichardEq.h_old;      % RHS^{(k)} + L h^{(k)}
        % ---------------------------------------------------------
        % 1. Resolver sistema linear para obter g (ponto fixo)
        % ---------------------------------------------------------
        gval = solve_system(M_L, RHS_L);
        
       
        % Resíduo do ponto fixo (usado apenas na aceleração, não no critério de parada)
        fval = gval - x;
        
        % ---------------------------------------------------------
        % 2. Aceleração de Anderson (se ativa)
        % ---------------------------------------------------------
        if mMax > 0 && iter >= AAstart
            % Armazenar valores atuais na janela
            if isempty(G)
                G = gval;
                F = fval;
                mAA = 1;
            else
                % Atualizar janela deslizante
                if size(G, 2) < mMax
                    G = [G, gval];
                    F = [F, fval];
                else
                    G = [G(:, 2:end), gval];
                    F = [F(:, 2:end), fval];
                end
                mAA = size(G, 2);
            end
            
            % Aplicar aceleração se houver pelo menos dois pontos
            if mAA >= 2
                % Matriz de diferenças dos resíduos
                dF = diff(F, 1, 2);
                
                % Resolver problema de mínimos quadrados:
                % encontrar alpha que minimiza ||f_k - dF * gamma||,
                % com alpha = [1; -gamma]
                [Q, R] = qr(dF, 0);
                gamma = R \ (Q' * fval);
                alpha = [1; -gamma];
                alpha = alpha / sum(alpha);   % normaliza para soma = 1
                
                % Nova aproximação combinada
                x_new = G * alpha;
                
                % Aplicar damping (relaxamento)
                % x = (1 - beta) * x_old + beta * x_new
                x = (1 - beta) * x + beta * x_new;
            else
                % Ainda não há histórico suficiente para aceleração
                x = gval;
            end
        else
            % Sem aceleração (iterações iniciais)
            x = gval;
        end
        parmRichardEq.hold=x;
        % ---------------------------------------------------------
        % 3. Avaliar o sistema no novo ponto x
        % ---------------------------------------------------------
        [M_new, RHS_new, ~] = evaluate_system(x, pmethod, env, parmRichardEq, ...
            preMPFAD, dt, time, source_wells);
        
        % Calcular resíduo original e erro
        RR = norm(M_new * x - RHS_new);
        if R0 == 0.0
            erro = RR;
        else
            erro = abs(RR / R0);
        end
        tabletol(iter+1, :) = [iter, erro];
        
        % ---------------------------------------------------------
        % 4. Verificar convergência (critério único: erro < nltol)
        % ---------------------------------------------------------
        if erro < nltol
            break;
        end
    end % for iter
    
    M_old=M_new;
    RHS_old=RHS_new;
    % Truncar histórico
    tabletol = tabletol(1:iter+1, :);
end

% ---------------------------------------------------------------------
% Função auxiliar: solve_system
% Resolve M * g = RHS usando GMRES com pré-condicionador ILU adaptativo
% ---------------------------------------------------------------------
function gval = solve_system(M, RHS)
    n = size(M, 1);
    restart = min(100, n);
    maxit = min(1000, n);
    tol = 1e-6;
    
    % Estimar condicionamento para decidir estratégia de ILU
    % (evita adicionar perturbação diagonal desnecessária)
    try
        c = condest(M);
        ill_conditioned = (c > 1e12);
    catch
        % Se condest falhar (matriz muito grande), assume mal condicionada
        ill_conditioned = true;
    end
    
    if ill_conditioned
        % Matriz mal condicionada: usar ILUTP com perturbação diagonal
        setup.type = 'ilutp';
        setup.droptol = 1e-6;
        setup.udiag = true;
        M_reg = M + 1e-8 * speye(n);
        [L, U] = ilu(M_reg, setup);
    else
        % Matriz bem condicionada: usar ILU(0) sem perturbação
        setup.type = 'nofill';
        [L, U] = ilu(M, setup);
    end
     
    % Resolver com GMRES
    [gval, flag, ~, ~] = gmres(M, RHS, restart, tol, maxit, L, U);
    
    if flag ~= 0
        % GMRES falhou: tentar solução direta para sistemas pequenos
        if n < 1000
            warning('GMRES não convergiu (flag=%d), usando solução direta.', flag);
            gval = M \ RHS;
        else
            warning('GMRES não convergiu (flag=%d) e sistema é grande. Solução pode ser imprecisa.', flag);
        end
    end
end

% ---------------------------------------------------------------------
% Função auxiliar: evaluate_system
% Avalia a matriz M e o vetor RHS no ponto x, conforme o método escolhido.
% Nota: O terceiro argumento de saída (x) é mantido por compatibilidade,
%       mas não é utilizado no código principal (descartado com '~').
% ---------------------------------------------------------------------
function [M_new, RHS_new, x] = evaluate_system(x, pmethod, env, parmRichardEq, ...
    preMPFAD, dt, time, source_wells)
    
    % Esta função deve ser implementada de acordo com o simulador específico.
    % As variáveis abaixo são meramente ilustrativas e DEVEM ser substituídas
    % pelos campos reais existentes em env, parmRichardEq, etc.
    %
    % Exemplo de implementação para cada método:
    
    switch pmethod
        case 'nlfvpp'
            % --- Exemplo NLFVPP (substituir pelas variáveis corretas) ---
            % As variáveis nflag, w, s, Con, etc. precisam estar definidas
            error('Método NLFVPP não implementado: faltam variáveis como nflag, w, s, Con, ...');
            
        case 'nlfvh'
            % --- Exemplo NLFVH ---
            error('Método NLFVH não implementado: faltam variáveis como parameter, mobility, etc.');
            
        case 'nlfvdmp'
            % --- Exemplo NLFVDMP ---
            error('Método NLFVDMP não implementado: faltam variáveis como weightDMP, mobility, etc.');
            
        case 'mpfad'
            % --- Exemplo MPFAD (baseado no código original) ---
            % Nota: esta implementação modifica parmRichardEq.h_old, o que pode ser intencional
            parmRichardEq.h_old = x;
            [env, parmRichardEq] = PLUG_kfunction(env, parmRichardEq, time);
            [preMPFAD] = ferncodes_Kde_Ded_Kt_Kn(env, parmRichardEq, preMPFAD);
            [preMPFAD, ~, ~] = ferncodes_Pre_LPEW_2_vect(preMPFAD, parmRichardEq, env);
            [M, I, ~] = ferncodes_globalmatrix(env, preMPFAD, parmRichardEq);
            [M_new, I] = addsource(sparse(M), I, source_wells, env);
            RHS_new = sourceterm(I, source_wells);
            
        otherwise
            % Método padrão (exemplo com NLFVH)
            error('Método desconhecido: %s', pmethod);
    end
    
    % Nota: o terceiro argumento de saída (x) é retornado inalterado para compatibilidade.
    % No código principal, este valor é ignorado (descartado com '~').
end