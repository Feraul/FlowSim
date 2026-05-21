%It is called by "preMPFA.m"

function [env] = ferncodes_calflag(env,parmRichardEq,time)
    m=1;
    
    nelem_nodes = size(env.geometry.coord,1);
    nelem_faces = size(env.geometry.bedge,1);

    % Inicializa
    nflag      = 5000 * ones(nelem_nodes,2);
    nflagface  = zeros(nelem_faces,2);

    % Flags de vértices e faces
    vertex_flag = env.geometry.bedge(:,4);
    face_flag   = env.geometry.bedge(:,5);

    % Índices dos vértices
    vertex_idx  = env.geometry.bedge(:,1);
    face_idx    = env.geometry.bedge(:,1:2);

    % Mapeamento lógico para vértices
    x_vertex = env.config.bcflag(:,1) == vertex_flag';   % (#bcflag × #faces)
    x_face   = env.config.bcflag(:,1) == face_flag';

    % Para cada face, encontrar a linha correta de bcflag
    [~, bc_row_vertex] = max(x_vertex, [], 1);   % 1×#faces
    [~, bc_row_face]   = max(x_face,   [], 1);

    % Extrair flags e valores
    bcflag_vertex = env.config.bcflag(bc_row_vertex, :);
    bcflag_face   = env.config.bcflag(bc_row_face, :);

    % ============================================================
    % CASO 341 e 341.1
    % ============================================================
    if env.config.numcase == 341 || env.config.numcase == 341.1

        % nflag para vértices
        nflag(vertex_idx,1) = bcflag_vertex(:,1);

        % *** VETORIZADO ***
        nflag(vertex_idx,2) = PLUG_bcfunction(vertex_idx, bc_row_vertex', time, env);

        % nflagface para faces
        nflagface(:,1) = bcflag_face(:,1);

        % *** VETORIZADO PARA FACES ***
        v1 = face_idx(:,1);
        v2 = face_idx(:,2);
        nflagface(:,2) = PLUG_bcfunction([v1 v2], bc_row_face', time, env);

        return
    end

    % ============================================================
    % CASOS 432 e 434
    % ============================================================
    if env.config.numcase == 432 || env.config.numcase == 434

        special = bcflag_vertex(:,1) == 101;
        normal  = ~special;

        % Caso especial
        nflag(vertex_idx(special),1) = bcflag_vertex(special,1);
        nflag(vertex_idx(special),2) = PLUG_bcfunction(vertex_idx(special), bc_row_vertex(special)', time, env);

        % Caso normal
        nflag(vertex_idx(normal),1) = bcflag_vertex(normal,1);
        nflag(vertex_idx(normal),2) = bcflag_vertex(normal,2);

        return
    end

    % ============================================================
    % CASO PADRÃO
    % ============================================================

    nflag(vertex_idx,1) = bcflag_vertex(:,1);
     mmmm= PLUG_bcfunction(vertex_idx, bc_row_vertex', time, env,parmRichardEq);
    nflag(vertex_idx,2)=mmmm(vertex_idx);
    nflagface(:,1) = bcflag_face(:,1);
    nflagface(:,2) = PLUG_bcfunction([face_idx(:,1) face_idx(:,2)], ...
        bc_row_face', time, env,parmRichardEq);

    env.config.nflag=nflag;
    env.config.nflagface=nflagface;

end
