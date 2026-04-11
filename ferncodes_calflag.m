%It is called by "preMPFA.m"

function [nflag,nflagface] = ferncodes_calflag(env,time)
    
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
    x_vertex = env.config.bcflag(:,1) == vertex_flag';   % matriz lógica (#bcflag × #faces)
    x_face   = env.config.bcflag(:,1) == face_flag';     % idem

    % Para cada face, encontrar a linha correta de bcflag
    [~, bc_row_vertex] = max(x_vertex, [], 1);   % 1×#faces
    [~, bc_row_face]   = max(x_face,   [], 1);

    % Extrair flags e valores
    bcflag_vertex = env.config.bcflag(bc_row_vertex, :);
    bcflag_face   = env.config.bcflag(bc_row_face, :);

    % ============================================================
    % CASO 341 e 341.1  (mesma lógica do seu código)
    % ============================================================
    if env.config.numcase == 341 || env.config.numcase == 341.1

        % nflag para vértices
        nflag(vertex_idx,1) = bcflag_vertex(:,1);
        nflag(vertex_idx,2) = arrayfun(@(v,r) PLUG_bcfunction(v,r,time,env), ...
                                       vertex_idx, bc_row_vertex');

        % nflagface para faces
        nflagface(:,1) = bcflag_face(:,1);
        nflagface(:,2) = arrayfun(@(v1,v2,r) PLUG_bcfunction([v1 v2],r,time,env), ...
                                  face_idx(:,1), face_idx(:,2), bc_row_face');

        return
    end

    % ============================================================
    % CASOS 432 e 434 (têm lógica especial)
    % ============================================================
    if env.config.numcase == 432 || env.config.numcase == 434

        special = bcflag_vertex(:,1) == 101;

        % Caso especial
        nflag(vertex_idx(special),1) = bcflag_vertex(special,1);
        nflag(vertex_idx(special),2) = arrayfun(@(v,r) PLUG_bcfunction(v,r,time,env), ...
                                                vertex_idx(special), bc_row_vertex(special)');

        % Caso normal
        normal = ~special;
        nflag(vertex_idx(normal),1) = bcflag_vertex(normal,1);
        nflag(vertex_idx(normal),2) = bcflag_vertex(normal,2);

        return
    end

    % ============================================================
    % CASO PADRÃO (todos os outros)
    % ============================================================
    nflag(vertex_idx,1) = bcflag_vertex(:,1);
    nflag(vertex_idx,2) = arrayfun(@(v,r) PLUG_bcfunction(v,r,time,env), vertex_idx, bc_row_vertex');
    nflagface(:,1)=bcflag_face(:,1);
    nflagface(:,2) = arrayfun(@(v1,v2,r) PLUG_bcfunction([v1 v2],r,time,env), ...
                                  face_idx(:,1), face_idx(:,2), bc_row_face');

end