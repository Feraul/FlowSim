function salvarStartDat()
    arquivoOrigem = 'C:\Users\flc59\Documents\FlowSim\Start.dat';
    pastaBase     = 'C:\Users\flc59\Documents\FlowSim\benchmarks';

    if ~isfile(arquivoOrigem)
        error('salvarStartDat:naoEncontrado', 'Start.dat nao encontrado em %s', arquivoOrigem);
    end

    conteudo = fileread(arquivoOrigem);

    % procura um caminho tipo "...Benchmonofasico1" ou "...BenchHydraulic439"
    % e extrai os digitos no FINAL do nome (antes de quebra de linha/espaco)
    tok = regexp(conteudo, '([A-Za-z_]+)(\d+)\s*(\r?\n|$)', 'tokens', 'once');

    if isempty(tok)
        error('salvarStartDat:numcaseNaoEncontrado', ...
            'Nao foi possivel extrair o numero do caso a partir do Start.dat.');
    end

    numCaso = str2double(tok{2});
    pastaCaso = fullfile(pastaBase, sprintf('Caso%d', numCaso));

    if ~exist(pastaCaso, 'dir')
        mkdir(pastaCaso);
    end

    destino = fullfile(pastaCaso, 'Start.dat');
    copyfile(arquivoOrigem, destino);

    fprintf('Start.dat (Caso%d) copiado para: %s\n', numCaso, destino);
end