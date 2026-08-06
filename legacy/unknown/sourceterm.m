function [mvector] = sourceterm(mvector, source_wells)
    if isstruct(source_wells) && isfield(source_wells,'source') &&...
            ~isempty(source_wells.source)
        source = source_wells.source;
    else
        % Caso439/Vauclin (ou qualquer caso sem fonte): source_wells = false, [], ou struct sem 'source'
        source = zeros(size(mvector));
    end
    mvector = mvector + source;
end