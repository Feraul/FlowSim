function [mvector]=sourceterm(mvector,source_wells)

   if isempty(source_wells) || (isstruct(source_wells) && ...
            (~isfield(source_wells,'source') || isempty(source_wells.source)))
        return;
    end
    source=source_wells.source;
    %The vector "mvector" is added to vector returned from function below.
    mvector = mvector + source;

end 