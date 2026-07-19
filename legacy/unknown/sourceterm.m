function [mvector]=sourceterm(mvector,source_wells)

   if isempty(source_wells) || (isstruct(source_wells) && ...
            (~isfield(source_wells,'wells') || isempty(source_wells.wells)))
        return;
    end
    source=source_wells.source;
    %The vector "mvector" is added to vector returned from function below.
    mvector = mvector + source;

end 