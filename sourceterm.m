function [mvector]=sourceterm(mvector,source_wells)
    source=source_wells.source;
    %The vector "mvector" is added to vector returned from function below.
    mvector = mvector + source;

end 