% definimos os pocos para cada do problema

function source_wells= defineWells(env,parmRichardEq)

wells =0;
numcase=env.config.numcase;

if numcase==341
    [P]=ferncodes_calcfonte;
elseif numcase==341.1
    [P]=ferncodes_calcfonte_1D;
else
    P=0;
end
% (2)Catch "source" came from "PLUG_sourcefunction"
source = PLUG_sourcefunction(P,env,0,parmRichardEq);

source_wells.wells=wells;
source_wells.source=source;
end