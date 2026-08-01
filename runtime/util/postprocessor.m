function postprocessor(env, step, time, options)
arguments
    env
    step
    time
    options.pressure         = []
    options.watersaturation  = []
    options.theta_n          = []
    options.k_theta          = []

    options.exact_sol        = []
    options.theta_exact      = []
    options.k_theta_exact    = []
    options.normk            = []
    options.satLabel         = ''   % permite sobrescrever o nome do campo (ex: 'Theta')
end

coord     = env.geometry.coord;
elem      = env.geometry.elem;
filepath  = env.mainpathfolders.path;
resfolder = env.mainpathfolders.resfolder;
auxnumcase = env.config.numcase;
if size(coord,2) == 2
    coord = [coord, zeros(size(coord,1),1)];
end
nelem   = size(elem,1);
npoints = size(coord,1);
triMask  = elem(:,4) == 0;
quadMask = ~triMask;
fname_vtk = fullfile(filepath, resfolder, ['res_00' num2str(step) '.vtk']);
fid = fopen(fname_vtk,'w');

fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'time=%g step=%d\n', time, step);
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

fprintf(fid, 'POINTS %d float\n', npoints);
fprintf(fid, '%.6E %.6E %.6E\n', coord');

cellsize = sum(3*triMask + 4*quadMask) + nelem;
fprintf(fid, 'CELLS %d %d\n', nelem, cellsize);
buf = cell(nelem,1);
if any(triMask)
    tri = elem(triMask,1:3) - 1;
    buf(triMask) = cellstr(num2str(tri, '3 %d %d %d'));
end
if any(quadMask)
    quad = elem(quadMask,1:4) - 1;
    buf(quadMask) = cellstr(num2str(quad, '4 %d %d %d %d'));
end
fprintf(fid, '%s\n', buf{:});

fprintf(fid, 'CELL_TYPES %d\n', nelem);
celltypes = zeros(nelem,1);
celltypes(triMask)  = 5;
celltypes(quadMask) = 9;
fprintf(fid, '%d\n', celltypes);

fprintf(fid, 'CELL_DATA %d\n', nelem);

if ~isempty(options.pressure)
    fprintf(fid, 'SCALARS Pressure float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%.6E\n', options.pressure);
end

if ~isempty(options.theta_n)
    fprintf(fid, 'SCALARS theta_numerico float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%.6E\n', options.theta_n);
end

if ~isempty(options.k_theta)
    fprintf(fid, 'SCALARS K_theta float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%.6E\n', options.k_theta);
end
%==========================================================================
if ~isempty(options.exact_sol)
    fprintf(fid, 'SCALARS exact_solution float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%.6E\n', options.exact_sol);
end
if ~isempty(options.theta_exact)
    fprintf(fid, 'SCALARS theta_exact_solution float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%.6E\n', options.theta_exact);
end

if ~isempty(options.k_theta_exact)
    fprintf(fid, 'SCALARS K_theta_exact_solution float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%.6E\n', options.k_theta_exact);
end
%==========================================================================
if ~isempty(options.watersaturation)
    if ~isempty(options.satLabel)
        satName = options.satLabel;
    elseif 200 < auxnumcase && auxnumcase < 300
        satName = 'Concentration';
    else
        satName = 'WaterContent';
    end
    fprintf(fid, 'SCALARS %s float 1\n', satName);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%.6E\n', options.watersaturation);
end

if ~isempty(options.normk)
    fprintf(fid, 'SCALARS FlowresultZ float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%.6E\n', options.normk);
end

if ~isempty(options.exact_sol)
    fprintf(fid, 'SCALARS ExactSolution float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%.6E\n', options.exact_sol);
end

fclose(fid);
end