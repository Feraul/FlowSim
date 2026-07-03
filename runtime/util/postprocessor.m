


function postprocessor(pressure, watersaturation,normk ,...
    time, env, step, parmRichardEq)

coord   = env.geometry.coord;
elem    = env.geometry.elem;
filepath = env.mainpathfolders.path;
resfolder = env.mainpathfolders.resfolder;

auxnumcase = env.config.numcase;
% if auxnumcase > 400
%     normk = parmRichardEq.normperm;
% else
%     normk = env.geometry.normperm;
% end

% Arquivo VTU
fname_vtu = fullfile(filepath, resfolder, ['res_00' num2str(step) '.vtu']);
fid = fopen(fname_vtu,'w');

% -------------------------------------------------
% Cabeçalho XML
% -------------------------------------------------
fprintf(fid, '<?xml version="1.0"?>\n');
fprintf(fid, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
fprintf(fid, '  <UnstructuredGrid>\n');
fprintf(fid, '    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n', size(coord,1), size(elem,1));

% -------------------------------------------------
% PONTOS
% -------------------------------------------------
fprintf(fid, '      <Points>\n');
fprintf(fid, '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
fprintf(fid, '          %.16E %.16E %.16E\n', coord');
fprintf(fid, '        </DataArray>\n');
fprintf(fid, '      </Points>\n');

% -------------------------------------------------
% CELLS
% -------------------------------------------------
% Conta triângulos e quadriláteros
triMask  = elem(:,4)==0;
quadMask = ~triMask;
nelem = size(elem,1);

% CELLS connectivity
fprintf(fid, '      <Cells>\n');

% Connectivity
fprintf(fid, '        <DataArray type="Int32" Name="connectivity" format="ascii">\n');
% Triângulos
if any(triMask)
    tri = elem(triMask,1:3)-1; % 0-based
    fprintf(fid, '          %d %d %d\n', tri');
end
% Quadriláteros
if any(quadMask)
    quad = elem(quadMask,1:4)-1;
    fprintf(fid, '          %d %d %d %d\n', quad');
end
fprintf(fid, '        </DataArray>\n');

% Offsets
fprintf(fid, '        <DataArray type="Int32" Name="offsets" format="ascii">\n');
offset = cumsum(3*triMask + 4*quadMask);
fprintf(fid, '          %d\n', offset);
fprintf(fid, '        </DataArray>\n');

% Cell types
fprintf(fid, '        <DataArray type="UInt8" Name="types" format="ascii">\n');
celltypes = zeros(nelem,1,'uint8');
celltypes(triMask) = 5; % VTK_TRIANGLE
celltypes(quadMask) = 9; % VTK_QUAD
fprintf(fid, '          %d\n', celltypes);
fprintf(fid, '        </DataArray>\n');

fprintf(fid, '      </Cells>\n');

% -------------------------------------------------
% CELL DATA
% -------------------------------------------------
fprintf(fid, '      <CellData Scalars="Pressure">\n');

% PRESSURE
fprintf(fid, '        <DataArray type="Float32" Name="Pressure" format="ascii">\n');
fprintf(fid, '          %.16E\n', pressure);
fprintf(fid, '        </DataArray>\n');

% Saturation / Concentration
if 200<auxnumcase && auxnumcase<300
    satName = 'Concentration';
else
    satName = 'ExactSolution';
end
fprintf(fid, '        <DataArray type="Float32" Name="%s" format="ascii">\n', satName);
fprintf(fid, '          %.16E\n', watersaturation);
fprintf(fid, '        </DataArray>\n');

% Norm permeability
fprintf(fid, '        <DataArray type="Float32" Name="FlowresultZ" format="ascii">\n');
fprintf(fid, '          %.16E\n', normk);
fprintf(fid, '        </DataArray>\n');

fprintf(fid, '      </CellData>\n');

% -------------------------------------------------
% Fecha XML
% -------------------------------------------------
fprintf(fid, '    </Piece>\n');
fprintf(fid, '  </UnstructuredGrid>\n');
fprintf(fid, '</VTKFile>\n');

fclose(fid);

end