function plotPlate3D(nodes, connectivity, w_nodes, colorField, scaleFactor, titleStr)
% plotPlate3D  3D surface plot of a deformed Mindlin plate
%
%   plotPlate3D(nodes, connectivity, w_nodes, colorField, scaleFactor, titleStr)
%
%   Inputs:
%     nodes        - (nNod x 2)  undeformed XY node coordinates
%     connectivity - (nEle x 4)  Q4 element connectivity (node IDs, 1-based)
%     w_nodes      - (nNod x 1)  transverse deflection at each node [m]
%     colorField   - (nNod x 1) OR (nEle x 1)  scalar to color the surface
%                    e.g. w_nodes (deflection) or element centroid moments
%                    Pass [] to color by deflection automatically.
%     scaleFactor  - displacement amplification factor (use 1 for true scale,
%                    larger values help visualise thin-plate deformation)
%                    Pass [] to auto-scale so max deflection = a/10.
%     titleStr     - title string for the figure ('' for no title)
%
%   The plot shows:
%     - Deformed surface, colored by colorField
%     - Element mesh edges drawn on the deformed surface
%     - Undeformed flat plate as a semi-transparent grey reference
%     - Colorbar and axis labels

nEle = size(connectivity, 1);
nNod = size(nodes, 1);

% --- Defaults ---
if isempty(colorField)
    colorField = w_nodes;
end

if isempty(scaleFactor)
    aLen        = max(nodes(:,1)) - min(nodes(:,1));
    wMax        = max(abs(w_nodes));
    scaleFactor = (aLen / 10) / max(wMax, eps);
end

% --- Deformed z-coordinates ---
z_def = w_nodes * scaleFactor;     % amplified transverse displacement

% --- Determine color mode (per-node or per-element) ---
perNode = (length(colorField) == nNod);

% --- Build patch vertex/face arrays ---
% Vertices: [x, y, z_deformed]
Verts = [nodes, z_def];

% Build CData
if perNode
    % Interpolate across each face from nodal values
    CData = colorField;           % (nNod x 1)  -> Gouraud shading
else
    CData = colorField;           % (nEle x 1)  -> flat face color
end

% --- Draw deformed surface ---
if perNode
    h = patch('Faces', connectivity, ...
              'Vertices', Verts, ...
              'FaceVertexCData', CData, ...
              'FaceColor', 'interp', ...
              'EdgeColor', [0.25 0.25 0.25], ...
              'LineWidth', 0.6, ...
              'FaceAlpha', 1.0);
else
    h = patch('Faces', connectivity, ...
              'Vertices', Verts, ...
              'FaceVertexCData', CData, ...
              'FaceColor', 'flat', ...
              'EdgeColor', [0.25 0.25 0.25], ...
              'LineWidth', 0.6, ...
              'FaceAlpha', 1.0);
end

hold on;

% --- Draw undeformed plate as transparent wireframe reference ---
Verts0 = [nodes, zeros(nNod, 1)];
patch('Faces', connectivity, ...
      'Vertices', Verts0, ...
      'FaceColor', [0.8 0.8 0.8], ...
      'EdgeColor', [0.5 0.5 0.5], ...
      'LineWidth', 0.5, ...
      'FaceAlpha', 0.15, ...
      'LineStyle', '--');

% --- Colorbar and labels ---
colormap('jet');
cb = colorbar;
if perNode && isequal(colorField, w_nodes)
    cb.Label.String = 'w  [m]';
else
    cb.Label.String = 'Field value';
end

xlabel('x  [m]', 'FontSize', 11);
ylabel('y  [m]', 'FontSize', 11);
zlabel(sprintf('w \\times %g  [m]', scaleFactor), 'FontSize', 11);

if ~isempty(titleStr)
    title(titleStr, 'FontSize', 12);
end

% --- 3D view settings ---
view([-35, 25]);
axis tight;
grid on;
box on;
lighting gouraud;
camlight('headlight');
set(gca, 'DataAspectRatio', [1, 1, max(range(z_def), eps) / (0.3 * range(nodes(:,1)))]);

% --- Annotation: scale factor ---
xLim = xlim; yLim = ylim; zLim = zlim;
text(xLim(1), yLim(2), zLim(2), ...
    sprintf('  scale \\times%g', scaleFactor), ...
    'FontSize', 9, 'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'bottom');

hold off;

end
