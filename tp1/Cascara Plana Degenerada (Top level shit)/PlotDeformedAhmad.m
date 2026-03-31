function PlotDeformedAhmad(mesh, U, scale)
    % mesh: estructura con .nodes y .connect
    % U: vector de desplazamientos (5 DOFs por nodo)
    % scale: factor de escala (ej. 50)
    
    nodes = mesh.nodes;
    connect = mesh.connect;
    
    % Extraer desplazamientos traslacionales (u, v, w)
    u = U(1:5:end);
    v = U(2:5:end);
    w = U(3:5:end);
    
    % Calcular nuevas posiciones
    nodes_def = nodes + scale * [u, v, w];
    
    figure('Color', 'w', 'Name', 'Estructura Deformada');
    
    % Graficamos la estructura original (transparente)
    patch('Faces', connect(:,1:4), 'Vertices', nodes, ...
          'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8], 'LineStyle', ':');
    hold on;
    
    % Graficamos la estructura deformada
    % Usamos el desplazamiento vertical w para el mapa de colores
    patch('Faces', connect(:,1:4), 'Vertices', nodes_def, ...
          'FaceVertexCData', abs(w), 'FaceColor', 'interp', 'EdgeColor', 'k');
    
    colorbar;
    colormap('jet');
    title(['Deformada (Factor de escala: ' num2str(scale) ')']);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(3); axis equal; grid on;
    
    % Añadir un punto rojo en el Punto B para identificarlo
    % (Buscamos de nuevo el nodo B en la lista deformada)
    phi_max = deg2rad(25); % Ajustar según tu geometría
    [~, idxB] = min(sum((nodes - [max(nodes(:,1)), max(nodes(:,2)), min(nodes(:,3))]).^2, 2));
    plot3(nodes_def(idxB,1), nodes_def(idxB,2), nodes_def(idxB,3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(nodes_def(idxB,1), nodes_def(idxB,2), nodes_def(idxB,3), '  Punto B', 'FontSize', 12, 'FontWeight', 'bold');
end