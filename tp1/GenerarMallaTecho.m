function mesh = GenerarMallaTecho(N,geometria,elemType,Check)
    
    if nargin < 4
        Check = false;
    end

    phi_max = geometria.phi_max;
    L_half = geometria.L_half;
    R = geometria.R;

    phi_val = linspace(0, deg2rad(phi_max), N+1); 
    y_val = linspace(0, L_half, N+1);             
    [PHI, Y] = meshgrid(phi_val, y_val);
    
    % Coordenadas Cartesianas 
    mesh.nodes = [R * sin(PHI(:)), Y(:), R * cos(PHI(:))];
    
    mesh.connect = [];
    numNodes = size(mesh.nodes, 1);
    nodesGrid = reshape(1:numNodes, N+1, N+1);
    
    switch elemType
        case 'Q4'
            % Generación de cuadriláteros de 4 nodos
            for j = 1:N % Recorre en dirección longitudinal (Y)
                for i = 1:N % Recorre en dirección del arco (Phi)
                    n1 = nodesGrid(j, i);
                    n2 = nodesGrid(j+1, i);
                    n3 = nodesGrid(j+1, i+1);
                    n4 = nodesGrid(j, i+1);
                    mesh.connect = [mesh.connect; n1, n2, n3, n4];
                end
            end
        case 'T3'
            % Generación de triángulos de 3 nodos (partiendo cada quad en dos)
            for j = 1:N
                for i = 1:N
                    n1 = nodesGrid(j, i);
                    n2 = nodesGrid(j+1, i);
                    n3 = nodesGrid(j+1, i+1);
                    n4 = nodesGrid(j, i+1);
                    % Triángulos 1 y 2
                    mesh.connect = [mesh.connect; n1, n2, n3; n1, n3, n4];
                end
            end
    end
    
    %Graficación de la malla
    if isequal(Check, true)
        figure('Color', 'w', 'Name', 'Malla del Cuarto de Techo');
        patch('Faces', mesh.connect, 'Vertices', mesh.nodes, ...
              'FaceColor', [0.8 0.9 1], 'EdgeColor', 'k', 'LineWidth', 1);
        axis equal; view(3); grid on; hold on;
        xlabel('X (Horizontal) [in]'); ylabel('Y (Longitudinal) [in]'); zlabel('Z (Vertical) [in]');
        
        % Numeración de Nodos (en azul)
        text(mesh.nodes(:,1), mesh.nodes(:,2), mesh.nodes(:,3), string(1:numNodes), ...
            'Color', 'b', 'FontSize', 8, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        
        % Numeración de Elementos (en rojo, al centroide)
        numElem = size(mesh.connect, 1);
        for e = 1:numElem
            nodes_e = mesh.connect(e, :);
            centroid = mean(mesh.nodes(nodes_e, :), 1);
            text(centroid(1), centroid(2), centroid(3)+1, string(e), ...
                'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
        end
        title(['Malla para el techo generada - Tipo: ', elemType]);
    end
end