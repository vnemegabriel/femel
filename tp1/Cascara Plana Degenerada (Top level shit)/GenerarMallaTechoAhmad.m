function mesh = GenerarMallaTechoAhmad(N, geometria, eleType,check)

    if nargin < 4
        check = false; 
    end % Por defecto no grafica
    
    R = geometria.R;
    L = geometria.L_half;
    phi_max = deg2rad(geometria.phi_max);


    % 1. Generar Nodos
    nDiv = 2 * N; % Doble de puntos para tener nodos en los puntos medios 
    phi_vals = linspace(0, phi_max, nDiv + 1);
    y_vals = linspace(0, L, nDiv + 1);
    [PHI, Y] = meshgrid(phi_vals, y_vals);
    
    nodes_raw = [R * sin(PHI(:)), Y(:), R * cos(PHI(:))];
    node_map = reshape(1:size(nodes_raw,1), nDiv + 1, nDiv + 1);
    connect_raw = [];
    
    % Generación de elementos saltando de a 2 nodos
    for j = 1:2:nDiv
        for i = 1:2:nDiv
            % Nodos de las esquinas (Orden Anti-horario para detJ positivo)
            n1 = node_map(i, j);
            n2 = node_map(i+2, j);
            n3 = node_map(i+2, j+2);
            n4 = node_map(i, j+2);
            
            switch eleType
                case 'AHMAD4' % Usamos solo las esquinas
                    connect_raw = [connect_raw; n1, n2, n3, n4];
                case 'AHMAD8'
                    m1 = node_map(i+1, j);   % medio 1-2
                    m2 = node_map(i+2, j+1); % medio 2-3
                    m3 = node_map(i+1, j+2); % medio 3-4
                    m4 = node_map(i, j+1);   % medio 4-1
                    connect_raw = [connect_raw; n1, n2, n3, n4, m1, m2, m3, m4];
                case 'AHMAD9'
                    m1 = node_map(i+1, j); m2 = node_map(i+2, j+1);
                    m3 = node_map(i+1, j+2); m4 = node_map(i, j+1);
                    c9 = node_map(i+1, j+1); % central
                    connect_raw = [connect_raw; n1, n2, n3, n4, m1, m2, m3, m4, c9];
            end
        end
    end

    % --- PROCESO DE LIMPIEZA DE NODOS (Crucial para evitar matriz singular) ---
    [used_nodes, ~, new_indices] = unique(connect_raw);
    mesh.nodes = nodes_raw(used_nodes, :);
    mesh.connect = reshape(new_indices, size(connect_raw));
    
    if check
        figure('Color', 'w','Name', 'Malla Generada');
        patch('Faces', mesh.connect(:,1:4), 'Vertices', mesh.nodes, ...
              'FaceColor', 'cyan', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
        hold on;
        plot3(mesh.nodes(:,1), mesh.nodes(:,2), mesh.nodes(:,3), 'r.', 'MarkerSize', 10);
        view(3); axis equal; grid on;
        title(['Malla de Elementos: ' eleType ' ( ' num2str(N) 'x' num2str(N) ' ) ']);
    end
end