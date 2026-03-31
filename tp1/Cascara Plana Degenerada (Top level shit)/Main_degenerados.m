%% Main_degenerados.m - Análisis de Cáscara 3D Degenerada 
clear; clc; close all;

%% 1. Parámetros Iniciales (Scordelis-Lo Benchmark)
geometria.R = 25;          
geometria.L_half = 25;     
geometria.t = 0.25;        % Espesor 
geometria.phi_max = 40;    % Ángulo de apertura [grados]

material.E = 4.32e8;       % 
material.nu = 0;           % Poisson

D1 = material.E / (1 - material.nu^2);    %Plane stress
G = material.E / (2 * (1 + material.nu));
material.G = G;
c = 5/6; % Factor de corrección de corte

D = [ D1,    material.nu*D1, 0,  0,   0;
      material.nu*D1, D1,    0,  0,   0;
      0,     0,     G,  0,   0;
      0,     0,     0,  c*G, 0;
      0,     0,     0,  0,   c*G ];
% Matriz Completa D 
material.D = D;

q_val = 90;                % Carga distribuida

ele_types = {'AHMAD4', 'AHMAD8', 'AHMAD9'};
mesh_sizes = [2, 4, 8, 16];
ndofnod = 5;                   % u, v, w, alpha, beta
integrationType = 'selective'; % Opciones: 'full', 'reduced', 'selective'

% Puntos de evaluación a lo largo del arco (phi de 0 a 40)
n_pts_eval = 16; 
phi_eval = linspace(0, geometria.phi_max, n_pts_eval);

resultados_wB = zeros(length(ele_types), length(mesh_sizes));

% Inicializar contenedores (Celdas para guardar los esfuerzos)
Res_Nx = cell(length(ele_types), length(mesh_sizes));
Res_My = cell(length(ele_types), length(mesh_sizes));
Res_Qy = cell(length(ele_types), length(mesh_sizes));

for k_type = 1:length(ele_types)
    eleType = ele_types{k_type};
    
    % Ajustar nnodele según el tipo de elemento
    switch eleType
        case 'AHMAD4', nnodele = 4;
        case 'AHMAD8', nnodele = 8;
        case 'AHMAD9', nnodele = 9;
    end
    
    %% 2. Generación/Carga de Malla
    for k_size = 1:length(mesh_sizes)
        N = mesh_sizes(k_size);

        mesh = GenerarMallaTechoAhmad(N, geometria, eleType);
        
        nnod = size(mesh.nodes, 1);
        nele = size(mesh.connect, 1);
        nDofTot = ndofnod * nnod;

        %% 3. Definición del Sistema Nodal (v1, v2, v3) adaptado
        nodalSys = zeros(3, 3, nnod);
        v_contri_sum = zeros(nnod, 3);
        v_contri_count = zeros(nnod, 1);
        
        % Definimos las coordenadas naturales de los nodos locales para cada tipo
        switch eleType
            case 'AHMAD4'
                ksi_loc = [-1  1  1 -1];
                eta_loc = [-1 -1  1  1];
            case 'AHMAD8'
                ksi_loc = [-1  1  1 -1  0  1  0 -1];
                eta_loc = [-1 -1  1  1 -1  0  1  0];
            case 'AHMAD9'
                ksi_loc = [-1  1  1 -1  0  1  0 -1  0];
                eta_loc = [-1 -1  1  1 -1  0  1  0  0];
        end
        
        for iele = 1:nele
            eleNodes = mesh.connect(iele, :);
            nodalCoords = mesh.nodes(eleNodes, :);
            
            for lnod = 1:nnodele
                inod = eleNodes(lnod);
                % Calculamos derivadas en el nodo local
                dN = shapefunsderAHMAD([ksi_loc(lnod), eta_loc(lnod)], eleType);
                
                % Jacobiano de la superficie media (2x3)
                jac = dN(:,:,1) * nodalCoords; 
                
                % Vector normal a la superficie media en ese nodo
                v3_ele = cross(jac(1,:), jac(2,:));

                v3_ele = v3_ele / norm(v3_ele);
                
                if v3_ele(3) < 0
                    v3_ele = -v3_ele;
                end
                
                % Acumulamos para promediar después
                v_contri_sum(inod, :) = v_contri_sum(inod, :) + v3_ele;
                
            end
        end
        
        % Normalización final y construcción de la terna
        for inod = 1:nnod
            if norm(v_contri_sum(inod, :)) > 0
                v3 = v_contri_sum(inod, :)' / norm(v_contri_sum(inod, :));
            else
                
                v3 = [0; 0; 1]; 
            end
            
            % Construcción de v1 y v2 (Terna ortonormal)
            v1 = cross([0 1 0]', v3);
            if norm(v1) < 1e-6
                v1 = cross(v3, [1 0 0]');
            end
            v1 = v1 / norm(v1);
            v2 = cross(v3, v1);

            v2 = v2 / norm(v2);
            
            nodalSys(:,:,inod) = [v1, v2, v3];
        end
        
        %% 4. Ensamblaje de la Matriz de Rigidez Global
        K = sparse(nDofTot, nDofTot);
        
        if nnodele > 4 
            n_full = 3; n_red = 2; % Para AHMAD8 y AHMAD9
        else
            n_full = 2; n_red = 1; % Para AHMAD4
        end
        
        % 4.2. Configurar el switch de integración
        switch lower(integrationType)
            case 'full'
                order_mf = n_full; order_s = n_full;
            case 'reduced'
                order_mf = n_red;  order_s = n_red;
            case 'selective'
                order_mf = n_full; order_s = n_red;
            otherwise
                error('Tipo de integración no válido.');
        end
        
        % Separamos la matriz constitutiva
        D_mf = material.D; 
        D_mf(4:5, :) = 0; 
        D_mf(:, 4:5) = 0; % Solo Membrana/Flexión
        
        D_s  = material.D; 
        D_s(1:3, :) = 0; 
        D_s(:, 1:3) = 0;  % Solo Corte
        
        % Reglas de Gauss
        [w_mf, g_mf, ng_mf] = gaussAHMAD(order_mf, 2); 
        [w_s,  g_s,  ng_s]  = gaussAHMAD(order_s, 2);
        
        for iele = 1:nele
            Ke = zeros(nnodele * 5, nnodele * 5);
            nodes_e = mesh.connect(iele, :);
            nodalCoords = mesh.nodes(nodes_e, :);
            v1_e = squeeze(nodalSys(:,1,nodes_e));
            v2_e = squeeze(nodalSys(:,2,nodes_e));
            v3_e = squeeze(nodalSys(:,3,nodes_e));
            t_e  = ones(1, nnodele) * geometria.t;
        
            % --- PARTE 1: Membrana y Flexión ---
            for ig = 1:ng_mf
                ksi = g_mf(ig,1); eta = g_mf(ig,2); zeta = g_mf(ig,3);
                [Ni, ~] = shapefunsAHMAD([ksi, eta], eleType);
                dN_nat = shapefunsderAHMAD([ksi, eta], eleType);
                
                jac = shelljac(Ni, dN_nat(:,:,1), zeta, nodalCoords, t_e, v3_e);
                detJ = abs(det(jac));
                invJ = inv(jac);
                
                % 1. Definir terna local (e1, e2, e3) en el punto de Gauss
                e1 = jac(1,:) / norm(jac(1,:)); 
                e3 = cross(e1, jac(2,:)); e3 = e3 / norm(e3);
                e2 = cross(e3, e1);
                R_gauss = [e1; e2; e3]; % Matriz de rotación 3x3
                
                % 2. Construir B Local (usando derivadas rotadas)
                B_loc = getBAHMAD_Local(Ni, dN_nat(:,:,1), invJ, R_gauss, zeta, t_e, v1_e, v2_e);
                
                Ke = Ke + B_loc' * D_mf * B_loc * detJ * w_mf(ig);
            end
        
            % --- PARTE 2: Corte Transversal ---
            for ig = 1:ng_s
                ksi = g_s(ig,1); eta = g_s(ig,2); zeta = g_s(ig,3);
                [Ni, ~] = shapefunsAHMAD([ksi, eta], eleType);
                dN_nat = shapefunsderAHMAD([ksi, eta], eleType);
                
                jac = shelljac(Ni, dN_nat(:,:,1), zeta, nodalCoords, t_e, v3_e);
                detJ = abs(det(jac));
                invJ = inv(jac);
                
                e1 = jac(1,:) / norm(jac(1,:)); 
                e3 = cross(e1, jac(2,:)); e3 = e3 / norm(e3);
                e2 = cross(e3, e1);
                R_gauss = [e1; e2; e3]; 
                
                B_loc = getBAHMAD_Local(Ni, dN_nat(:,:,1), invJ, R_gauss, zeta, t_e, v1_e, v2_e);
                
                Ke = Ke + B_loc' * D_s * B_loc * detJ * w_s(ig);
            end
        
            % Ensamblaje
            dofs = reshape(bsxfun(@plus, (nodes_e-1)*5, (1:5)'), 1, []);
            K(dofs, dofs) = K(dofs, dofs) + Ke;
        end
        
        %% 5. Vector de Cargas Globales (Carga Gravitatoria)
        F = zeros(nDofTot, 1);
        
        [w_surf, g_surf, ng_surf] = gaussAHMAD(order_mf, 1); 
        
        for iele = 1:nele
            nodes_e = mesh.connect(iele, :);
            nodalCoords = mesh.nodes(nodes_e, :);
            v3_e = squeeze(nodalSys(:,3,nodes_e));
            t_e  = ones(1, nnodele) * geometria.t;
            
            Fe = zeros(nnodele * ndofnod, 1);

            for ig = 1:ng_surf
                ksi = g_surf(ig,1); 
                eta = g_surf(ig,2); 
                zeta = 0; % La carga actúa sobre la superficie media
                
                % Evaluamos funciones de forma y sus derivadas naturales
                [Ni, ~] = shapefunsAHMAD([ksi, eta], eleType);
                dN_nat = shapefunsderAHMAD([ksi, eta], eleType);
                
                % Jacobiano en la superficie media para obtener el diferencial de área
                jac = shelljac(Ni, dN_nat(:,:,1), zeta, nodalCoords, t_e, v3_e);
                
                % dA es la norma del producto vectorial de los vectores tangentes (ksi, eta)
                dA = norm(cross(jac(1,:), jac(2,:)));
                
                % Distribución de la carga a los nodos (Consistente)
                % q_val actúa en la dirección Z global (DOF 3 de cada nodo)
                for i = 1:nnodele
                    idx_z = (i - 1) * ndofnod + 3; % El 3er DOF es la componente Z
                    % CORRECCIÓN: Eliminamos el factor /2
                    Fe(idx_z) = Fe(idx_z) - Ni(i) * q_val * dA * w_surf(ig);
                end
            end

            % Ensamblaje Global Vectorizado (similar a la matriz K)
            dofs = reshape(bsxfun(@plus, (nodes_e-1)*ndofnod, (1:ndofnod)'), 1, []);
            F(dofs) = F(dofs) + Fe/2;
        end
        
        %% 6. Condiciones de Borde
        bc = true(nnod, ndofnod); 
        tol = 1e-3;
        
        % 1. BORDE APOYADO (Diafragma rígido en Y = 0)
        % Restringimos desplazamiento en X (1) y en Z (3). El cilindro no se hunde ni se ensancha ahí.
        idx_apoyo = mesh.nodes(:,2) < tol; 
        bc(idx_apoyo, [1 3]) = false; 
        
        % 2. SIMETRÍA LONGITUDINAL (Corte en el centro del galpón, Y = 25)
        % Restringimos desplazamiento en Y (2) y la rotación beta (5)
        idx_sim_long = mesh.nodes(:,2) > geometria.L_half - tol;
        bc(idx_sim_long, [2 5]) = false;
        
        % 3. SIMETRÍA TRANSVERSAL / CUMBRERA (Corte en el pico del techo, X = 0)
        % Restringimos desplazamiento en X (1) y la rotación alpha (4)
        idx_sim_trans = mesh.nodes(:,1) < tol;
        bc(idx_sim_trans, [1 4]) = false;
        
        dofs_libres = find(bc');
        
        %% 7. Resolución del Sistema
        % Inicializamos el vector de desplazamientos globales (nnod * 5, 1)
        U = zeros(nDofTot, 1);
        
        U(dofs_libres) = K(dofs_libres, dofs_libres) \ F(dofs_libres);
        
        
        %% 8. Post-procesamiento 
        
        % Coordenadas del Punto B (Borde libre en el plano de simetría Y=25)
        phi_max = deg2rad(geometria.phi_max);
        xB = geometria.R * sin(phi_max);
        yB = geometria.L_half;
        zB = geometria.R * cos(phi_max);
        
        % Buscar el nodo más cercano a esas coordenadas
        distancias = sum((mesh.nodes - [xB, yB, zB]).^2, 2);
        [~, nodeB] = min(distancias);
        
        resultados_wB(k_type, k_size) = abs(U((nodeB-1)*ndofnod + 3));
        
        factor_escala = 10; 
        %PlotDeformedAhmad(mesh, U, factor_escala); %Descomentar para ver
        %la geometria deformada
        
        [Nx_vec, My_vec, Qy_vec] = CalcularEsfuerzosAHMAD(phi_eval, mesh, U, nodalSys, material, geometria, eleType);
        
        Res_Nx{k_type, k_size} = Nx_vec;
        Res_My{k_type, k_size} = My_vec;
        Res_Qy{k_type, k_size} = Qy_vec;
        
    end  %Termina el bucle de las distintas densidades de Malla
end %Termina el bucle de los distintos tipos de elementos

%% 9. Generación de Graficos
% Creamos una figura grande para acomodar los 4 gráficos cómodamente
figure('Color', 'w', 'Name', 'Resultados Completos: Scordelis-Lo');


% SUBPLOT 1: Convergencia de mallas
subplot(2, 2, 1);
hold on;
estilos = {'-b+', '-rx', '-gs'}; 
nombres = {'AHMAD4 (Lineal)', 'AHMAD8 (Serendipity)', 'AHMAD9 (Lagrangiano)'};

for k = 1:length(ele_types)
    plot(mesh_sizes, resultados_wB(k, :), estilos{k}, ...
        'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', nombres{k});
end

linea_ref = yline(0.3024, '--k', 'Valor de Referencia (0.3024)', ...
    'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'FontSize', 10);
linea_ref.DisplayName = 'Referencia Oñate';

grid on; box on;
set(gca, 'XScale', 'linear', 'XTick', mesh_sizes);
xlabel('Densidad de Malla (N)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Desplazamiento vertical |W_B| (in)', 'FontSize', 11, 'FontWeight', 'bold');
title('Convergencia: Techo de Scordelis-Lo', 'FontSize', 12);
legend('Location', 'southeast', 'FontSize', 9);
ylim([0, 0.35]); 
xlim([mesh_sizes(1)-1, mesh_sizes(end)+1]);
hold off;

% SUBPLOTS 2, 3 y 4: Esfuerzos N_x, M_y, Q_y
ks = 3; % Graficamos los esfuerzos para la malla N=8 (índice 3)
n_size = mesh_sizes(ks); 

esfuerzos_data = {Res_Nx, Res_My, Res_Qy};
titulos = {'Esfuerzo de Membrana N_x', 'Momento Flector M_y', 'Esfuerzo de Corte Q_y'};
unidades = {'lb/in', 'lb*in/in', 'lb/in'};

% Mapeo de subplots: 2 (Arriba Derecha), 3 (Abajo Izquierda), 4 (Abajo Derecha)
posiciones_subplot = [2, 3, 4]; 

for e = 1:3 
    subplot(2, 2, posiciones_subplot(e));
    hold on;
    
    % Graficamos las curvas para los 3 elementos
    plot(phi_eval, esfuerzos_data{e}{1, ks}, '-b+', 'LineWidth', 1.2, 'DisplayName', 'AHMAD4 (Lineal)');
    plot(phi_eval, esfuerzos_data{e}{2, ks}, '-rx', 'LineWidth', 1.2, 'DisplayName', 'AHMAD8 (Serendipity)');
    plot(phi_eval, esfuerzos_data{e}{3, ks}, '-go', 'LineWidth', 1.2, 'DisplayName', 'AHMAD9 (Lagrangiano)');
    
    grid on; box on;
    xlabel('Ángulo \phi (grados)', 'FontSize', 11);
    ylabel([titulos{e} ' [' unidades{e} ']'], 'FontSize', 11);
    title({titulos{e}, ['Malla de ' num2str(n_size) 'x' num2str(n_size) ' (Y=25)']}, 'FontSize', 11);
    
    % Leyenda un poco más pequeña para no tapar las curvas
    legend('Location', 'best', 'FontSize', 8);
    
    % Línea neutra visual
    yline(0, 'k', 'Alpha', 0.3);
    hold off;
end

%% 10. Reporte Profesional de Resultados y Errores en Consola
w_ref = 0.3024; % ft
ancho_linea = 95;
linea_doble = repmat('=', 1, ancho_linea);
linea_simple = repmat('-', 1, ancho_linea);

fprintf('\n%s\n', linea_doble);
fprintf('%65s\n', '      REPORTE DE ANÁLISIS DE ELEMENTOS FINITOS: TECHO DE SCORDELIS-LO');
fprintf('%s\n\n', linea_doble);

% --- 10.1. Resumen de Parámetros del Modelo ---
fprintf('PARÁMETROS GEOMÉTRICOS Y MATERIALES:\n');
fprintf(' Radio (R): %6.2f in    | Longitud (L): %6.2f in   | Espesor (t): %6.2f in\n', geometria.R, geometria.L_half*2, geometria.t);
fprintf(' Módulo E : %6.2e psi | Poisson (v) : %6.2f      | Carga (q)  : %6.2f psi (z)\n', material.E, material.nu, q_val);
fprintf('%s\n\n', linea_simple);

% --- 10.2. Estudio de Convergencia de Desplazamientos ---
fprintf('ESTUDIO DE CONVERGENCIA INDEPENDIENTE DE LA MALLA - DESPLAZAMIENTO VERTICAL |W_B|\n');
fprintf('Valor Teórico de Referencia: %.4f ft\n', w_ref);
fprintf('%s\n', linea_simple);

% Encabezado dinámico de columnas
fprintf('%-22s', 'Tipo de Elemento');
for n = mesh_sizes
    fprintf('|   N = %-2d (Error%%) ', n);
end
fprintf('\n%s\n', linea_simple);

% Filas de datos de convergencia
for i = 1:length(ele_types)
    fprintf('%-22s', ele_types{i});
    for j = 1:length(mesh_sizes)
        val = resultados_wB(i,j);
        error_rel = abs(val - w_ref)/w_ref * 100;
        fprintf('|  %6.4f (%5.1f%%)  ', val, error_rel);
    end
    fprintf('\n');
end
fprintf('%s\n\n', linea_simple);

% --- 10.3. Metodología de Recuperación de Tensiones ---
fprintf('METODOLOGÍA DE RECUPERACIÓN DE ESFUERZOS (Evaluación sobre Y = %g ft):\n', geometria.L_half);
fprintf(' * Membrana (N_x) y Flexión (M_y): Evaluados en el borde geométrico (eta = 1.0).\n');
fprintf(' * Corte Transversal (Q_y): Evaluado en los Puntos de Barlow (Superconvergentes)\n');
fprintf('   para mitigar el fenómeno de Corte Parásito (Shear Locking).\n');
fprintf('   -> Elementos Lineales (AHMAD4)    : Evaluado en ksi = 0.000\n');
fprintf('   -> Elementos Cuadráticos (AHMAD8/9): Evaluado en ksi = 0.577 (1/sqrt(3))\n');
fprintf('%s\n', linea_simple);

% --- 10.4. Resumen Cuantitativo de Esfuerzos Extremos ---
% Extraemos los valores máximos absolutos de la malla N=8 para dar datos de diseño
idx_malla_resumen = 3; % Índice correspondiente a N=8
fprintf('RESUMEN DE ESFUERZOS EXTREMOS DE DISEÑO (Para la malla de N = %d):\n', mesh_sizes(idx_malla_resumen));
fprintf('%-22s | %-20s | %-20s | %-20s\n', 'Tipo de Elemento', 'Max |N_x| (lb/ft)', 'Max |M_y| (lb*ft/ft)', 'Max |Q_y| (lb/ft)');
fprintf('%s\n', linea_simple);

for i = 1:length(ele_types)
    % Extraemos los vectores de la celda para la malla N=8
    Nx_vector = esfuerzos_data{1}{i, idx_malla_resumen};
    My_vector = esfuerzos_data{2}{i, idx_malla_resumen};
    Qy_vector = esfuerzos_data{3}{i, idx_malla_resumen};
    
    % Calculamos el máximo absoluto
    max_Nx = max(abs(Nx_vector));
    max_My = max(abs(My_vector));
    max_Qy = max(abs(Qy_vector));
    
    fprintf('%-22s | %20.2f | %20.2f | %20.2f\n', ele_types{i}, max_Nx, max_My, max_Qy);
end
fprintf('%s\n\n', linea_doble);
