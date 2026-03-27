%TP1: Análisis de Elementos Cascara

% Parte 1: Techo de Scordelis-Lo con elementos planos
clear,clc, close all

%% 1. Parámetros Iniciales
geometria.R = 25;          % Radio [in]
geometria.L_half = 25;     % Mitad del largo por simetría [in]
geometria.phi_max = 40;    % Ángulo de apertura [grados]
geometria.t = 0.25;        % Espesor [in]

material(1).E = 4.32e8;      % Módulo elástico [psi] 
material(1).nu = 0;          % Coeficiente de Poisson

q_psi = 90;                 %En psi

[Cb, Cs] = getElementConstitutive(material(1), geometria.t);
properties(1).Cb = Cb;
properties(1).Cs = Cs;

%% 2. Generación de malla 
%Se va a generar una malla cuadrada de N x N elementos para que se pueda ir cambiando facilmente
meshSizes = [2, 4, 8, 16]; 
elemType = 'Q4';                      % Elegir entre Q4 o T3
wb_results = zeros(size(meshSizes));
ndofnod = 6;

mesh = cell(1, length(meshSizes));

for k = 1:length(meshSizes)
    N = meshSizes(k);
    mesh{k} = GenerarMallaTecho(N,geometria,elemType,true);
    mesh{k}.numNodes = size(mesh{k}.nodes, 1);
    mesh{k}.numElem = size(mesh{k}.connect, 1);

    mesh{k}.properties(1).Cb = Cb;
    mesh{k}.properties(1).Cs = Cs;

    %% 3. BC
    mesh{k}.bc = true(mesh{k}.numNodes, ndofnod);         %true = dof libre
    tol = 1e-3;
    
    % 1. Pared en Y = 0
    % Restringimos u(1) y w(3) porque el benchmark dice "Rigid Diaphragm"
    idx_pared = mesh{k}.nodes(:,2) < tol;
    mesh{k}.bc(idx_pared, [1 3]) = false;
    
    % 2. Simetría Longitudinal (Y = 25)
    % Restringimos v (2), theta_x (4) y theta_z (6)
    idx_sim_long = mesh{k}.nodes(:,2) > geometria.L_half - tol;
    mesh{k}.bc(idx_sim_long, [2 4 6]) = false;
    
    % 3. Simetría Transversal / Cumbrera (X = 0)
    % Restringimos u (1), theta_y (5) y theta_z (6)
    idx_sim_trans = mesh{k}.nodes(:,1) < tol;
    mesh{k}.bc(idx_sim_trans, [1 5 6]) = false;


    %% 4. Vector de Cargas Globales (Chequear)  
    mesh{k}.F = zeros(mesh{k}.numNodes * 6, 1); 
         
    for e = 1:mesh{k}.numElem
        nodes_e = mesh{k}.connect(e, :);
        n_en = length(nodes_e);           % Número de nodos del elemento (3 o 4)
        p = mesh{k}.nodes(nodes_e, :); 
        
        % Cálculo del área según el tipo de elemento
        if n_en == 3 % Triángulo (T3)
            Area_e = 0.5 * norm(cross(p(2,:)-p(1,:), p(3,:)-p(1,:)));
        elseif n_en == 4 % Cuadrilátero (Q4)
            A1 = 0.5 * norm(cross(p(2,:)-p(1,:), p(3,:)-p(1,:)));
            A2 = 0.5 * norm(cross(p(3,:)-p(1,:), p(4,:)-p(1,:)));
            Area_e = A1 + A2;
        end
        
        % Distribución de carga consistente 
        Fe_nodal = (q_psi * Area_e) / n_en;
        
        for i = 1:n_en
            node_idx = nodes_e(i);
            % Aplicamos la carga en el DOF 3 (Z global)
            global_dof_z = (node_idx - 1) * 6 + 3;
            mesh{k}.F(global_dof_z) = mesh{k}.F(global_dof_z) - Fe_nodal;
        end
    end
    
    %% 5. Ensamblaje de la Matriz de Rigidez Global y Resolución
    % Pre-asignación de memoria para eficiencia (Tripletas I, J, V)
    % Cada elemento Q4 aporta 24x24 = 576 entradas a la matriz global
    entradasPorElemento = (size(mesh{k}.connect,2) * ndofnod)^2;
    I = zeros(mesh{k}.numElem * entradasPorElemento, 1);
    J = zeros(mesh{k}.numElem * entradasPorElemento, 1);
    V = zeros(mesh{k}.numElem * entradasPorElemento, 1);
    
    count = 0;
    
    for e = 1:mesh{k}.numElem
        % a. Obtener conectividad y coordenadas 3D del elemento
        nodes_e = mesh{k}.connect(e, :);
        cords_nodos = mesh{k}.nodes(nodes_e, :);
        
        % b. Calcular matriz de rigidez elemental (Local -> Global)
        % Usamos integración 'selective' para evitar el shear locking
        Ke = CalcularRigidezQLLL(cords_nodos, geometria, material(1), 'selective');
        
        % c. Determinar los Grados de Libertad (DOFs) globales del elemento
        numNodosElem = length(nodes_e);
        dofs_e = zeros(1, numNodosElem * ndofnod);
        for i = 1:numNodosElem
            idx_inicio = (i-1)*ndofnod + 1;
            dofs_e(idx_inicio : idx_inicio + 5) = (nodes_e(i)-1)*ndofnod + (1:6);
        end
        
        % d. Llenado de tripletas para la matriz sparse
        [grid_J, grid_I] = meshgrid(dofs_e, dofs_e);
        pos_inicio = count + 1;
        pos_fin = count + (numNodosElem * ndofnod)^2;
        
        I(pos_inicio:pos_fin) = grid_I(:);
        J(pos_inicio:pos_fin) = grid_J(:);
        V(pos_inicio:pos_fin) = Ke(:);
        
        count = pos_fin;
    end
    
    % e. Creación de la matriz sparse
    K_global = sparse(I, J, V, mesh{k}.numNodes * ndofnod, mesh{k}.numNodes * ndofnod);
    
    % f. Aplicación de condiciones de contorno (Reducción del sistema)
    % Convertimos la matriz 'bc' (numNodes x 6) a un vector lineal de DOFs
    bc_lineal = mesh{k}.bc'; 
    dofs_libres = find(bc_lineal(:)); % Identifica los 'true' (libres)
    
    % g. Resolución del sistema K * U = F
    U = zeros(mesh{k}.numNodes * ndofnod, 1);
    U(dofs_libres) = K_global(dofs_libres, dofs_libres) \ mesh{k}.F(dofs_libres);
    
    mesh{k}.U = U;

    % h. Extracción de resultados (Punto de control)
    % Buscamos el desplazamiento vertical (Z) en el centro del borde libre
    % En nuestro cuarto de modelo: X = R*sin(40), Y = 25
    [~, nodo_control] = max(mesh{k}.nodes(:,1) + mesh{k}.nodes(:,2)); 
    wb_results(k) = U((nodo_control-1)*ndofnod + 3);

end

%% 6. Gráfico de Convergencia 
sol_ref = -0.3024; % Valor de referencia

figure('Color', 'w', 'Name', 'Convergencia');
p = plot(meshSizes, wb_results, 'k--+', 'LineWidth', 1.2, 'MarkerSize', 8); 
hold on;
yline(sol_ref, 'r--', 'Solución Oñate', 'LineWidth', 1, 'LabelHorizontalAlignment', 'left');
grid on;
set(gca, 'YDir', 'reverse'); % 
set(gca, 'XLim', [2 16]);    %
xlabel('Densidad de Malla (N x N)');
ylabel('Desplazamiento Vertical w_b [in]');
title('Convergencia de Wb');
ylim([-0.55 -0.20]); 
xlim([0 18]); 
legend('Elemento QLLL', 'Referencia', 'Location', 'NorthEast');

%% 7. Post-procesamiento: Recuperación para malla 8x8
k = 3; % Malla 8x8
% Evaluamos en el centroide de cada elemento (0,0) para mayor precisión
evalPoints = [0; 0]; 

[N_res, M_res, Q_res] = CalcularEsfuerzosCompletos(mesh{k}, geometria, material(1), mesh{k}.U, evalPoints);

%% 8. Extracción de Esfuerzos sobre la línea A-B (Y = 25)
tol = 1e-3;
L_target = geometria.L_half;

% Inicializar vectores para el gráfico
theta_plot = [];
Nx_plot = [];
My_plot = [];
Qy_plot = [];

dy_elem = geometria.L_half / meshSizes(k);

for e = 1:mesh{k}.numElem
    coords_e = mesh{k}.nodes(mesh{k}.connect(e, :), :);
    centroide_3D = mean(coords_e, 1);
    
    % Filtro robusto: ¿Está el centroide en la última franja de elementos?
    if abs(centroide_3D(2) - L_target) < (dy_elem / 2 + tol)
        th = atan2(centroide_3D(1), centroide_3D(3));
        theta_plot = [theta_plot; rad2deg(th)];
        
        Nx_plot = [Nx_plot; N_res(1, 1, e)]; % Nx es el índice 1
        My_plot = [My_plot; M_res(2, 1, e)]; % My es el índice 2
        Qy_plot = [Qy_plot; Q_res(2, 1, e)]; % Qy es el índice 2
    end
end

% Ordenar para evitar el en el gráfico
[theta_plot, sortIdx] = sort(theta_plot);
Nx_plot = Nx_plot(sortIdx);
My_plot = My_plot(sortIdx);
Qy_plot = Qy_plot(sortIdx);


%% 9. Gráficos de Resultados

% --- Gráfico de Nx ---
figure('Color', 'w', 'Name', 'Distribución de Nx en A-B');
plot(theta_plot, Nx_plot/1e4, 'k--d', 'LineWidth', 1.2); % Divido por 1e4 como en Oñate
grid on; title('Distribución de N_x'' (x10^4)'); xlabel('\theta'); ylabel('N_x'' [lb/in]');

% --- Gráfico de My' ---
figure('Color', 'w', 'Name', 'Distribución de My en A-B');
plot(theta_plot, My_plot, 'k--o', 'LineWidth', 1.2, 'MarkerSize', 6);
grid on;
xlabel('Ángulo \theta [grados]');
ylabel('Momento M_y'' [lb*in/in]');
title('Distribución de M_y'' a lo largo de A-B');


% --- Gráfico de Qy' ---
figure('Color', 'w', 'Name', 'Distribución de Qy en A-B');
plot(theta_plot, Qy_plot, 'k--s', 'LineWidth', 1.2, 'MarkerSize', 6);
grid on;
xlabel('Ángulo \theta [grados]');
ylabel('Esfuerzo de Corte Q_y'' [lb/in]');
title('Distribución de Q_y'' a lo largo de A-B');
