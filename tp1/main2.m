%TP1: Análisis de Elementos Cáscara
% Parte 1: Techo de Scordelis-Lo con elementos planos
% Elementos: QLLL (§6.7.1 Oñate) y QLQL (§6.7.4 Oñate)
clear, clc, close all

%% 1. Parámetros Iniciales
geometria.R       = 25;     % Radio [in]
geometria.L_half  = 25;     % Mitad del largo por simetría [in]
geometria.phi_max = 40;     % Ángulo de apertura [grados]
geometria.t       = 0.25;   % Espesor [in]

material.E  = 4.32e8;       % Módulo elástico [psi]
material.nu = 0;            % Coeficiente de Poisson

q_psi     = 90;             % Carga distribuida [psi]
sol_ref   = -0.3024;        % wb referencia Oñate [in]
meshSizes = [2, 4, 8, 16];
elemType  = 'Q4';
ndofnod   = 6;

[Cb, Cs] = getElementConstitutive(material, geometria.t);

%% 2. Pre-generación de mallas (compartidas por ambos elementos)
mesh = cell(1, length(meshSizes));
for k = 1:length(meshSizes)
    N = meshSizes(k);
    mesh{k} = GenerarMallaTecho(N, geometria, elemType, false);
    mesh{k}.numNodes = size(mesh{k}.nodes, 1);
    mesh{k}.numElem  = size(mesh{k}.connect, 1);
end

%% =========================================================================
%% PARTE A — QLLL
%% =========================================================================
fprintf('=== Análisis QLLL ===\n');
wb_QLLL   = zeros(size(meshSizes));
mesh_QLLL = mesh;

for k = 1:length(meshSizes)
    fprintf('  Malla %dx%d ...\n', meshSizes(k), meshSizes(k));
    mesh_QLLL{k} = aplicarBC(mesh_QLLL{k}, geometria, ndofnod);
    mesh_QLLL{k} = aplicarCargas(mesh_QLLL{k}, q_psi, ndofnod);

    entradasPE = (size(mesh_QLLL{k}.connect,2) * ndofnod)^2;
    I = zeros(mesh_QLLL{k}.numElem * entradasPE, 1);
    J = zeros(mesh_QLLL{k}.numElem * entradasPE, 1);
    V = zeros(mesh_QLLL{k}.numElem * entradasPE, 1);
    count = 0;

    for e = 1:mesh_QLLL{k}.numElem
        nodes_e     = mesh_QLLL{k}.connect(e, :);
        cords_nodos = mesh_QLLL{k}.nodes(nodes_e, :);
        Ke = CalcularRigidezQLLL(cords_nodos, geometria, material, 'selective');
        [I, J, V, count] = ensamblarTriplete(I, J, V, count, Ke, nodes_e, ndofnod);
    end

    [wb_QLLL(k), mesh_QLLL{k}] = resolverYExtraer(mesh_QLLL{k}, I, J, V, ndofnod);
end

%% =========================================================================
%% PARTE B — QLQL
%% =========================================================================
fprintf('=== Análisis QLQL ===\n');
wb_QLQL   = zeros(size(meshSizes));
mesh_QLQL = mesh;

for k = 1:length(meshSizes)
    fprintf('  Malla %dx%d ...\n', meshSizes(k), meshSizes(k));
    mesh_QLQL{k} = aplicarBC(mesh_QLQL{k}, geometria, ndofnod);
    mesh_QLQL{k} = aplicarCargas(mesh_QLQL{k}, q_psi, ndofnod);

    entradasPE = (size(mesh_QLQL{k}.connect,2) * ndofnod)^2;
    I = zeros(mesh_QLQL{k}.numElem * entradasPE, 1);
    J = zeros(mesh_QLQL{k}.numElem * entradasPE, 1);
    V = zeros(mesh_QLQL{k}.numElem * entradasPE, 1);
    count = 0;

    for e = 1:mesh_QLQL{k}.numElem
        nodes_e     = mesh_QLQL{k}.connect(e, :);
        cords_nodos = mesh_QLQL{k}.nodes(nodes_e, :);
        % QLQL con integración selectiva de corte (igual que QLLL)
        Ke = CalcularRigidezQLQL(cords_nodos, geometria, material, 'selective');
        [I, J, V, count] = ensamblarTriplete(I, J, V, count, Ke, nodes_e, ndofnod);
    end

    [wb_QLQL(k), mesh_QLQL{k}] = resolverYExtraer(mesh_QLQL{k}, I, J, V, ndofnod);
end

%% =========================================================================
%% 3. Gráfico de Convergencia — QLLL vs QLQL
%% =========================================================================
figure('Color', 'w', 'Name', 'Convergencia QLLL vs QLQL');
plot(meshSizes, wb_QLLL, 'k--+', 'LineWidth', 1.2, 'MarkerSize', 8); hold on;
plot(meshSizes, wb_QLQL, 'b--s', 'LineWidth', 1.2, 'MarkerSize', 8);
yline(sol_ref, 'r--', 'Referencia Oñate', 'LineWidth', 1, 'LabelHorizontalAlignment', 'left');
grid on;
set(gca, 'YDir', 'reverse');
xlabel('Densidad de Malla (N × N)');
ylabel('Desplazamiento Vertical w_b [in]');
title('Convergencia de w_b — Techo de Scordelis-Lo');
ylim([-0.55 -0.20]);
xlim([0 18]);
legend('QLLL', 'QLQL', 'Referencia', 'Location', 'NorthEast');

fprintf('\n--- Convergencia w_b [in] ---\n');
fprintf('%-8s %10s %10s %10s\n', 'Malla', 'QLLL', 'QLQL', 'Ref');
for k = 1:length(meshSizes)
    fprintf('%-8d %10.4f %10.4f %10.4f\n', meshSizes(k), wb_QLLL(k), wb_QLQL(k), sol_ref);
end

%% =========================================================================
%% 4. Post-procesamiento — Esfuerzos malla 8×8 (k=3)
%% =========================================================================
k_post     = 3;   % índice malla 8×8
evalPoints = [0; 0];

% QLLL
[N_qlll, M_qlll, Q_qlll] = CalcularEsfuerzosCompletos( ...
    mesh_QLLL{k_post}, geometria, material, mesh_QLLL{k_post}.U, evalPoints);
[theta_qlll, Nx_qlll, My_qlll, Qy_qlll] = extraerLineaAB( ...
    mesh_QLLL{k_post}, geometria, meshSizes(k_post), N_qlll, M_qlll, Q_qlll);

% QLQL
[N_qlql, M_qlql, Q_qlql] = CalcularEsfuerzosCompletos( ...
    mesh_QLQL{k_post}, geometria, material, mesh_QLQL{k_post}.U, evalPoints);
[theta_qlql, Nx_qlql, My_qlql, Qy_qlql] = extraerLineaAB( ...
    mesh_QLQL{k_post}, geometria, meshSizes(k_post), N_qlql, M_qlql, Q_qlql);

%% =========================================================================
%% 5. Gráficos de Esfuerzos comparativos en línea A-B
%% =========================================================================

figure('Color', 'w', 'Name', 'Nx en A-B');
plot(theta_qlll, Nx_qlll/1e4, 'k--d', 'LineWidth', 1.2, 'MarkerSize', 6); hold on;
plot(theta_qlql, Nx_qlql/1e4, 'b--s', 'LineWidth', 1.2, 'MarkerSize', 6);
grid on;
xlabel('Ángulo \theta [grados]');
ylabel('N_x'' (×10^4) [lb/in]');
title('Distribución de N_x'' a lo largo de A-B (malla 8×8)');
legend('QLLL', 'QLQL', 'Location', 'Best');

figure('Color', 'w', 'Name', 'My en A-B');
plot(theta_qlll, My_qlll, 'k--d', 'LineWidth', 1.2, 'MarkerSize', 6); hold on;
plot(theta_qlql, My_qlql, 'b--s', 'LineWidth', 1.2, 'MarkerSize', 6);
grid on;
xlabel('Ángulo \theta [grados]');
ylabel('M_y'' [lb·in/in]');
title('Distribución de M_y'' a lo largo de A-B (malla 8×8)');
legend('QLLL', 'QLQL', 'Location', 'Best');

figure('Color', 'w', 'Name', 'Qy en A-B');
plot(theta_qlll, Qy_qlll, 'k--d', 'LineWidth', 1.2, 'MarkerSize', 6); hold on;
plot(theta_qlql, Qy_qlql, 'b--s', 'LineWidth', 1.2, 'MarkerSize', 6);
grid on;
xlabel('Ángulo \theta [grados]');
ylabel('Q_y'' [lb/in]');
title('Distribución de Q_y'' a lo largo de A-B (malla 8×8)');
legend('QLLL', 'QLQL', 'Location', 'Best');


%% =========================================================================
%  Funciones locales
%% =========================================================================

function mesh_k = aplicarBC(mesh_k, geometria, ndofnod)
    tol = 1e-3;
    mesh_k.bc = true(mesh_k.numNodes, ndofnod);
    % Pared rígida Y=0: u(1), w(3)
    mesh_k.bc(mesh_k.nodes(:,2) < tol, [1 3]) = false;
    % Simetría longitudinal Y=L_half: v(2), θx(4), θz(6)
    mesh_k.bc(mesh_k.nodes(:,2) > geometria.L_half - tol, [2 4 6]) = false;
    % Simetría transversal X=0: u(1), θy(5), θz(6)
    mesh_k.bc(mesh_k.nodes(:,1) < tol, [1 5 6]) = false;
end

function mesh_k = aplicarCargas(mesh_k, q_psi, ndofnod)
    mesh_k.F = zeros(mesh_k.numNodes * ndofnod, 1);
    for e = 1:mesh_k.numElem
        nodes_e = mesh_k.connect(e, :);
        n_en    = length(nodes_e);
        p       = mesh_k.nodes(nodes_e, :);
        if n_en == 3
            Area_e = 0.5 * norm(cross(p(2,:)-p(1,:), p(3,:)-p(1,:)));
        else
            A1 = 0.5 * norm(cross(p(2,:)-p(1,:), p(3,:)-p(1,:)));
            A2 = 0.5 * norm(cross(p(3,:)-p(1,:), p(4,:)-p(1,:)));
            Area_e = A1 + A2;
        end
        Fe_nodal = (q_psi * Area_e) / n_en;
        for i = 1:n_en
            dof_z = (nodes_e(i)-1)*ndofnod + 3;
            mesh_k.F(dof_z) = mesh_k.F(dof_z) - Fe_nodal;
        end
    end
end

function [I, J, V, count] = ensamblarTriplete(I, J, V, count, Ke, nodes_e, ndofnod)
    numN   = length(nodes_e);
    dofs_e = zeros(1, numN * ndofnod);
    for i = 1:numN
        s = (i-1)*ndofnod + 1;
        dofs_e(s:s+5) = (nodes_e(i)-1)*ndofnod + (1:6);
    end
    [gJ, gI] = meshgrid(dofs_e, dofs_e);
    n2 = (numN*ndofnod)^2;
    I(count+1:count+n2) = gI(:);
    J(count+1:count+n2) = gJ(:);
    V(count+1:count+n2) = Ke(:);
    count = count + n2;
end

function [wb, mesh_k] = resolverYExtraer(mesh_k, I, J, V, ndofnod)
    nDOF        = mesh_k.numNodes * ndofnod;
    K_global    = sparse(I, J, V, nDOF, nDOF);
    bc_lineal   = mesh_k.bc';
    dofs_libres = find(bc_lineal(:));
    U           = zeros(nDOF, 1);
    U(dofs_libres) = K_global(dofs_libres, dofs_libres) \ mesh_k.F(dofs_libres);
    mesh_k.U    = U;
    [~, nctrl]  = max(mesh_k.nodes(:,1) + mesh_k.nodes(:,2));
    wb          = U((nctrl-1)*ndofnod + 3);
end

function [theta_plot, Nx_plot, My_plot, Qy_plot] = extraerLineaAB( ...
        mesh_k, geometria, N, N_res, M_res, Q_res)
    tol      = 1e-3;
    L_target = geometria.L_half;
    dy_elem  = geometria.L_half / N;
    theta_plot = []; Nx_plot = []; My_plot = []; Qy_plot = [];
    for e = 1:mesh_k.numElem
        c3D = mean(mesh_k.nodes(mesh_k.connect(e,:), :), 1);
        if abs(c3D(2) - L_target) < (dy_elem/2 + tol)
            theta_plot = [theta_plot; rad2deg(atan2(c3D(1), c3D(3)))];
            Nx_plot    = [Nx_plot;    N_res(1,1,e)];
            My_plot    = [My_plot;    M_res(2,1,e)];
            Qy_plot    = [Qy_plot;    Q_res(2,1,e)];
        end
    end
    [theta_plot, idx] = sort(theta_plot);
    Nx_plot = Nx_plot(idx); My_plot = My_plot(idx); Qy_plot = Qy_plot(idx);
end