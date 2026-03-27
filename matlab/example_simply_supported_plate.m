% example_simply_supported_plate.m
%
% Simply-supported square plate under uniform transverse load q0.
%
% Boundary conditions (soft simply-supported, SSSS):
%   w = 0 on all four edges (phi_x, phi_y left free).
%
% Reference solution (Navier series, Timoshenko & Woinowsky-Krieger):
%   w_center = alpha * q0 * a^4 / D
%   alpha    = 0.00406  (infinite series, converged value)
%   D        = E*t^3 / (12*(1 - nu^2))   [plate bending stiffness]
%
% The script solves the problem for increasing mesh refinements and plots
% the convergence of the center deflection toward the analytical solution.

clear; clc; close all;

% =========================================================================
% 1. PROBLEM PARAMETERS
% =========================================================================
a   = 1.0;     % plate side length [m]
t   = 0.01;    % plate thickness   [m]  (t/a = 0.01 -> thin plate)
q0  = 1e3;     % uniform load      [N/m^2]

material.E  = 210e9;
material.nu = 0.3;
material.G  = material.E / (2*(1 + material.nu));

integrationMethod = 'selective';   % 'full' | 'selective' | 'reduced'
nDofNod = 3;                       % DOFs per node: [w, phi_x, phi_y]
npgLoad = 2;                       % Gauss points per direction for load

% =========================================================================
% 2. ANALYTICAL REFERENCE
% =========================================================================
D        = material.E * t^3 / (12*(1 - material.nu^2));
w_exact  = 0.00406 * q0 * a^4 / D;
fprintf('Plate bending stiffness D = %.4e N.m\n', D);
fprintf('Analytical center deflection (Navier): w = %.6e m\n\n', w_exact);

% =========================================================================
% 3. CONSTITUTIVE MATRICES  (same for all mesh refinements)
% =========================================================================
[Cb, Cs] = getElementConstitutive(material, t);

% =========================================================================
% 4. MESH REFINEMENTS AND CONVERGENCE LOOP
% =========================================================================
meshSizes = [2, 4, 8, 16];   % number of elements per side
w_fem     = zeros(size(meshSizes));

for iMesh = 1:length(meshSizes)

    n = meshSizes(iMesh);    % elements per side

    % --- 4a. Generate structured Q4 mesh --------------------------------
    nNodX = n + 1;
    nNodY = n + 1;
    nNod  = nNodX * nNodY;
    nEle  = n * n;

    [X, Y] = meshgrid(linspace(0, a, nNodX), linspace(0, a, nNodY));
    nodes  = [X(:), Y(:)];    % (nNod x 2)

    % Q4 connectivity: counter-clockwise, bottom-left origin
    connectivity = zeros(nEle, 4);
    iEle = 0;
    for jj = 1:n
        for ii = 1:n
            iEle = iEle + 1;
            n1 = (jj-1)*nNodX + ii;
            n2 = (jj-1)*nNodX + ii + 1;
            n3 =  jj   *nNodX + ii + 1;
            n4 =  jj   *nNodX + ii;
            connectivity(iEle,:) = [n1, n2, n3, n4];
        end
    end

    % --- 4b. DOF assignment ---------------------------------------------
    nDofTot  = nDofNod * nNod;
    nodeDofs = reshape(1:nDofTot, nDofNod, nNod)';   % (nNod x 3)

    nNodEle  = 4;
    nDofEle  = nDofNod * nNodEle;
    elementDofs = zeros(nEle, nDofEle);
    for ie = 1:nEle
        eleNodes          = connectivity(ie,:);
        elementDofs(ie,:) = reshape(nodeDofs(eleNodes,:)', 1, []);
    end

    % Build elements struct (required by getCharacteristicStresses)
    elements.connectivity = connectivity;
    elements.property     = ones(nEle, 1);      % all elements same property
    elements.dof          = elementDofs;

    % Single property entry
    properties(1).Cb = Cb;
    properties(1).Cs = Cs;

    % --- 4c. Stiffness matrix assembly ----------------------------------
    elementStiffness = cell(nEle, 1);
    rowIndex         = cell(nEle, 1);
    columnIndex      = cell(nEle, 1);

    for ie = 1:nEle
        eleDofs  = elements.dof(ie,:)';
        nodesEle = nodes(connectivity(ie,:), :);
        [elementStiffness{ie}, ~, ~, rowIndex{ie}, columnIndex{ie}] = ...
            getStiffnessMatrix_Plates(nodesEle, eleDofs, Cb, Cs, integrationMethod);
    end

    K = assembleSystem(elementStiffness, rowIndex, columnIndex);

    % --- 4d. Load vector assembly ---------------------------------------
    R = zeros(nDofTot, 1);

    for ie = 1:nEle
        eleDofs  = elements.dof(ie,:)';
        nodesEle = nodes(connectivity(ie,:), :);
        fz = q0 * ones(nNodEle, 1);
        mx = zeros(nNodEle, 1);
        my = zeros(nNodEle, 1);
        [Rele, ~] = applyLoadOnElement(nodesEle, eleDofs, fz, mx, my, npgLoad);
        R(eleDofs) = R(eleDofs) + Rele;
    end

    % --- 4e. Boundary conditions: simply supported on all four edges ----
    %   Soft SSSS: w = 0 on x=0, x=a, y=0, y=a
    [fd1, ~] = fixBoundary(nodes, 1, 0, 1, nDofNod);
    [fd2, ~] = fixBoundary(nodes, 1, a, 1, nDofNod);
    [fd3, ~] = fixBoundary(nodes, 2, 0, 1, nDofNod);
    [fd4, ~] = fixBoundary(nodes, 2, a, 1, nDofNod);

    fixedDofs = unique([fd1; fd2; fd3; fd4]);
    freeDofs  = setdiff((1:nDofTot)', fixedDofs);

    % --- 4f. Solve ------------------------------------------------------
    D_vec = solveSystem(K, R, freeDofs);

    % --- 4g. Extract center deflection ----------------------------------
    %   Find node closest to (a/2, a/2)
    dist      = sqrt((nodes(:,1) - a/2).^2 + (nodes(:,2) - a/2).^2);
    [~, iCtr] = min(dist);
    w_fem(iMesh) = D_vec(nDofNod*(iCtr-1) + 1);

    fprintf('Mesh %dx%d  |  w_center = %.6e m  |  error = %.2f%%\n', ...
        n, n, w_fem(iMesh), 100*abs(w_fem(iMesh) - w_exact)/w_exact);
end

% =========================================================================
% 5. POST-PROCESSING  (finest mesh)
% =========================================================================
fprintf('\n--- Post-processing on %dx%d mesh ---\n', n, n);

% Moment resultants at element centroids
evaluationPoints = [0, 0];
[M_res, Q_res] = getCharacteristicStresses(nodes, elements, properties, D_vec, evaluationPoints);

Mx_ele = squeeze(M_res(1, 1, :));   % Mx at centroid of each element
My_ele = squeeze(M_res(2, 1, :));   % My
Qx_ele = squeeze(Q_res(1, 1, :));   % Qx

% Center element index (element closest to plate center)
eleCenter = zeros(nEle, 2);
for ie = 1:nEle
    eleCenter(ie,:) = mean(nodes(connectivity(ie,:), :));
end
dist_ele    = sqrt((eleCenter(:,1) - a/2).^2 + (eleCenter(:,2) - a/2).^2);
[~, iCtrEle] = min(dist_ele);
fprintf('Mx at center element = %.4e N.m/m\n', Mx_ele(iCtrEle));
fprintf('My at center element = %.4e N.m/m\n', My_ele(iCtrEle));

% Analytical Mx at center (Navier):  Mx_center = beta * q0 * a^2
%   beta ≈ 0.0479 for nu = 0.3 (Timoshenko)
Mx_exact = 0.0479 * q0 * a^2;
fprintf('Analytical Mx at center ≈ %.4e N.m/m\n', Mx_exact);

% =========================================================================
% 6. PLOTS
% =========================================================================

w_nodes = D_vec(1:nDofNod:end);   % w at every node  (nNod x 1)

% --- Convergence plot ---
figure('Name','Convergence');
semilogx(meshSizes, w_fem*1e6, 'o-b', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
yline(w_exact*1e6, '--r', 'LineWidth', 1.5, 'Label', 'Analytical (Navier)');
xlabel('Elements per side');
ylabel('Center deflection  w  [\mum]');
title('Convergence — Simply-supported plate under uniform load');
legend('FEM (selective integration)', 'Navier solution', 'Location','southeast');
grid on;

% --- 3D deformed plate — colored by deflection (auto scale factor) ---
figure('Name','Deformed plate — deflection');
plotPlate3D(nodes, connectivity, w_nodes, w_nodes, [], ...
    sprintf('Deformed plate  w  [m] — %d\\times%d mesh (amplified)', n, n));

% --- 3D deformed plate — colored by Mx (element centroid values) ---
figure('Name','Deformed plate — Mx');
plotPlate3D(nodes, connectivity, w_nodes, Mx_ele, [], ...
    sprintf('Deformed plate colored by M_x  [N\\cdotm/m] — %d\\times%d mesh', n, n));
colorbar;   % update label
cb = findobj(gcf, 'Type', 'ColorBar');
cb.Label.String = 'M_x  [N\cdotm/m]';

% --- 2D top-view Mx contour (element centroid values) ---
figure('Name','Moment Mx — top view');
for ie = 1:nEle
    xv = nodes(connectivity(ie,:), 1);
    yv = nodes(connectivity(ie,:), 2);
    patch(xv, yv, Mx_ele(ie)*ones(4,1), Mx_ele(ie), 'EdgeColor', 'none');
end
colormap('jet'); colorbar;
xlabel('x [m]'); ylabel('y [m]');
title('Bending moment  M_x  [N\cdotm/m] at element centroids');
axis equal tight;

% --- Mesh visualization ---
figure('Name','Mesh');
meshplot(connectivity, nodes, 'k', false);
title(sprintf('Q4 mesh — %d \\times %d elements', n, n));
xlabel('x [m]'); ylabel('y [m]');
