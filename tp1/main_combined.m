%% main_combined.m
%  Análisis Unificado: QLLL, QLQL y Degenerados (AHMAD4/8/9)
%  Benchmark: Techo de Scordelis-Lo — referencia QLLL 40×40
%
%  Referencias:
%   - Belytschko, T., Stolarski, H., Liu, W.K., Carpenter, N. and Ong, J.S.J.
%     Stress projection for membrane and shear locking in shell finite elements.
%     Comput. Methods Appl. Mech. Engrg., 51(1), 221–258, 1985.
%   - Oñate, E. Structural Analysis with the Finite Element Method.
%     Linear Statics Vol.2: Beams, Plates and Shells.
clear; clc; close all;

%% 0. Rutas de acceso
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir, 'Cascara Plana Degenerada (Top level shit)'));

%% 1. Parámetros comunes
geometria.R       = 25;      % Radio [in]
geometria.L_half  = 25;      % Mitad del largo [in]
geometria.phi_max = 40;      % Ángulo de apertura [°]
geometria.t       = 0.25;    % Espesor [in]

material.E  = 4.32e8;        % Módulo elástico [psi]
material.nu = 0;             % Poisson

% Matriz constitutiva 5×5 para elementos degenerados (Ahmad)
D1         = material.E / (1 - material.nu^2);
G          = material.E / (2*(1 + material.nu));
c          = 5/6;                   % Factor de corrección de Reissner-Mindlin
material.G = G;
material.D = [D1,             material.nu*D1, 0,   0,    0  ;
              material.nu*D1, D1,             0,   0,    0  ;
              0,              0,              G,   0,    0  ;
              0,              0,              0,   c*G,  0  ;
              0,              0,              0,   0,    c*G];

q_val     = 90;               % Carga gravitatoria [psi]
meshSizes = [2, 4, 8, 16];
refSize   = 40;               % Malla de referencia QLLL 40×40
ndofnod6  = 6;                % DOFs/nodo — QLLL/QLQL

% Puntos angulares para evaluación de esfuerzos degenerados
n_pts_eval = 16;
phi_eval   = linspace(0, geometria.phi_max, n_pts_eval);

%% 2. Referencia: QLLL 40×40
fprintf('Calculando referencia QLLL 40x40...\n');
[meshRef, sol_ref_signed] = resolverMalla(refSize, geometria, material, q_val, ndofnod6, 'QLLL');
sol_ref = abs(sol_ref_signed);   % Positivo (desplazamiento hacia abajo)

[theta_Ref, Nx_Ref, My_Ref, Qy_Ref] = extraerLinea(meshRef, geometria, material, refSize);

%% 3. Elementos no degenerados: QLLL y QLQL
fprintf('Calculando QLLL y QLQL...\n');
mesh_L  = cell(1, length(meshSizes));
mesh_Q  = cell(1, length(meshSizes));
wb_QLLL = zeros(size(meshSizes));
wb_QLQL = zeros(size(meshSizes));

for k = 1:length(meshSizes)
    [mesh_L{k}, wb_QLLL(k)] = resolverMalla(meshSizes(k), geometria, material, q_val, ndofnod6, 'QLLL');
    [mesh_Q{k}, wb_QLQL(k)] = resolverMalla(meshSizes(k), geometria, material, q_val, ndofnod6, 'QLQL');
end
wb_QLLL = abs(wb_QLLL);
wb_QLQL = abs(wb_QLQL);

% Esfuerzos en la malla 8×8
ks8 = find(meshSizes == 8);
[theta_L8, Nx_L8, My_L8, Qy_L8] = extraerLinea(mesh_L{ks8}, geometria, material, meshSizes(ks8));
[theta_Q8, Nx_Q8, My_Q8, Qy_Q8] = extraerLinea(mesh_Q{ks8}, geometria, material, meshSizes(ks8));

%% 4. Elementos degenerados: AHMAD4, AHMAD8, AHMAD9
fprintf('Calculando elementos degenerados...\n');
ele_types = {'AHMAD4', 'AHMAD8', 'AHMAD9'};
nTypes    = length(ele_types);
nSizes    = length(meshSizes);

wb_Ah  = zeros(nTypes, nSizes);
Res_Nx = cell(nTypes, nSizes);
Res_My = cell(nTypes, nSizes);
Res_Qy = cell(nTypes, nSizes);

for k_type = 1:nTypes
    eleType = ele_types{k_type};
    for k_size = 1:nSizes
        N = meshSizes(k_size);
        [wB, mesh_ah, U_ah, nodalSys_ah] = resolverMallaAhmad(N, geometria, material, q_val, eleType);
        wb_Ah(k_type, k_size) = wB;
        [Res_Nx{k_type,k_size}, Res_My{k_type,k_size}, Res_Qy{k_type,k_size}] = ...
            CalcularEsfuerzosAHMAD(phi_eval, mesh_ah, U_ah, nodalSys_ah, material, geometria, eleType);
    end
end

%% 5. Figura 1 — Convergencia unificada de w_B
figure('Color','w','Name','Convergencia unificada w_B');
hold on;

% No degenerados
plot(meshSizes, wb_QLLL, 'k--+',  'LineWidth',1.4,'MarkerSize',9,'DisplayName','QLLL');
plot(meshSizes, wb_QLQL, 'b--d',  'LineWidth',1.4,'MarkerSize',9,'DisplayName','QLQL');

% Degenerados
estilos_ah = {'-r^', '-gs', '-mo'};
nombres_ah = {'AHMAD4', 'AHMAD8', 'AHMAD9'};
for k = 1:nTypes
    plot(meshSizes, wb_Ah(k,:), estilos_ah{k}, ...
        'LineWidth',1.4,'MarkerSize',8,'DisplayName',nombres_ah{k});
end

% Referencia QLLL 40×40
yline(sol_ref, '--', 'Color',[0.7 0 0], 'LineWidth',1.6, ...
    'Label', sprintf('QLLL 40\\times40: %.4f in', sol_ref), ...
    'LabelHorizontalAlignment','left', ...
    'DisplayName','Ref. QLLL 40\times40');

grid on; box on;
set(gca,'XTick',meshSizes);
xlim([meshSizes(1)-1, meshSizes(end)+1]);
xlabel('Densidad de malla (N\timesN)','FontSize',11);
ylabel('|w_B| [in]','FontSize',11);
title('Convergencia de desplazamiento vertical en B — Scordelis-Lo','FontSize',12);
legend('Location','best','FontSize',10);
hold off;

%% 6. Figuras 2-4 — Distribución de esfuerzos a lo largo de A-B (malla 8×8)
stress_nd  = {{Nx_L8, My_L8, Qy_L8}, {Nx_Q8, My_Q8, Qy_Q8}};
stress_ref = {Nx_Ref, My_Ref, Qy_Ref};
theta_nd   = {theta_L8, theta_Q8};
Res_all    = {Res_Nx, Res_My, Res_Qy};

titulos_s  = {'Esfuerzo de Membrana N_{x''}', ...
              'Momento Flector M_{y''}', ...
              'Esfuerzo de Corte Q_{y''}'};
ylabels_s  = {'N_{x''} [lb/in]', ...
              'M_{y''} [lb\cdotin/in]', ...
              'Q_{y''} [lb/in]'};
figNames_s = {'Nx — linea A-B', 'My — linea A-B', 'Qy — linea A-B'};
estilos_nd = {'k--+', 'b--d'};
nombres_nd = {'QLLL', 'QLQL'};

for e = 1:3
    figure('Color','w','Name',figNames_s{e});
    hold on;

    % QLLL y QLQL — 8×8
    for k = 1:2
        plot(theta_nd{k}, stress_nd{k}{e}, estilos_nd{k}, ...
            'LineWidth',1.2,'MarkerSize',7, ...
            'DisplayName',[nombres_nd{k} ' (8\times8)']);
    end

    % Degenerados — 8×8
    for k = 1:nTypes
        plot(phi_eval, Res_all{e}{k, ks8}, estilos_ah{k}, ...
            'LineWidth',1.2,'MarkerSize',6, ...
            'DisplayName',[nombres_ah{k} ' (8\times8)']);
    end

    % Referencia QLLL 40×40
    plot(theta_Ref, stress_ref{e}, 'k-', 'LineWidth',2, ...
        'DisplayName','QLLL 40\times40 (Ref.)');

    yline(0,'k:','LineWidth',0.8,'HandleVisibility','off');
    grid on; box on;
    xlabel('\phi [°]','FontSize',11);
    ylabel(ylabels_s{e},'FontSize',11);
    title([titulos_s{e} ' — linea A-B  (8\times8 vs. 40\times40)'],'FontSize',12);
    legend('Location','best','FontSize',9);
    hold off;
end

%% 6b. (OPCIONAL) Figuras separadas — solo QLLL y QLQL vs. Referencia
%  Para activar: eliminar las líneas "%{" y "%}" que rodean el bloque.
%{
% ---- Convergencia: solo no degenerados ----
figure('Color','w','Name','Convergencia QLLL y QLQL');
hold on;
plot(meshSizes, wb_QLLL, 'k--+', 'LineWidth',1.4,'MarkerSize',9,'DisplayName','QLLL');
plot(meshSizes, wb_QLQL, 'b--d', 'LineWidth',1.4,'MarkerSize',9,'DisplayName','QLQL');
yline(sol_ref, '--', 'Color',[0.7 0 0], 'LineWidth',1.6, ...
    'Label',sprintf('QLLL 40\\times40: %.4f in', sol_ref), ...
    'LabelHorizontalAlignment','left', ...
    'DisplayName','Ref. QLLL 40\times40');
grid on; box on;
set(gca,'XTick',meshSizes);
xlim([meshSizes(1)-1, meshSizes(end)+1]);
xlabel('Densidad de malla (N\timesN)','FontSize',11);
ylabel('|w_B| [in]','FontSize',11);
title('Convergencia w_B — QLLL y QLQL','FontSize',12);
legend('Location','best','FontSize',10);
hold off;

% ---- Esfuerzos 8x8: solo no degenerados ----
for e = 1:3
    figure('Color','w','Name',['(No deg.) ' figNames_s{e}]);
    hold on;
    for k = 1:2
        plot(theta_nd{k}, stress_nd{k}{e}, estilos_nd{k}, ...
            'LineWidth',1.2,'MarkerSize',7, ...
            'DisplayName',[nombres_nd{k} ' (8\times8)']);
    end
    plot(theta_Ref, stress_ref{e}, 'k-', 'LineWidth',2, ...
        'DisplayName','QLLL 40\times40 (Ref.)');
    yline(0,'k:','LineWidth',0.8,'HandleVisibility','off');
    grid on; box on;
    xlabel('\phi [°]','FontSize',11);
    ylabel(ylabels_s{e},'FontSize',11);
    title([titulos_s{e} ' — QLLL y QLQL  (8\times8 vs. 40\times40)'],'FontSize',12);
    legend('Location','best','FontSize',9);
    hold off;
end
%}

%% 6c. (OPCIONAL) Figuras separadas — solo Degenerados vs. Referencia
%  Para activar: eliminar las líneas "%{" y "%}" que rodean el bloque.
%{
% ---- Convergencia: solo degenerados ----
figure('Color','w','Name','Convergencia AHMAD4, AHMAD8, AHMAD9');
hold on;
for k = 1:nTypes
    plot(meshSizes, wb_Ah(k,:), estilos_ah{k}, ...
        'LineWidth',1.4,'MarkerSize',8,'DisplayName',nombres_ah{k});
end
yline(sol_ref, '--', 'Color',[0.7 0 0], 'LineWidth',1.6, ...
    'Label',sprintf('QLLL 40\\times40: %.4f in', sol_ref), ...
    'LabelHorizontalAlignment','left', ...
    'DisplayName','Ref. QLLL 40\times40');
grid on; box on;
set(gca,'XTick',meshSizes);
xlim([meshSizes(1)-1, meshSizes(end)+1]);
xlabel('Densidad de malla (N\timesN)','FontSize',11);
ylabel('|w_B| [in]','FontSize',11);
title('Convergencia w_B — Elementos Degenerados (Ahmad)','FontSize',12);
legend('Location','best','FontSize',10);
hold off;

% ---- Esfuerzos 8x8: solo degenerados ----
for e = 1:3
    figure('Color','w','Name',['(Deg.) ' figNames_s{e}]);
    hold on;
    for k = 1:nTypes
        plot(phi_eval, Res_all{e}{k, ks8}, estilos_ah{k}, ...
            'LineWidth',1.2,'MarkerSize',6, ...
            'DisplayName',[nombres_ah{k} ' (8\times8)']);
    end
    plot(theta_Ref, stress_ref{e}, 'k-', 'LineWidth',2, ...
        'DisplayName','QLLL 40\times40 (Ref.)');
    yline(0,'k:','LineWidth',0.8,'HandleVisibility','off');
    grid on; box on;
    xlabel('\phi [°]','FontSize',11);
    ylabel(ylabels_s{e},'FontSize',11);
    title([titulos_s{e} ' — Degenerados  (8\times8 vs. 40\times40)'],'FontSize',12);
    legend('Location','best','FontSize',9);
    hold off;
end
%}

%% 7. Reporte en consola
w_ref_lit = 0.3024;    % Valor bibliográfico de referencia [in]
linea_d   = repmat('=', 1, 100);
linea_s   = repmat('-', 1, 100);

fprintf('\n%s\n', linea_d);
fprintf('  REPORTE UNIFICADO: TECHO DE SCORDELIS-LO — QLLL, QLQL, AHMAD4/8/9\n');
fprintf('%s\n\n', linea_d);
fprintf('  R = %.0f in  |  L = %.0f in  |  t = %.2f in  |  E = %.2e psi  |  nu = %.1f  |  q = %.0f psi\n\n', ...
    geometria.R, geometria.L_half*2, geometria.t, material.E, material.nu, q_val);
fprintf('  Referencia numérica  (QLLL 40x40) : %.4f in\n', sol_ref);
fprintf('  Referencia bibliográfica (Oñate)   : %.4f in\n\n', w_ref_lit);

fprintf('%s\n', linea_s);
fprintf('  CONVERGENCIA  |w_B| [in]  (error %% respecto QLLL 40x40)\n');
fprintf('%s\n', linea_s);

header = sprintf('  %-18s', 'Elemento');
for n = meshSizes
    header = [header, sprintf('|   N=%-2d  (error%%)  ', n)]; %#ok<AGROW>
end
fprintf('%s\n%s\n', header, linea_s);

all_names = [nombres_nd, nombres_ah];
all_wb    = [wb_QLLL; wb_QLQL; wb_Ah];   % 5 × nSizes

for i = 1:5
    fprintf('  %-18s', all_names{i});
    for k = 1:nSizes
        err = (all_wb(i,k) - sol_ref) / sol_ref * 100;
        fprintf('|  %.4f  (%+5.1f%%)  ', all_wb(i,k), err);
    end
    fprintf('\n');
end
fprintf('%s\n\n', linea_d);

%% ================================================================ %%
%%  FUNCIONES LOCALES
%% ================================================================ %%

% ---- resolverMalla (QLLL / QLQL) --------------------------------
function [msh, wb] = resolverMalla(N, geometria, material, q_psi, ndofnod, elemType)
    msh = GenerarMallaTecho(N, geometria, 'Q4', false);
    msh.numNodes = size(msh.nodes, 1);
    msh.numElem  = size(msh.connect, 1);

    tol = 1e-3;
    msh.bc = true(msh.numNodes, ndofnod);
    msh.bc(msh.nodes(:,2) < tol,                      [1 3])   = false;
    msh.bc(msh.nodes(:,2) > geometria.L_half - tol,   [2 4 6]) = false;
    msh.bc(msh.nodes(:,1) < tol,                      [1 5 6]) = false;

    msh.F = zeros(msh.numNodes * ndofnod, 1);
    for e = 1:msh.numElem
        nd = msh.connect(e,:);
        p  = msh.nodes(nd,:);
        A  = 0.5 * norm(cross(p(2,:)-p(1,:), p(3,:)-p(1,:)));
        if length(nd) == 4
            A = A + 0.5 * norm(cross(p(3,:)-p(1,:), p(4,:)-p(1,:)));
        end
        fne = q_psi * A / length(nd);
        for i = 1:length(nd)
            msh.F((nd(i)-1)*ndofnod + 3) = msh.F((nd(i)-1)*ndofnod + 3) - fne;
        end
    end

    nDOF_e = size(msh.connect, 2) * ndofnod;
    I = zeros(msh.numElem * nDOF_e^2, 1);
    J = zeros(msh.numElem * nDOF_e^2, 1);
    V = zeros(msh.numElem * nDOF_e^2, 1);
    cnt = 0;
    for e = 1:msh.numElem
        nd   = msh.connect(e,:);
        crds = msh.nodes(nd,:);
        switch upper(elemType)
            case 'QLLL', Ke = CalcularRigidezQLLL(crds, geometria, material, 'selective');
            case 'QLQL', Ke = CalcularRigidezQLQL(crds, geometria, material);
        end
        dofs = zeros(1, nDOF_e);
        for i = 1:length(nd)
            s = (i-1)*ndofnod + 1;
            dofs(s:s+5) = (nd(i)-1)*ndofnod + (1:6);
        end
        [gJ, gI] = meshgrid(dofs, dofs);
        r = cnt+1 : cnt+nDOF_e^2;
        I(r) = gI(:);  J(r) = gJ(:);  V(r) = Ke(:);
        cnt  = cnt + nDOF_e^2;
    end

    K  = sparse(I, J, V, msh.numNodes*ndofnod, msh.numNodes*ndofnod);
    fl = find(msh.bc');
    U  = zeros(msh.numNodes*ndofnod, 1);
    U(fl) = K(fl,fl) \ msh.F(fl);
    msh.U = U;

    [~, nc] = max(msh.nodes(:,1) + msh.nodes(:,2));
    wb = U((nc-1)*ndofnod + 3);
end

% ---- extraerLinea (QLLL / QLQL) ---------------------------------
function [theta_out, Nx_out, My_out, Qy_out] = extraerLinea(msh, geometria, material, N)
    evalPoints = [0; 0];
    [N_res, M_res, Q_res] = CalcularEsfuerzosCompletos(msh, geometria, material, msh.U, evalPoints);

    tol   = 1e-3;
    L_tgt = geometria.L_half;
    dy    = geometria.L_half / N;

    theta_out = [];  Nx_out = [];  My_out = [];  Qy_out = [];
    for e = 1:msh.numElem
        ctr = mean(msh.nodes(msh.connect(e,:),:), 1);
        if abs(ctr(2) - L_tgt) < (dy/2 + tol)
            theta_out = [theta_out; rad2deg(atan2(ctr(1), ctr(3)))]; %#ok<AGROW>
            Nx_out    = [Nx_out;    N_res(1,1,e)];                   %#ok<AGROW>
            My_out    = [My_out;    M_res(2,1,e)];                   %#ok<AGROW>
            Qy_out    = [Qy_out;    Q_res(2,1,e)];                   %#ok<AGROW>
        end
    end
    [theta_out, s] = sort(theta_out);
    Nx_out = Nx_out(s);  My_out = My_out(s);  Qy_out = Qy_out(s);
end

% ---- resolverMallaAhmad (AHMAD4 / AHMAD8 / AHMAD9) --------------
function [wB, mesh, U, nodalSys] = resolverMallaAhmad(N, geometria, material, q_val, eleType)
    ndofnod = 5;
    switch eleType
        case 'AHMAD4', nnodele = 4;
        case 'AHMAD8', nnodele = 8;
        case 'AHMAD9', nnodele = 9;
    end
    integrationType = 'selective';

    mesh    = GenerarMallaTechoAhmad(N, geometria, eleType);
    nnod    = size(mesh.nodes,  1);
    nele    = size(mesh.connect, 1);
    nDofTot = ndofnod * nnod;

    % --- Coordenadas naturales de los nodos de referencia por tipo ---
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

    % --- Sistema nodal v1, v2, v3 ---
    v_contri_sum = zeros(nnod, 3);
    for iele = 1:nele
        eleNodes    = mesh.connect(iele,:);
        nodalCoords = mesh.nodes(eleNodes,:);
        for lnod = 1:nnodele
            inod   = eleNodes(lnod);
            dN     = shapefunsderAHMAD([ksi_loc(lnod), eta_loc(lnod)], eleType);
            jac    = dN(:,:,1) * nodalCoords;
            v3_ele = cross(jac(1,:), jac(2,:));
            v3_ele = v3_ele / norm(v3_ele);
            if v3_ele(3) < 0, v3_ele = -v3_ele; end
            v_contri_sum(inod,:) = v_contri_sum(inod,:) + v3_ele;
        end
    end

    nodalSys = zeros(3, 3, nnod);
    for inod = 1:nnod
        if norm(v_contri_sum(inod,:)) > 0
            v3 = v_contri_sum(inod,:)' / norm(v_contri_sum(inod,:));
        else
            v3 = [0; 0; 1];
        end
        v1 = cross([0; 1; 0], v3);
        if norm(v1) < 1e-6, v1 = cross(v3, [1; 0; 0]); end
        v1 = v1 / norm(v1);
        v2 = cross(v3, v1);
        v2 = v2 / norm(v2);
        nodalSys(:,:,inod) = [v1, v2, v3];
    end

    % --- Órdenes de integración ---
    if nnodele > 4, n_full = 3; n_red = 2; else, n_full = 2; n_red = 1; end
    switch lower(integrationType)
        case 'full',      order_mf = n_full; order_s = n_full;
        case 'reduced',   order_mf = n_red;  order_s = n_red;
        case 'selective', order_mf = n_full; order_s = n_red;
    end

    D_mf = material.D; D_mf(4:5,:) = 0; D_mf(:,4:5) = 0;   % Membrana/flexión
    D_s  = material.D; D_s(1:3,:)  = 0; D_s(:,1:3)  = 0;   % Corte

    [w_mf, g_mf, ng_mf] = gaussAHMAD(order_mf, 2);
    [w_s,  g_s,  ng_s ] = gaussAHMAD(order_s,  2);

    % --- Ensamblaje de K ---
    K = sparse(nDofTot, nDofTot);
    for iele = 1:nele
        Ke          = zeros(nnodele*5, nnodele*5);
        nodes_e     = mesh.connect(iele,:);
        nodalCoords = mesh.nodes(nodes_e,:);
        v1_e = squeeze(nodalSys(:,1,nodes_e));
        v2_e = squeeze(nodalSys(:,2,nodes_e));
        v3_e = squeeze(nodalSys(:,3,nodes_e));
        t_e  = ones(1, nnodele) * geometria.t;

        for ig = 1:ng_mf
            ksi = g_mf(ig,1); eta = g_mf(ig,2); zeta = g_mf(ig,3);
            [Ni,~]  = shapefunsAHMAD([ksi, eta], eleType);
            dN_nat  = shapefunsderAHMAD([ksi, eta], eleType);
            jac     = shelljac(Ni, dN_nat(:,:,1), zeta, nodalCoords, t_e, v3_e);
            detJ    = abs(det(jac));  invJ = inv(jac);
            e1 = jac(1,:)/norm(jac(1,:));
            e3 = cross(e1, jac(2,:)); e3 = e3/norm(e3);
            e2 = cross(e3, e1);
            R_gauss = [e1; e2; e3];
            B_loc   = getBAHMAD_Local(Ni, dN_nat(:,:,1), invJ, R_gauss, zeta, t_e, v1_e, v2_e);
            Ke      = Ke + B_loc' * D_mf * B_loc * detJ * w_mf(ig);
        end

        for ig = 1:ng_s
            ksi = g_s(ig,1); eta = g_s(ig,2); zeta = g_s(ig,3);
            [Ni,~]  = shapefunsAHMAD([ksi, eta], eleType);
            dN_nat  = shapefunsderAHMAD([ksi, eta], eleType);
            jac     = shelljac(Ni, dN_nat(:,:,1), zeta, nodalCoords, t_e, v3_e);
            detJ    = abs(det(jac));  invJ = inv(jac);
            e1 = jac(1,:)/norm(jac(1,:));
            e3 = cross(e1, jac(2,:)); e3 = e3/norm(e3);
            e2 = cross(e3, e1);
            R_gauss = [e1; e2; e3];
            B_loc   = getBAHMAD_Local(Ni, dN_nat(:,:,1), invJ, R_gauss, zeta, t_e, v1_e, v2_e);
            Ke      = Ke + B_loc' * D_s * B_loc * detJ * w_s(ig);
        end

        dofs = reshape(bsxfun(@plus, (nodes_e-1)*5, (1:5)'), 1, []);
        K(dofs, dofs) = K(dofs, dofs) + Ke;
    end

    % --- Vector de cargas (carga gravitatoria) ---
    F = zeros(nDofTot, 1);
    [w_surf, g_surf, ng_surf] = gaussAHMAD(order_mf, 1);
    for iele = 1:nele
        nodes_e     = mesh.connect(iele,:);
        nodalCoords = mesh.nodes(nodes_e,:);
        v3_e = squeeze(nodalSys(:,3,nodes_e));
        t_e  = ones(1, nnodele) * geometria.t;
        Fe   = zeros(nnodele*ndofnod, 1);
        for ig = 1:ng_surf
            ksi = g_surf(ig,1); eta = g_surf(ig,2); zeta = 0;
            [Ni,~]  = shapefunsAHMAD([ksi, eta], eleType);
            dN_nat  = shapefunsderAHMAD([ksi, eta], eleType);
            jac     = shelljac(Ni, dN_nat(:,:,1), zeta, nodalCoords, t_e, v3_e);
            dA      = norm(cross(jac(1,:), jac(2,:)));
            for i = 1:nnodele
                idx_z     = (i-1)*ndofnod + 3;
                Fe(idx_z) = Fe(idx_z) - Ni(i) * q_val * dA * w_surf(ig);
            end
        end
        dofs = reshape(bsxfun(@plus, (nodes_e-1)*ndofnod, (1:ndofnod)'), 1, []);
        F(dofs) = F(dofs) + Fe/2;
    end

    % --- Condiciones de borde ---
    bc  = true(nnod, ndofnod);
    tol = 1e-3;
    bc(mesh.nodes(:,2) < tol,                     [1 3]) = false;  % Diafragma rígido Y=0
    bc(mesh.nodes(:,2) > geometria.L_half - tol,  [2 5]) = false;  % Simetría longitudinal Y=L/2
    bc(mesh.nodes(:,1) < tol,                     [1 4]) = false;  % Simetría cumbrera X=0
    dofs_libres = find(bc');

    % --- Solución ---
    U = zeros(nDofTot, 1);
    U(dofs_libres) = K(dofs_libres, dofs_libres) \ F(dofs_libres);

    % --- Punto B (borde libre en Y = L_half) ---
    phi_max_rad = deg2rad(geometria.phi_max);
    xB = geometria.R * sin(phi_max_rad);
    yB = geometria.L_half;
    zB = geometria.R * cos(phi_max_rad);
    [~, nodeB] = min(sum((mesh.nodes - [xB, yB, zB]).^2, 2));
    wB = abs(U((nodeB-1)*ndofnod + 3));
end
