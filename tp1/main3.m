%% TP1: Análisis de Elementos Cáscara — QLLL, QLQL vs TLQL
%  Scordelis-Lo roof benchmark
clear, clc, close all

%% 1. Parámetros
geometria.R       = 25;     % Radio [in]
geometria.L_half  = 25;     % Mitad del largo (simetría) [in]
geometria.phi_max = 40;     % Ángulo de apertura [grados]
geometria.t       = 0.25;   % Espesor [in]

material.E  = 4.32e8;   % Módulo elástico [psi]
material.nu = 0;        % Poisson

q_psi   = 90;           % Carga gravitatoria [psi]
ndofnod = 6;            % DOFs por nodo

meshSizes = [2, 4, 8, 16];
refSize = 40; % Malla de referencia 

%% 2. Bucle de mallas
% a. QLLL
wb_QLLL = zeros(size(meshSizes));
mesh_L  = cell(1, length(meshSizes));
for k = 1:length(meshSizes)
    [mesh_L{k}, wb_QLLL(k)] = resolverMalla(meshSizes(k), geometria, material, q_psi, ndofnod, 'QLLL');
end

% b. QLQL
wb_QLQL = zeros(size(meshSizes));
mesh_Q  = cell(1, length(meshSizes));
for k = 1:length(meshSizes)
    [mesh_Q{k}, wb_QLQL(k)] = resolverMalla(meshSizes(k), geometria, material, q_psi, ndofnod, 'QLQL');
end

%% 3. Convergencia de w_B
[meshRef,sol_ref] = resolverMalla(refSize, geometria, material, q_psi, ndofnod, 'QLLL'); % Malla de referencia
 
figure('Color','w', 'Name','Convergencia w_B');
plot(meshSizes, wb_QLLL, 'k--+', 'LineWidth',1.4, 'MarkerSize',9, 'DisplayName','QLLL');
hold on;
plot(meshSizes, wb_QLQL, 'b--d', 'LineWidth',1.4, 'MarkerSize',9, 'DisplayName','QLQL');

yline(sol_ref, 'r--', 'LineWidth',1.2, ...
      'Label',sprintf('Ref. Oñate: %.4f ft', sol_ref), ...
      'LabelHorizontalAlignment','left', ...
      'DisplayName','Referencia');
grid on;
set(gca,'YDir','reverse');
ylim([-0.55 -0.20]);  xlim([0 18]);
xlabel('Densidad de malla (N\timesN)','FontSize',11);
ylabel('w_B [in]','FontSize',11);
title('Convergencia de w_B','FontSize',12);
legend('Location','East','FontSize',10);

%% 6. Recuperación de esfuerzos — malla 8×8
mesh8 = find(meshSizes == 8);

[theta_L, Nx_L, My_L, Qy_L] = extraerLinea(mesh_L{mesh8}, geometria, material, meshSizes(mesh8));
[theta_Q, Nx_Q, My_Q, Qy_Q] = extraerLinea(mesh_Q{mesh8}, geometria, material, meshSizes(mesh8));
[theta_Ref, Nx_Ref, My_Ref, Qy_Ref] = extraerLinea(meshRef, geometria, material, refSize);

%% 7. Gráficos de distribución A-B
% Nx'
figure('Color','w','Name','Nx a lo largo de A-B');
plot(theta_L, Nx_L/1e4, 'k--+','LineWidth',1.2,'MarkerSize',7,'DisplayName','QLLL (8x8)'); hold on;
plot(theta_Q, Nx_Q/1e4, 'b--d','LineWidth',1.2,'MarkerSize',7,'DisplayName','QLQL (8x8)');
plot(theta_Ref, Nx_Ref/1e4, 'r--','LineWidth',1.2,'MarkerSize',7,'DisplayName','QLLL (40x40)');
grid on;
xlabel('\theta [°]','FontSize',11);  ylabel('N_{x''} \times10^{-4} [lb/in]','FontSize',11);
title('Distribución de N_{x''} — línea A-B','FontSize',12);
legend('Location','NorthWest','FontSize',10);

% My'
figure('Color','w','Name','My a lo largo de A-B');
plot(theta_L, My_L, 'k--+','LineWidth',1.2,'MarkerSize',7,'DisplayName','QLLL (8x8)'); hold on;
plot(theta_Q, My_Q, 'b--d','LineWidth',1.2,'MarkerSize',7,'DisplayName','QLQL (8x8)');
plot(theta_Ref, My_Ref, 'r--','LineWidth',1.2,'MarkerSize',7,'DisplayName','QLLL (40x40)');
grid on;
xlabel('\theta [°]','FontSize',11);  ylabel('M_{y''} [lb·in/in]','FontSize',11);
title('Distribución de M_{y''} — línea A-B','FontSize',12);
legend('Location','SouthEast','FontSize',10);

% Qy'
figure('Color','w','Name','Qy a lo largo de A-B');
plot(theta_L, Qy_L, 'k--+','LineWidth',1.2,'MarkerSize',7,'DisplayName','QLLL (8x8)'); hold on;
plot(theta_Q, Qy_Q, 'b--d','LineWidth',1.2,'MarkerSize',7,'DisplayName','QLQL (8x8)');
plot(theta_Ref, Qy_Ref, 'r--','LineWidth',1.2,'MarkerSize',7,'DisplayName','QLLL (40x40)');
yline(0,'k:','LineWidth',0.8,'HandleVisibility','off');
grid on; ylim([0 350]);
xlabel('\theta [°]','FontSize',11);  ylabel('Q_{y''} [lb/in]','FontSize',11);
title('Distribución de Q_{y''} — línea A-B','FontSize',12);
legend('Location','North','FontSize',10);
%% 8. Tabla de diferencias punto a punto (malla 8x8)

fprintf('\n  Diferencia de wb (porcentaje respecto a referencia):\n');
perc_diff_QLLL = (wb_QLLL - sol_ref)./sol_ref * 100;
perc_diff_QLQL = (wb_QLQL - sol_ref)./sol_ref * 100;
for k = 1:length(meshSizes)
    fprintf('    Malla %2dx%-2d: QLLL - Ref = %+.4f%%, QLQL - Ref = %+.4f%%\n', ...
        meshSizes(k), meshSizes(k), perc_diff_QLLL(k), perc_diff_QLQL(k));
end
%% ================================================================ %%
%%  FUNCIONES LOCALES
%% ================================================================ %%

function [msh, wb] = resolverMalla(N, geometria, material, q_psi, ndofnod, elemType)
    % 1. Determinar tipo de elemento geométrico para la malla
    if strcmpi(elemType, 'TLQL')
        geomType = 'T3';
    else
        geomType = 'Q4';
    end

    % 2. Generar malla
    msh = GenerarMallaTecho(N, geometria, geomType, false);
    msh.numNodes = size(msh.nodes, 1);
    msh.numElem  = size(msh.connect, 1);

    % 3. Condiciones de Borde (Compatibles con DOFs globales)
    tol = 1e-3;
    msh.bc = true(msh.numNodes, ndofnod);
    msh.bc(msh.nodes(:,2) < tol,                      [1 3])   = false; % Pared rígida (Y=0)
    msh.bc(msh.nodes(:,2) > geometria.L_half - tol,   [2 4 6]) = false; % Simetría Y (Y=L/2)
    msh.bc(msh.nodes(:,1) < tol,                      [1 5 6]) = false; % Simetría cumbrera (X=0)

    % 4. Vector de Cargas F (Adaptable a 3 o 4 nodos por área tributaria)
    msh.F = zeros(msh.numNodes * ndofnod, 1);
    for e = 1:msh.numElem
        nd = msh.connect(e,:);
        p  = msh.nodes(nd,:);
        
        % Área de un triángulo (1-2-3)
        A = 0.5 * norm(cross(p(2,:)-p(1,:), p(3,:)-p(1,:)));
        % Si es Q4, sumamos el área del segundo triángulo (1-3-4)
        if length(nd) == 4
            A = A + 0.5 * norm(cross(p(3,:)-p(1,:), p(4,:)-p(1,:)));
        end
        
        fne = q_psi * A / length(nd);
        for i = 1:length(nd)
            msh.F((nd(i)-1)*ndofnod + 3) = msh.F((nd(i)-1)*ndofnod + 3) - fne;
        end
    end

    % 5. Ensamblaje K global
    nDOF_e = size(msh.connect, 2) * ndofnod; % 18 (T3) o 24 (Q4)
    I = zeros(msh.numElem * nDOF_e^2, 1);
    J = zeros(msh.numElem * nDOF_e^2, 1);
    V = zeros(msh.numElem * nDOF_e^2, 1);
    cnt = 0;

    for e = 1:msh.numElem
        nd   = msh.connect(e,:);
        crds = msh.nodes(nd,:);

        % Asignación de elemento (Totalmente modular)
        switch upper(elemType)
            case 'QLLL'
                Ke = CalcularRigidezQLLL(crds, geometria, material, 'selective');
            case 'QLQL'
                Ke = CalcularRigidezQLQL(crds, geometria, material);
        end

        % Mapeo de DOFs para N nodos
        dofs = zeros(1, nDOF_e);
        for i = 1:length(nd)
            s = (i-1)*ndofnod + 1;
            dofs(s:s+5) = (nd(i)-1)*ndofnod + (1:6);
        end

        [gJ, gI] = meshgrid(dofs, dofs);
        r = cnt+1 : cnt+nDOF_e^2;
        I(r) = gI(:);  
        J(r) = gJ(:);  
        V(r) = Ke(:);
        cnt  = cnt + nDOF_e^2;
    end

    % 6. Solución del sistema esparso
    K = sparse(I, J, V, msh.numNodes*ndofnod, msh.numNodes*ndofnod);
    fl = find(msh.bc');
    U = zeros(msh.numNodes*ndofnod, 1);
    U(fl) = K(fl,fl) \ msh.F(fl);
    msh.U = U;

    % 7. Extracción de w_B (punto B, nodo central inferior)
    [~, nc] = max(msh.nodes(:,1) + msh.nodes(:,2)); 
    wb = U((nc-1)*ndofnod + 3);
end

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
            theta_out = [theta_out; rad2deg(atan2(ctr(1), ctr(3)))]; 
            Nx_out    = [Nx_out;    N_res(1,1,e)]; 
            My_out    = [My_out;    M_res(2,1,e)]; 
            Qy_out    = [Qy_out;    Q_res(2,1,e)]; 
        end
    end

    [theta_out, s] = sort(theta_out);
    Nx_out = Nx_out(s);  My_out = My_out(s);  Qy_out = Qy_out(s);
end