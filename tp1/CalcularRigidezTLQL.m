function [Ke_global] = CalcularRigidezTLQL(nodes3D, geometry, material)

    E = material.E;  nu = material.nu;  t = geometry.t;
    G = E / (2*(1+nu));

    % 1. Sistema Local Consistente (3 nodos)
    [T_full, nodesLoc] = getTransformationMatrixTLQL(nodes3D);
    x = nodesLoc(:,1); y = nodesLoc(:,2);
    
    % Área geométrica del triángulo
    Area = 0.5 * abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)));

    % 2. Matrices Constitutivas
    Dm = (E*t / (1-nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    Db = (E*t^3 / (12*(1-nu^2))) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    Ds = (5/6) * G * t * eye(2);

    Ke_local = zeros(18, 18);
    
    dN_nat = [-1, 1, 0; 
              -1, 0, 1];
    J = dN_nat * nodesLoc(:,1:2);
    dN_xy = J \ dN_nat; % Derivadas físicas dN/dx y dN/dy

    Bm = zeros(3, 18); Bb = zeros(3, 18);
    for k = 1:3
        col = 6*(k-1);
        
        % --- Membrana (CST) ---
        Bm(1, col+1) = dN_xy(1,k);   % du/dx
        Bm(2, col+2) = dN_xy(2,k);   % dv/dy
        Bm(3, col+1) = dN_xy(2,k);   % du/dy
        Bm(3, col+2) = dN_xy(1,k);   % dv/dx
        
        % --- Flexión (Curvaturas Constantes) ---
        Bb(1, col+5) =  dN_xy(1,k);  % d(theta_y)/dx
        Bb(2, col+4) = -dN_xy(2,k);  % -d(theta_x)/dy
        Bb(3, col+5) =  dN_xy(2,k);  % d(theta_y)/dy
        Bb(3, col+4) = -dN_xy(1,k);  % -d(theta_x)/dx
    end
    
    % Integración de Membrana y Flexión (1 punto es exacto para CST)
    Ke_local = Ke_local + (Bm' * Dm * Bm) * Area;
    Ke_local = Ke_local + (Bb' * Db * Bb) * Area;

    % Corte
    pts = [0.5, 0.0;   % Borde 1-2
           0.5, 0.5;   % Borde 2-3
           0.0, 0.5];  % Borde 3-1
    w_gauss = Area / 3;
    
    for i = 1:3
        xi = pts(i,1); eta = pts(i,2);
        
        % Matriz Bs reconstruida según la cinemática del TLQL para 
        % imponer compatibilidad tangencial en los mid-nodes.
        Bs = getBs_TLQL(xi, eta, J);
        
        % Ensamblaje Gaussiano (3 puntos integran exactamente la cuadrática)
        Ke_local = Ke_local + (Bs' * Ds * Bs) * w_gauss;
    end

    % =====================================================================
    % RIGIDEZ DE DRILLING (6to Grado de Libertad)
    % =====================================================================
    K_drill = zeros(18, 18);
    drill_stiff = 1e-4 * G * t * Area; 
    for k = 1:3
        idx = 6*(k-1) + 6;
        K_drill(idx, idx) = drill_stiff;
    end
    Ke_local = Ke_local + K_drill;

    % 4. Transformación al Sistema Global
    Ke_global = T_full' * Ke_local * T_full;
end

%% ========================================================================
% FUNCIONES LOCALES DE UTILIDAD
%% ========================================================================

function Bs = getBs_TLQL(xi, eta, J)
    % Calcula la matriz de deformación por corte asumida 
    % Evalúa las componentes covariantes exactas en los 3 tie-points (A,B,C)
    
    % A = borde 1-2 (xi=0.5, eta=0)
    BgA = get_covariant(0.5, 0.0, 1, J); 
    
    % B = borde 3-1 (xi=0, eta=0.5)
    BgB = get_covariant(0.0, 0.5, 2, J);
    
    % C = borde 2-3 (xi=0.5, eta=0.5)
    BgC_xi  = get_covariant(0.5, 0.5, 1, J);
    BgC_eta = get_covariant(0.5, 0.5, 2, J);
    
    % Factor de relajación TLQL (Constante c transversal)
    c = BgC_eta - BgC_xi - BgB + BgA;
    
    % Interpolación asumida de las deformaciones covariantes
    Bg_xi  = BgA - eta * c;
    Bg_eta = BgB + xi  * c;
    Bg = [Bg_xi; Bg_eta];
    
    % Transformar a coordenadas físicas (x, y)
    Bs = J \ Bg;
end

function Bg_dir = get_covariant(xi, eta, dir, J)
    % Evalúa el esfuerzo covariante puro de una dirección dada
    N = [1 - xi - eta, xi, eta];
    dN_nat = [-1, 1, 0; 
              -1, 0, 1];
              
    Bg_dir = zeros(1, 18);
    for k = 1:3
        col = 6*(k-1);
        if dir == 1 % Componente a lo largo de xi
            Bg_dir(col+3) = dN_nat(1,k);            % dw/dxi
            Bg_dir(col+4) = -N(k) * J(1,2);         % -theta_x * dy/dxi (CORREGIDO)
            Bg_dir(col+5) =  N(k) * J(1,1);         % +theta_y * dx/dxi (CORREGIDO)
        else        % Componente a lo largo de eta
            Bg_dir(col+3) = dN_nat(2,k);            % dw/deta
            Bg_dir(col+4) = -N(k) * J(2,2);         % -theta_x * dy/deta (CORREGIDO)
            Bg_dir(col+5) =  N(k) * J(2,1);         % +theta_y * dx/deta (CORREGIDO)
        end
    end
end

function [T_full, nodesLoc] = getTransformationMatrixTLQL(nodes3D)
    % Define el sistema coordenado local para un triángulo de 3 nodos
    P1 = nodes3D(1, :); P2 = nodes3D(2, :); P3 = nodes3D(3, :);
    
    v1 = P2 - P1;                    v1 = v1 / norm(v1);
    v3 = cross(v1, P3 - P1);         v3 = v3 / norm(v3);
    v2 = cross(v3, v1);              v2 = v2 / norm(v2);
    
    % Matriz de rotación elemental (3x3) y para nodo completo (6x6)
    T3 = [v1; v2; v3];
    T6 = blkdiag(T3, T3);
    
    % Matriz de transformación global a local (18x18)
    T_full = blkdiag(T6, T6, T6);
    
    % Nodos mapeados al plano local 2D
    nodesLoc = zeros(3, 3);
    for i = 1:3
        nodesLoc(i, :) = (nodes3D(i, :) - P1) * T3';
    end
end