function [Ke_global] = CalcularRigidezQLQL(nodes3D, geometry, material)
    % CalcularRigidezQLQL — Rigidez de cáscara plana QLQL (Oñate §6.7.4)
    %
    % FORMULACIÓN:
    %   w  : bilineal Q4 (igual que QLLL)
    %   θ  : cuadrática incompleta — 4 DOFs de esquina (θx,θy) + 4 DOFs
    %        jerárquicos de lado (Δθs5…Δθs8, uno por arista, Ec. 6.105)
    %   γ  : campo de corte lineal sustituto (mismo que QLLL, Ec. 6.86)
    %   Integración: 2×2 Gauss full para todos los términos
    %   Los 4 DOFs jerárquicos se condensan estáticamente → Ke final 24×24
    %
    % Signatura compatible con CalcularRigidezQLLL (integrationType omitido
    % porque QLQL siempre usa 2×2 full y no necesita integración selectiva).
    %
    % nodes3D  : [4×3] coordenadas globales de los 4 nodos
    % geometry : struct con campo .t
    % material : struct con campos .E y .nu

    E  = material.E;
    nu = material.nu;
    t  = geometry.t;

    % 1. Sistema local y transformación global (24×24) --------------------
    [T, nodesLocal] = getTransformationMatrix(nodes3D);

    % 2. Matrices constitutivas -------------------------------------------
    Dm = (E*t         / (1-nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    Db = (E*t^3/12    / (1-nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    G  = E / (2*(1+nu));
    Ds = (5/6)*G*t * eye(2);

    % 3. Vectores unitarios tangentes de cada arista (2D local) -----------
    % Numeración Fig. 6.19: arista s5=1-2, s6=2-3, s7=3-4, s8=4-1
    sideNodes = [1,2; 2,3; 3,4; 4,1];
    e_side = zeros(4,2);
    for s = 1:4
        v = nodesLocal(sideNodes(s,2),:) - nodesLocal(sideNodes(s,1),:);
        e_side(s,:) = v / norm(v);
    end

    % 4. Integración 2×2 Gauss --------------------------------------------
    [gp, gw] = getGaussPoints(2);
    Area_elem = 0;

    % Sistema expandido: 28 DOFs locales
    %   1-24 → DOFs estándar de nodos [u,v,w,θx,θy,θz] × 4
    %   25-28 → Δθs5,Δθs6,Δθs7,Δθs8 (jerárquicos, sin DOF global)
    Kee = zeros(28, 28);

    for i = 1:2
        for j = 1:2
            xi  = gp(i);  eta = gp(j);
            wj  = gw(i) * gw(j);

            [Bm, Bb_exp, detJ] = getMatricesMB_QLQL(nodesLocal, xi, eta, e_side);
            [Bs,         ~   ] = getMatrixS(nodesLocal, xi, eta);

            wJ = detJ * wj;
            Area_elem = Area_elem + wJ;

            % Membrana: u,v de 4 nodos → DOFs locales [1,2, 7,8, 13,14, 19,20]
            idxM = [1,2, 7,8, 13,14, 19,20];
            Kee(idxM, idxM) = Kee(idxM, idxM) + Bm' * Dm * Bm * wJ;

            % Flexión expandida: θx,θy esquinas + Δθs lados
            %   DOFs locales esquinas : [4,5, 10,11, 16,17, 22,23]
            %   DOFs jerárquicos      : [25, 26, 27, 28]
            idxB = [4,5, 10,11, 16,17, 22,23, 25,26,27,28];
            Kee(idxB, idxB) = Kee(idxB, idxB) + Bb_exp' * Db * Bb_exp * wJ;

            % Corte: w,θx,θy de 4 nodos → DOFs locales [3,4,5, 9,10,11, 15,16,17, 21,22,23]
            % (Los DOFs jerárquicos no contribuyen a γ en QLQL)
            idxS = [3,4,5, 9,10,11, 15,16,17, 21,22,23];
            Kee(idxS, idxS) = Kee(idxS, idxS) + Bs' * Ds * Bs * wJ;
        end
    end

    % 5. Rigidez de drilling (igual que QLLL) -----------------------------
    k_drill = 1e-3 * E * t * Area_elem;
    idxD = [6, 12, 18, 24];
    for d = 1:4
        Kee(idxD(d), idxD(d)) = Kee(idxD(d), idxD(d)) + k_drill;
    end

    % 6. Condensación estática de los 4 DOFs jerárquicos ------------------
    % Partición: [Kaa Kab; Kba Kbb] con fuerzas externas nulas en DOFs b
    % → Ke_condensada = Kaa - Kab * Kbb^{-1} * Kba
    Kaa = Kee(1:24, 1:24);
    Kab = Kee(1:24, 25:28);
    Kbb = Kee(25:28, 25:28);
    Ke_local = Kaa - Kab * (Kbb \ Kab');

    % 7. Transformación al sistema global ---------------------------------
    Ke_global = T' * Ke_local * T;
end


%% =========================================================================
%  Funciones auxiliares (privadas al archivo)
%% =========================================================================

function [gp, gw] = getGaussPoints(n)
    if n == 1
        gp = 0;  gw = 2;
    elseif n == 2
        gp = [-1/sqrt(3), 1/sqrt(3)];
        gw = [1, 1];
    end
end

function [N, dN] = shapeQ4(xi, eta)
    N  = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    dN = 0.25 * [-(1-eta),  (1-eta), (1+eta), -(1+eta);
                 -(1-xi),  -(1+xi),  (1+xi),   (1-xi) ];
end

function [Bm, Bb_exp, detJ] = getMatricesMB_QLQL(nodesLoc, xi, eta, e_side)
    % Retorna:
    %   Bm     : [3×8]  — membrana (idéntico a QLLL)
    %   Bb_exp : [3×12] — flexión expandida
    %            cols 1-8  → contribución de θx,θy en 4 esquinas
    %            cols 9-12 → contribución de Δθs5,Δθs6,Δθs7,Δθs8

    [~, dN_nat] = shapeQ4(xi, eta);
    J    = dN_nat * nodesLoc;
    detJ = det(J);
    dN   = J \ dN_nat;   % [2×4] derivadas físicas de N_i

    % Membrana
    Bm = zeros(3, 8);
    for i = 1:4
        c = (i-1)*2 + 1;
        Bm(:, c:c+1) = [dN(1,i),  0;
                         0,        dN(2,i);
                         dN(2,i),  dN(1,i)];
    end

    % Flexión — esquinas: kappa = [∂θy/∂x; -∂θx/∂y; ∂θy/∂y - ∂θx/∂x]
    % (idéntico a QLLL, Bb de 8 DOFs)
    Bb_corner = zeros(3, 8);
    for i = 1:4
        c = (i-1)*2 + 1;
        Bb_corner(:, c:c+1) = [0,        dN(1,i);
                                -dN(2,i), 0;
                                -dN(1,i), dN(2,i)];
    end

    % Flexión — DOFs jerárquicos de lado (Ec. 6.105 Oñate)
    %
    % Funciones de forma jerárquicas (valor = 1 en nodo de lado):
    %   M5(ξ,η) = (1-ξ²)(1-η)/2   lado s5: η=-1, midpoint (0,-1) → M5=1 ✓
    %   M6(ξ,η) = (1+ξ)(1-η²)/2   lado s6: ξ=+1, midpoint (+1,0) → M6=1 ✓
    %   M7(ξ,η) = (1-ξ²)(1+η)/2   lado s7: η=+1, midpoint (0,+1) → M7=1 ✓
    %   M8(ξ,η) = (1-ξ)(1-η²)/2   lado s8: ξ=-1, midpoint (-1,0) → M8=1 ✓
    %
    % Derivadas naturales [∂/∂ξ; ∂/∂η]:
    dM_nat = zeros(2,4);
    dM_nat(:,1) = [    -xi*(1-eta);    -(1-xi^2)/2 ];   % dM5
    dM_nat(:,2) = [ (1-eta^2)/2;       -eta*(1+xi) ];   % dM6
    dM_nat(:,3) = [    -xi*(1+eta);     (1-xi^2)/2 ];   % dM7
    dM_nat(:,4) = [ -(1-eta^2)/2;      -eta*(1-xi) ];   % dM8

    % Derivadas físicas
    dM = J \ dM_nat;   % [2×4]

    % Δθs_k: rotación tangencial escalar → θ = Mk·Δθs_k·e_side_k
    %   θx = Mk·ex,  θy = Mk·ey
    % Curvaturas (κx=∂θy/∂x, κy=-∂θx/∂y, κxy=∂θy/∂y-∂θx/∂x):
    %   κx  =  ey·∂Mk/∂x
    %   κy  = -ex·∂Mk/∂y
    %   κxy =  ey·∂Mk/∂y - ex·∂Mk/∂x
    Bb_side = zeros(3, 4);
    for s = 1:4
        ex_s = e_side(s,1);
        ey_s = e_side(s,2);
        Bb_side(:,s) = [  ey_s * dM(1,s);
                          -ex_s * dM(2,s);
                           ey_s * dM(2,s) - ex_s * dM(1,s)];
    end

    Bb_exp = [Bb_corner, Bb_side];   % [3×12]
end

function [Bs, detJ] = getMatrixS(nodesLoc, xi, eta)
    % Campo de corte lineal sustituto — idéntico al de QLLL (Oñate §6.7.1)
    % γ = [∂w/∂x + θy;  ∂w/∂y - θx]
    [N, dN_nat] = shapeQ4(xi, eta);
    J    = dN_nat * nodesLoc;
    detJ = det(J);
    dN   = J \ dN_nat;

    Bs = zeros(2, 12);
    for i = 1:4
        c = (i-1)*3 + 1;
        Bs(:, c:c+2) = [dN(1,i),   0,    N(i);
                         dN(2,i), -N(i),  0   ];
    end
end