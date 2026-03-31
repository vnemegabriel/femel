function [Ke_global] = CalcularRigidezQLQL(nodes3D, geometry, material)
% CalcularRigidezQLQL  Elemento shell plano QLQL

    E = material.E;  nu = material.nu;  t = geometry.t;
    G = E / (2*(1+nu));

    % 1. Sistema Local Consistente
    [T_full, nodesLoc] = getTransformationMatrix(nodes3D);
    x = nodesLoc(:,1); y = nodesLoc(:,2);

    % 2. Matrices Constitutivas
    Dm = (E*t / (1-nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    Db = (E*t^3 / (12*(1-nu^2))) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    Ds = (5/6) * G * t * eye(2);

    Ke_local = zeros(24, 24);

    % 3. Integración Numérica (2x2 Gauss Completa)
    gp = [-1/sqrt(3), 1/sqrt(3)];
    w_gauss = [1, 1];

    Bg1_1 = zeros(1, 24); Bg1_3 = zeros(1, 24);
    Bg2_2 = zeros(1, 24); Bg2_4 = zeros(1, 24);

    % Tying point 1 (Lado 1-2: xi=0, eta=-1)
    Bg1_1(3) = -0.5; Bg1_1(9) = 0.5;
    Bg1_1(4) = -0.25*(y(2)-y(1)); Bg1_1(10) = -0.25*(y(2)-y(1));
    Bg1_1(5) =  0.25*(x(2)-x(1)); Bg1_1(11) =  0.25*(x(2)-x(1));

    % Tying point 3 (Lado 4-3: xi=0, eta=1) -> Sentido de 4 a 3
    Bg1_3(21) = -0.5; Bg1_3(15) = 0.5;
    Bg1_3(22) = -0.25*(y(3)-y(4)); Bg1_3(16) = -0.25*(y(3)-y(4));
    Bg1_3(23) =  0.25*(x(3)-x(4)); Bg1_3(17) =  0.25*(x(3)-x(4));

    % Tying point 2 (Lado 2-3: xi=1, eta=0) -> Sentido de 2 a 3
    Bg2_2(9) = -0.5; Bg2_2(15) = 0.5;
    Bg2_2(10) = -0.25*(y(3)-y(2)); Bg2_2(16) = -0.25*(y(3)-y(2));
    Bg2_2(11) =  0.25*(x(3)-x(2)); Bg2_2(17) =  0.25*(x(3)-x(2));

    % Tying point 4 (Lado 1-4: xi=-1, eta=0) -> Sentido de 1 a 4
    Bg2_4(3) = -0.5; Bg2_4(21) = 0.5;
    Bg2_4(4) = -0.25*(y(4)-y(1)); Bg2_4(22) = -0.25*(y(4)-y(1));
    Bg2_4(5) =  0.25*(x(4)-x(1)); Bg2_4(23) =  0.25*(x(4)-x(1));

    [N0, dN_nat0] = shapeQ4(0, 0);
    J0 = dN_nat0 * nodesLoc;
    detJ0 = det(J0);
    dN_xy0 = J0 \ dN_nat0;

    % Ensamblaje sobre puntos de Gauss
    for i = 1:2
        for j = 1:2
            xi = gp(i);
            eta = gp(j);
            
            [~, dN_nat] = shapeQ4(xi, eta);
            J = dN_nat * nodesLoc;
            detJ = det(J);
            dN_xy = J \ dN_nat;
            
            Bm = zeros(3, 24);
            Bb = zeros(3, 24);
            
            for k = 1:4
                col = 6*(k-1);
                
                % Componentes normales evaluadas normalmente
                Bm(1, col+1) = dN_xy(1,k);
                Bm(2, col+2) = dN_xy(2,k);
                % Componente de corte evaluada en el centroide 
                Bm(3, col+1) = dN_xy0(2,k); 
                Bm(3, col+2) = dN_xy0(1,k); 
                
                % Flexión (Completa)
                Bb(1, col+5) =  dN_xy(1,k);
                Bb(2, col+4) = -dN_xy(2,k);
                Bb(3, col+5) =  dN_xy(2,k); 
                Bb(3, col+4) = -dN_xy(1,k);
            end
            
            % Corte Transversal: Campo sustituto interpolado MITC4
            Bg1 = 0.5*(1 - eta)*Bg1_1 + 0.5*(1 + eta)*Bg1_3;
            Bg2 = 0.5*(1 - xi)*Bg2_4  + 0.5*(1 + xi)*Bg2_2;
            Bg = [Bg1; Bg2];
            
            Bs = J \ Bg;
            
            factor = w_gauss(i) * w_gauss(j) * detJ;
            Ke_local = Ke_local + (Bm' * Dm * Bm) * factor;
            Ke_local = Ke_local + (Bb' * Db * Bb) * factor;
            Ke_local = Ke_local + (Bs' * Ds * Bs) * factor;
        end
    end

    % 4. RIGIDEZ FICTICIA DE DRILLING CONSISTENTE
    alpha_drill = 10e-3; % Valor recomendado por Oñate
    
    B_drill = zeros(1, 24);
    for k = 1:4
        col = 6*(k-1);
        % Acoplamiento con la rotación macroscópica de membrana en el centroide
        B_drill(col+1) =  0.5 * dN_xy0(2,k); 
        B_drill(col+2) = -0.5 * dN_xy0(1,k); 
        B_drill(col+6) =  N0(k);             
    end
    
    K_drill_mat = (alpha_drill * G * t) * (B_drill' * B_drill) * (4 * detJ0);
    Ke_local = Ke_local + K_drill_mat;

    % 5. Transformación al sistema global
    Ke_global = T_full' * Ke_local * T_full;
end

function [N, dN_nat] = shapeQ4(xi, eta)
    N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    dN_nat = 0.25 * [-(1-eta),  (1-eta),  (1+eta), -(1+eta);
                     -(1-xi),  -(1+xi),   (1+xi),   (1-xi)];
end