function [N_el, M_el, Q_el] = RecuperarEsfuerzosQLQL(elem_data, material, geometry, ...
                                                       U_global, dofs_e, xi, eta)
    E  = material.E;
    nu = material.nu;
    t  = geometry.t;

    % Constitutivas
    Dm = (E*t      /(1-nu^2)) * [1,nu,0; nu,1,0; 0,0,(1-nu)/2];
    Db = (E*t^3/12 /(1-nu^2)) * [1,nu,0; nu,1,0; 0,0,(1-nu)/2];
    G  = E/(2*(1+nu));
    Ds = (5/6)*G*t*eye(2);

    T        = elem_data.T;
    nodesLoc = elem_data.nodesLoc;
    e_sides  = elem_data.e_sides;
    RecovMat = elem_data.RecovMat;   % [4×24] en marco local

    % DOFs locales de corner
    d_global = U_global(dofs_e);       % [24×1]
    d_local  = T * d_global;           % [24×1] en marco local del elemento

    % Reconstrucción de DOFs internos Δθs
    d_int = RecovMat * d_local;        % [4×1]

    % Cinemática en (xi, eta)
    [N, dN_nat]     = shapeQ4_loc(xi, eta);
    J               = dN_nat * nodesLoc;
    dN              = J \ dN_nat;
    [phi, dphi_nat] = bubbleFuncs(xi, eta);
    dphi            = J \ dphi_nat;

    % Membrana Bm [3×8]
    Bm = zeros(3,8);
    for k = 1:4
        c = (k-1)*2+1;
        Bm(:,c:c+1) = [dN(1,k),0; 0,dN(2,k); dN(2,k),dN(1,k)];
    end
    N_el = Dm * Bm * d_local([1,2,7,8,13,14,19,20]);

    % Flexión Bb [3×12]
    Bb = zeros(3,12);
    for k = 1:4
        c = (k-1)*2+1;
        Bb(:,c:c+1) = [0,dN(1,k); -dN(2,k),0; -dN(1,k),dN(2,k)];
    end
    for k = 1:4
        ex = e_sides(k,1); ey = e_sides(k,2);
        c  = 8+k;
        Bb(1,c) =  dphi(1,k)*ey;
        Bb(2,c) = -dphi(2,k)*ex;
        Bb(3,c) = -dphi(1,k)*ex + dphi(2,k)*ey;
    end
    M_el = Db * Bb * [d_local([4,5,10,11,16,17,22,23]); d_int];

    % Corte Bs [2×16]
    Bs = zeros(2,16);
    for k = 1:4
        c = (k-1)*3+1;
        Bs(:,c:c+2) = [dN(1,k),0,N(k); dN(2,k),-N(k),0];
    end
    for k = 1:4
        ex = e_sides(k,1); ey = e_sides(k,2);
        c  = 12+k;
        Bs(1,c) =  phi(k)*ey;
        Bs(2,c) = -phi(k)*ex;
    end
    Q_el = Ds * Bs * [d_local([3,4,5,9,10,11,15,16,17,21,22,23]); d_int];
end


%% ================================================================ %%
%% Funciones locales (copias privadas para independencia del archivo)
%% ================================================================ %%

function [N,dN] = shapeQ4_loc(xi,eta)
    N  = 0.25*[(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)];
    dN = 0.25*[-(1-eta),(1-eta),(1+eta),-(1+eta);
               -(1-xi),-(1+xi),(1+xi),(1-xi)];
end

function [phi,dphi_nat] = bubbleFuncs(xi,eta)
    phi      = zeros(1,4);
    dphi_nat = zeros(2,4);
    phi(1)        = (1-xi^2)*(1-eta);
    dphi_nat(1,1) = -2*xi*(1-eta);
    dphi_nat(2,1) = -(1-xi^2);
    phi(2)        = (1-eta^2)*(1+xi);
    dphi_nat(1,2) =  (1-eta^2);
    dphi_nat(2,2) = -2*eta*(1+xi);
    phi(3)        = (1-xi^2)*(1+eta);
    dphi_nat(1,3) = -2*xi*(1+eta);
    dphi_nat(2,3) =  (1-xi^2);
    phi(4)        = (1-eta^2)*(1-xi);
    dphi_nat(1,4) = -(1-eta^2);
    dphi_nat(2,4) = -2*eta*(1-xi);
end