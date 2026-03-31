function B = getBAHMAD_Local(Ni, dN_nat, invJ, R_gauss, zeta, tEle, v1Ele, v2Ele)
    nnodele = length(Ni);
    B = zeros(5, nnodele * 5);
    
    % 1. Derivadas cartesianas globales dN/dX
    dN_glob = invJ * [dN_nat; zeros(1, nnodele)]; 
    
    % 2. Terna local en el punto de Gauss (l, m, n son e1, e2, e3)
    l = R_gauss(1,:); m = R_gauss(2,:); n = R_gauss(3,:);
    
    % 3. Inversa del Jacobiano rotada para d(zeta)/dX'
    invJ_loc = R_gauss * invJ;

    for i = 1:nnodele
        % --- PARTE TRASLACIONAL (u, v, w) ---
        dNi_dxloc = l * dN_glob(:,i);
        dNi_dyloc = m * dN_glob(:,i);
        dNi_dzloc = n * dN_glob(:,i); 
        
        idx = (i-1)*5;
        B(1, idx+(1:3)) = l * dNi_dxloc;
        B(2, idx+(1:3)) = m * dNi_dyloc;
        B(3, idx+(1:3)) = l * dNi_dyloc + m * dNi_dxloc;
        B(4, idx+(1:3)) = l * dNi_dzloc + n * dNi_dxloc;
        B(5, idx+(1:3)) = m * dNi_dzloc + n * dNi_dyloc;

        % --- PARTE ROTACIONAL (alpha, beta) ---
        % f = brazo de palanca (zeta * t/2)
        f = 0.5 * tEle(i);
        
        % Gradientes de (Ni * zeta) en ejes locales
        dN_zeta_dx = dNi_dxloc * zeta + Ni(i) * invJ_loc(1,3);
        dN_zeta_dy = dNi_dyloc * zeta + Ni(i) * invJ_loc(2,3);
        dN_zeta_dz = dNi_dzloc * zeta + Ni(i) * invJ_loc(3,3);

        % Vectores directores nodales
        V1 = v1Ele(:,i); V2 = v2Ele(:,i);
        
        % PROYECCIONES CRUCIALES:
        % Alpha (Giro sobre V2) -> Produce desplazamientos en V1
        v1l = l*V1; 
        v1m = m*V1;
        v1n = n*V1;

        % Beta (Giro sobre V1) -> Produce desplazamientos en -V2
        v2l = l*V2; 
        v2m = m*V2; 
        v2n = n*V2;

        % Llenado de B para alpha (columna 4)
        B(1, idx+4) = f * (v1l * dN_zeta_dx);
        B(2, idx+4) = f * (v1m * dN_zeta_dy);
        B(3, idx+4) = f * (v1l * dN_zeta_dy + v1m * dN_zeta_dx);
        B(4, idx+4) = f * (v1l * dN_zeta_dz + v1n * dN_zeta_dx);
        B(5, idx+4) = f * (v1m * dN_zeta_dz + v1n * dN_zeta_dy);

        % Llenado de B para beta (columna 5) con signo NEGATIVO por convención horaria
        B(1, idx+5) = f * (-v2l * dN_zeta_dx);
        B(2, idx+5) = f * (-v2m * dN_zeta_dy);
        B(3, idx+5) = f * (-v2l * dN_zeta_dy - v2m * dN_zeta_dx);
        B(4, idx+5) = f * (-v2l * dN_zeta_dz - v2n * dN_zeta_dx);
        B(5, idx+5) = f * (-v2m * dN_zeta_dz - v2n * dN_zeta_dy);
    end
end