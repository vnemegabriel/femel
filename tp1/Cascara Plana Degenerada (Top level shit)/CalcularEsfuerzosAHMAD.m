function [Nx, My, Qy] = CalcularEsfuerzosAHMAD(phi_vec, mesh, U, nodalSys, material, geometria, eleType)
    t = geometria.t;
    D = material.D;
    n_pts = length(phi_vec);
    Nx = zeros(1, n_pts); My = zeros(1, n_pts); Qy = zeros(1, n_pts);
    
    wz = [1, 1]; 
    gz = [-1/sqrt(3), 1/sqrt(3)]; 
    
    for p = 1:n_pts
        phi_target = phi_vec(p);
        y_target = geometria.L_half; 
        
        iele = 0;
        for e = 1:size(mesh.connect, 1)
            nodes_e = mesh.connect(e, :);
            coords = mesh.nodes(nodes_e, :);
            if any(coords(:,2) > y_target - 1e-3)
                phi_nodes = rad2deg(atan2(coords(:,1), coords(:,3)));
                if phi_target >= (min(phi_nodes) - 1e-3) && phi_target <= (max(phi_nodes) + 1e-3)
                    iele = e; break; 
                end
            end
        end
        if iele == 0, continue; end
        
        nodes_e = mesh.connect(iele, :);
        coords = mesh.nodes(nodes_e, :);
        phi_nodes = rad2deg(atan2(coords(:,1), coords(:,3)));
        
        % --- CORRECCIÓN CRÍTICA DE MAPEO GEOMÉTRICO ---
        % PHI corresponde a la dirección eta (Circunferencial)
        eta = 2 * (phi_target - min(phi_nodes)) / (max(phi_nodes) - min(phi_nodes)) - 1;
        
        % Y=25 corresponde a la dirección ksi (Longitudinal)
        ksi_f = 1.0; % Borde físico para Flexión/Membrana
        
        % Puntos de Barlow longitudinales para evitar Corte Parásito
        if strcmp(eleType, 'AHMAD4')
            ksi_c = 0.0; 
        else
            ksi_c = 1 / sqrt(3); 
        end
        
        v1_e = squeeze(nodalSys(:, 1, nodes_e));
        v2_e = squeeze(nodalSys(:, 2, nodes_e));
        v3_e = squeeze(nodalSys(:, 3, nodes_e));
        t_vec = ones(1, length(nodes_e)) * t;
        Ue = U(reshape(bsxfun(@plus, (nodes_e-1)*5, (1:5)'), 1, []));
        
        % --- INTEGRACIÓN EN EL ESPESOR ---
        for igz = 1:length(gz)
            zeta = gz(igz);
            
            % 1. Flexión y Membrana (En ksi_f)
            [Ni_f, ~] = shapefunsAHMAD([ksi_f, eta], eleType);
            dN_f = shapefunsderAHMAD([ksi_f, eta], eleType);
            jac_f = shelljac(Ni_f, dN_f(:,:,1), zeta, coords, t_vec, v3_e);
            
            e1f = jac_f(1,:)/norm(jac_f(1,:));
            e3f = cross(jac_f(1,:), jac_f(2,:))/norm(cross(jac_f(1,:), jac_f(2,:)));
            B_f = getBAHMAD_Local(Ni_f, dN_f(:,:,1), inv(jac_f), [e1f; cross(e3f, e1f); e3f], zeta, t_vec, v1_e, v2_e);
            sigma_f = D * B_f * Ue;
            
            % 2. Corte Transversal (En ksi_c, Punto de Barlow)
            [Ni_c, ~] = shapefunsAHMAD([ksi_c, eta], eleType);
            dN_c = shapefunsderAHMAD([ksi_c, eta], eleType);
            jac_c = shelljac(Ni_c, dN_c(:,:,1), zeta, coords, t_vec, v3_e);
            
            e1c = jac_c(1,:)/norm(jac_c(1,:));
            e3c = cross(jac_c(1,:), jac_c(2,:))/norm(cross(jac_c(1,:), jac_c(2,:)));
            B_c = getBAHMAD_Local(Ni_c, dN_c(:,:,1), inv(jac_c), [e1c; cross(e3c, e1c); e3c], zeta, t_vec, v1_e, v2_e);
            sigma_c = D * B_c * Ue;
            
            dz = t/2; 
            Nx(p) = Nx(p) + sigma_f(1) * dz * wz(igz);
            My(p) = My(p) - sigma_f(2) * (zeta * t/2) * dz * wz(igz);
            
            % MAGIA PURA: Ahora sigma_c(3) recorrerá el arco correctamente y sin serrucho
            Qy(p) = Qy(p) + sigma_c(3) * dz * wz(igz); 
        end
    end
end