function [N, M, Q] = CalcularEsfuerzosCompletos(mesh, geometry, material, D_global, evalPoints)
    % Recuperación de esfuerzos Nx, My, Qy para elementos de cáscara plana
    nEle = size(mesh.connect, 1);
    nEval = size(evalPoints, 2);
    t = geometry.t; 
    E = material.E; 
    nu = material.nu;
    
    % Matrices constitutivas
    Dm = (E*t/(1-nu^2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
    Db = (E*t^3/(12*(1-nu^2))) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
    G = E / (2*(1+nu));
    Ds = (5/6) * G * t * eye(2);
    
    N = zeros(3, nEval, nEle); 
    M = zeros(3, nEval, nEle); 
    Q = zeros(2, nEval, nEle);
    
    dNNodes = shapefunsder(evalPoints, 'Q4');
    
    for iEle = 1:nEle
        % 1. Obtener la matriz de transformación del elemento (YA ES 24x24)
        nodesEle3D = mesh.nodes(mesh.connect(iEle, :), :);
        [Te, nodesLoc2D] = getTransformationMatrix(nodesEle3D);
        
        % 2. Extraer desplazamientos globales del elemento (24x1)
        nodes_e = mesh.connect(iEle, :);
        idx_global = zeros(24, 1);
        for i = 1:4
            idx_global((i-1)*6 + (1:6)) = (nodes_e(i)-1)*6 + (1:6);
        end
        
        % Aseguramos que sea un vector columna antes de multiplicar
        Ue_global = D_global(idx_global);
        Ue_global = Ue_global(:); 
        
        % 3. Transformar a local (24x24 * 24x1 = 24x1)
        Ue_local = Te * Ue_global;
        
        for p = 1:nEval
            jac = dNNodes(:,:,p) * nodesLoc2D;
            dNBxy = jac \ dNNodes(:,:,p);
            
            [Ni, ~] = shapefuns(evalPoints(:,p), 'Q4');
            
            Bm = zeros(3, 24); Bb = zeros(3, 24); Bs = zeros(2, 24);
            
            for i = 1:4
                col = (i-1)*6;
                % Membrana
                Bm(1, col+1) = dNBxy(1,i); 
                Bm(2, col+2) = dNBxy(2,i);
                Bm(3, col+1) = dNBxy(2,i); 
                Bm(3, col+2) = dNBxy(1,i);
                
                % Flexión (kappa = [dth_y/dx; -dth_x/dy; ...])
                Bb(1, col+5) =  dNBxy(1,i); 
                Bb(2, col+4) = -dNBxy(2,i);
                Bb(3, col+5) =  dNBxy(2,i); 
                Bb(3, col+4) = -dNBxy(1,i);
                
                % Corte (gamma = [dw/dx + th_y; dw/dy - th_x])
                Bs(1, col+3) = dNBxy(1,i); 
                Bs(1, col+5) = Ni(i);      
                Bs(2, col+3) = dNBxy(2,i); 
                Bs(2, col+4) = -Ni(i);     
            end
            
            N(:,p,iEle) = Dm * Bm * Ue_local;
            M(:,p,iEle) = Db * Bb * Ue_local;
            Q(:,p,iEle) = Ds * Bs * Ue_local;
        end
    end
end

function [Ni, numNodos] = shapefuns(upg, elementType)
    % Evaluamos N para el elemento Q4 en cada punto de Gauss
    if strcmpi(elementType, 'Q4')
        numNodos = 4;
        npg = size(upg, 2);
        Ni = zeros(numNodos, 1, npg);
        
        for p = 1:npg
            xi = upg(1, p);
            eta = upg(2, p);
            % N_i = 1/4 * (1 +- xi) * (1 +- eta)
            Ni(:, 1, p) = 0.25 * [ (1 - xi)*(1 - eta);
                                   (1 + xi)*(1 - eta);
                                   (1 + xi)*(1 + eta);
                                   (1 - xi)*(1 + eta) ];
        end
    else
        error('Elemento no soportado');
    end
end

function dN = shapefunsder(upg, elementType)
    % Evaluamos las derivadas naturales de N para el elemento Q4
    if strcmpi(elementType, 'Q4')
        npg = size(upg, 2);
        % dN tendrá tamaño (2 derivadas, 4 nodos, npg puntos de Gauss)
        dN = zeros(2, 4, npg);
        
        for p = 1:npg
            xi = upg(1, p);
            eta = upg(2, p);
            
            % dN/dxi
            dN(1, :, p) = 0.25 * [ -(1 - eta), (1 - eta), (1 + eta), -(1 + eta) ];
            % dN/deta
            dN(2, :, p) = 0.25 * [ -(1 - xi), -(1 + xi), (1 + xi),  (1 - xi) ];
        end
    else
        error('Elemento no soportado');
    end
end