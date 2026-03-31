function [w, g, ng] = gaussAHMAD(n_puntos_plano, n_puntos_espesor)
    % n_puntos_plano: ej. 3 para 3x3
    % n_puntos_espesor: ej. 2 para zeta
    
    % Gauss 1D para obtener bases
    [w1, g1] = gauss_base(n_puntos_plano);
    [wZ, gZ] = gauss_base(n_puntos_espesor);
    
    ng = (n_puntos_plano^2) * n_puntos_espesor;
    w = zeros(ng, 1);
    g = zeros(ng, 3);
    
    count = 0;
    for i = 1:n_puntos_plano
        for j = 1:n_puntos_plano
            for k = 1:n_puntos_espesor
                count = count + 1;
                g(count, :) = [g1(i), g1(j), gZ(k)];
                w(count) = w1(i) * w1(j) * wZ(k);
            end
        end
    end
end

function [w, g] = gauss_base(n)
    switch n
        case 1
            g = 0; w = 2;
        case 2
            g = [-1, 1]/sqrt(3); w = [1, 1];
        case 3
            g = [-sqrt(0.6), 0, sqrt(0.6)]; w = [5/9, 8/9, 5/9];
    end
end