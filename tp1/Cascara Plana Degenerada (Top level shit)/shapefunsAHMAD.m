function [Ni,N] = shapefunsAHMAD(pointArray,eleType)  
    ngauss = size(pointArray, 1);
    Id = eye(3);

    switch upper(eleType)
    
        case 'AHMAD9'
            N  = zeros(3,3*9,ngauss);
            Ni = zeros(1,9,ngauss);

            for igauss = 1:ngauss
                ksi = pointArray(igauss,1);
                eta = pointArray(igauss,2);
    
                N9 =      (1 - ksi^2)*(1 - eta^2);
                N8 = 0.50*(1 - ksi  )*(1 - eta^2) - 0.5*N9;
                N7 = 0.50*(1 - ksi^2)*(1 + eta  ) - 0.5*N9;
                N6 = 0.50*(1 + ksi  )*(1 - eta^2) - 0.5*N9;
                N5 = 0.50*(1 - ksi^2)*(1 - eta  ) - 0.5*N9;
                N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8 + 0.5*N9);
                N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7 + 0.5*N9);
                N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6 + 0.5*N9);
                N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8 + 0.5*N9);
    
                Ni(1,:,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8 N9];
                N (:,:,igauss) = [ N1*Id, N2*Id, N3*Id, ...
                    N4*Id, N5*Id, N6*Id, ...
                    N7*Id, N8*Id, N9*Id ];
            end
    
        case 'AHMAD8'
            N  = zeros(3,3*8,ngauss);
            Ni = zeros(1,  8,ngauss);

            for igauss = 1:ngauss
                ksi = pointArray(igauss,1);
                eta = pointArray(igauss,2);
    
                N8 = 0.50*(1 - ksi  )*(1 - eta^2);
                N7 = 0.50*(1 - ksi^2)*(1 + eta  );
                N6 = 0.50*(1 + ksi  )*(1 - eta^2);
                N5 = 0.50*(1 - ksi^2)*(1 - eta  );
                N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8);
                N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7);
                N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6);
                N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8);
    
                Ni(1,:,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8];
                N (:,:,igauss) = [ N1*Id, N2*Id, N3*Id, ...
                    N4*Id, N5*Id, N6*Id, ...
                    N7*Id, N8*Id ];
            end
    
        case 'AHMAD4'
            N  = zeros(3,3*4,ngauss);
            Ni = zeros(1,  4,ngauss);

            for igauss = 1:ngauss
                ksi = pointArray(igauss,1);
                eta = pointArray(igauss,2);
    
                N4 = 0.25*(1 - ksi)*(1 + eta);
                N3 = 0.25*(1 + ksi)*(1 + eta);
                N2 = 0.25*(1 + ksi)*(1 - eta);
                N1 = 0.25*(1 - ksi)*(1 - eta);
    
                Ni(1,:,igauss) = [N1 N2 N3 N4];
                N (:,:,igauss) = [ N1*Id, N2*Id, N3*Id, N4*Id ];
            end
    
    end
end