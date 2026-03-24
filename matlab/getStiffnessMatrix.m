function K = getStiffnessMatrix(nodesElem,eleType,t,E,nPointsGauss) 

npg_xi = nPointsGauss(1);
npg_eta = nPointsGauss(2);

[weights_xi,lgp_xi] = gauss1D(npg_xi);
[weights_eta,lgp_eta] = gauss1D(npg_eta);

nDofNod = 2; %2D
nNodEle = size(nodesElem,1);
nDofEle = nDofNod*nNodEle;

switch eleType

    case 'Q4'
        K = zeros(nDofEle);

        for ipg_xi = 1:npg_xi
            for ipg_eta = 1:npg_eta

            xi = lgp_xi(ipg_xi);
            eta = lgp_eta(ipg_eta);

            dN_iso = shapefunsDer([xi eta],'Q4');
            
            Jacobian = dN_iso*nodesElem;

            dNxy = Jacobian\dN_iso; %Jacobian is square matrix --> dim(dNxy) ==  dim(dN_iso)

%-----------OLD -- Doesnt compute the inverse directly

            % dx_dxi = dN_dxi*nodes(:,1);
            % dx_deta = dN_deta*nodes(:,1);
            % dy_dxi = dN_dxi*nodes(:,2);
            % dy_deta = dN_deta*nodes(:,2);

            % Jacobian = [dx_dxi dy_dxi;
            %           dx_deta dy_deta];
         

            % J_inv =inv(Jacobian);
            % 
            % J_hat = [J_inv zeros(2);
            %          zeros(2) J_inv];
            % 
            % A_0 = [1 0 0 0;
            %        0 0 0 1;
            %        0 1 1 0]; %Defines coefficient of derivatives that compose the strain
            % 
            % A = A_0 * J_hat; % 3X4
            % 
            % G = [dN1_dxi 0 dN2_dxi 0 dN3_dxi 0 dN4_dxi 0;
            %  dN1_deta 0 dN2_deta 0 dN3_deta 0 dN4_deta 0;
            %  0 dN1_dxi 0 dN2_dxi 0 dN3_dxi 0 dN4_dxi;
            %  0 dN1_deta 0 dN2_deta 0 dN3_deta 0 dN4_deta]; % 4X8

            % B = A*G;

%-----------OLD 
        
            %Strain displacement B matrix 
            B = zeros(3,nDofEle);
            B(1, 1:2:end) = dNxy(1,:);
            B(2, 2:2:end) = dNxy(2,:);
            B(3, 1:2:end) = dNxy(2,:);
            B(3, 2:2:end) = dNxy(1,:);

            %Integrand
            K = K + weights_xi(ipg_xi)*weights_eta(ipg_eta)*t*B'*E*B*det(Jacobian); 
            end
        end
    case 'Q8'
        K = zeros(nDofEle);

        for ipg_xi = 1:npg_xi
            for ipg_eta = 1:npg_eta

            xi = lgp_xi(ipg_xi);
            eta = lgp_eta(ipg_eta);

            dN_iso = shapefunsDer([xi eta],'Q8');
            
            Jacobian = dN_iso*nodesElem;

            dNxy = Jacobian\dN_iso;

            %Strain displacement B matrix 
            B = zeros(3,nDofEle);
            B(1, 1:2:end) = dNxy(1,:);
            B(2, 2:2:end) = dNxy(2,:);
            B(3, 1:2:end) = dNxy(2,:);
            B(3, 2:2:end) = dNxy(1,:);

            %Integrand
            K = K + weights_xi(ipg_xi)*weights_eta(ipg_eta)*t*B'*E*B*det(Jacobian); 
            end
        end
end

end


