function [Rele, eleDofs] = applyLoadOnElement(nodesEle, eleDofs, fz, mx, my, npg)
% applyLoadOnElement  Element load vector for Mindlin plate distributed load
%
%   [Rele, eleDofs] = applyLoadOnElement(nodesEle, eleDofs, fz, mx, my, npg)
%
%   Inputs:
%     nodesEle - (nNodEle x 2) physical node coordinates [x, y]
%     eleDofs  - (nDofEle x 1) global DOF indices for this element
%     fz       - (nNodEle x 1) nodal transverse load values  [N/m^2]
%     mx       - (nNodEle x 1) nodal moment load about x     [N/m]
%     my       - (nNodEle x 1) nodal moment load about y     [N/m]
%     npg      - number of Gauss points per direction (scalar)
%
%   Outputs:
%     Rele    - (nDofEle x 1) element load vector
%     eleDofs - (nDofEle x 1) global DOF indices (passed through)
%
%   The distributed load is approximated via shape functions:
%     f_e = int( N^T * N * q ) dA
%   where N is (3 x nDofEle) and q = [fz1;mx1;my1; fz2;mx2;my2; ...] (nDofEle x 1)
%
%   DOF order per node: [w, phi_x, phi_y]
%   Hardcoded for Q4 elements.

[wpg, upg, npg] = gauss(repmat(npg, 2, 1));  % overwrites npg with total point count

Ni = shapefuns(upg,    'Q4');   % (1 x nNodEle x npg)  precomputed
dN = shapefunsDer(upg, 'Q4');   % (2 x nNodEle x npg)  precomputed

nDofEle = numel(eleDofs);
Rele    = zeros(nDofEle, 1);

q = reshape([fz, mx, my]', [], 1);   % unwrap 3 column vectors -> (nDofEle x 1)

for ipg = 1:npg

    jac = dN(:,:,ipg) * nodesEle(:,1:2);   % 2x2 Jacobian

    N = zeros(3, nDofEle);
    N(1, 1:3:end) = Ni(1,:,ipg);   % fz  -> w     DOFs
    N(2, 2:3:end) = Ni(1,:,ipg);   % mx  -> phi_x DOFs
    N(3, 3:3:end) = Ni(1,:,ipg);   % my  -> phi_y DOFs

    Rele = Rele + N' * N * q * wpg(ipg) * det(jac);
end

end
