function [Kele, KeleB, KeleS, rowIndex, columnIndex] = getStiffnessMatrix_Plates(nodesEle, eleDofs, Cb, Cs, integrationMethod)
% getStiffnessMatrix_Plates  Mindlin plate element stiffness matrix
%
%   [Kele, KeleB, KeleS, rowIndex, columnIndex] = ...
%       getStiffnessMatrix_Plates(nodesEle, eleDofs, Cb, Cs, integrationMethod)
%
%   Inputs:
%     nodesEle          - (nNodEle x 2) physical node coordinates [x, y]
%     eleDofs           - (nDofEle x 1) global DOF indices for this element
%     Cb                - (3x3) bending constitutive matrix (getElementConstitutive)
%     Cs                - (2x2) shear constitutive matrix  (getElementConstitutive)
%     integrationMethod - 'full'      : 2x2 bending + 2x2 shear
%                         'selective' : 2x2 bending + 1x1 shear  (avoids shear locking)
%                         'reduced'   : 1x1 bending + 1x1 shear
%
%   Outputs:
%     Kele        - (nDofEle x nDofEle) total element stiffness  (KeleB + KeleS)
%     KeleB       - (nDofEle x nDofEle) bending contribution
%     KeleS       - (nDofEle x nDofEle) shear contribution
%     rowIndex    - (nDofEle x nDofEle) row    DOF index matrix (for sparse assembly)
%     columnIndex - (nDofEle x nDofEle) column DOF index matrix (for sparse assembly)
%
%   DOF order per node: [w, phi_x, phi_y]  (3 DOF/node)
%   Hardcoded for Q4 elements (bilinear shape functions).

switch integrationMethod
    case 'full'
        npgBending = 2;
        npgShear   = 2;
    case 'selective'
        npgBending = 2;
        npgShear   = 1;
    case 'reduced'
        npgBending = 1;
        npgShear   = 1;
    otherwise
        error('Integration method ''%s'' not supported.', integrationMethod);
end

nDofEle = numel(eleDofs);

KeleB = zeros(nDofEle);
KeleS = zeros(nDofEle);

% --- Gauss quadratures (called once, outside the loops) ---
[wpgB, upgB, npgB] = gauss(repmat(npgBending, 2, 1));
[wpgS, upgS, npgS] = gauss(repmat(npgShear,   2, 1));

% --- Shape function derivatives at all Gauss points (precomputed) ---
dNB = shapefunsDer(upgB, 'Q4');   % (2 x nNodEle x npgB)
dNS = shapefunsDer(upgS, 'Q4');   % (2 x nNodEle x npgS)
NiS = shapefuns(upgS, 'Q4');      % (1 x nNodEle x npgS)

% --- Bending contribution ---
for ipgB = 1:npgB

    jac    = dNB(:,:,ipgB) * nodesEle(:,1:2);   % 2x2 Jacobian
    dNBxy  = jac \ dNB(:,:,ipgB);               % 2 x nNodEle  (physical derivs)

    Bb = zeros(3, nDofEle);
    Bb(1, 2:3:end) = dNBxy(1,:);   % kappa_x  = d(phi_x)/dx
    Bb(2, 3:3:end) = dNBxy(2,:);   % kappa_y  = d(phi_y)/dy
    Bb(3, 2:3:end) = dNBxy(2,:);   % kappa_xy = d(phi_x)/dy + d(phi_y)/dx
    Bb(3, 3:3:end) = dNBxy(1,:);

    KeleB = KeleB + Bb' * Cb * Bb * wpgB(ipgB) * det(jac);
end

% --- Shear contribution ---
for ipgS = 1:npgS

    jac   = dNS(:,:,ipgS) * nodesEle(:,1:2);    % 2x2 Jacobian
    dNSxy = jac \ dNS(:,:,ipgS);                % 2 x nNodEle

    Bs = zeros(2, nDofEle);
    Bs(1, 1:3:end) =  dNSxy(1,:);    % gamma_xz: dw/dx
    Bs(2, 1:3:end) =  dNSxy(2,:);    % gamma_yz: dw/dy
    Bs(1, 2:3:end) = -NiS(1,:,ipgS); % gamma_xz: -phi_x
    Bs(2, 3:3:end) = -NiS(1,:,ipgS); % gamma_yz: -phi_y

    KeleS = KeleS + Bs' * Cs * Bs * wpgS(ipgS) * det(jac);
end

Kele = KeleB + KeleS;

% --- Index matrices for sparse assembly (see assembleSystem) ---
rowIndex    = repmat(eleDofs,  1, size(eleDofs,1));
columnIndex = repmat(eleDofs', size(eleDofs,1), 1);

end
