function [M, Q] = getCharacteristicStresses(nodes, elements, properties, D, evaluationPoints)
% getCharacteristicStresses  Mindlin plate stress resultants at evaluation points
%
%   [M, Q] = getCharacteristicStresses(nodes, elements, properties, D, evaluationPoints)
%
%   Inputs:
%     nodes           - (nNod x 2)       global node coordinates [x, y]
%     elements        - struct with fields:
%                         .connectivity  (nEle x nNodEle) element node IDs (1-based)
%                         .property      (nEle x 1)       property index per element
%                         .dof           (nEle x nDofEle) global DOF indices per element
%     properties      - struct array, each entry has:
%                         .Cb  (3x3) bending constitutive matrix
%                         .Cs  (2x2) shear constitutive matrix
%     D               - (nDofTot x 1)   global displacement vector
%     evaluationPoints- (nEvalPoints x 2) isoparametric coordinates [xi, eta]
%                       e.g. [0 0] for centroid, or Gauss points for superconvergence
%
%   Outputs:
%     M  - (3 x nEvalPoints x nEle)  bending moment resultants
%            M(1,:,:) = Mx,   M(2,:,:) = My,   M(3,:,:) = Mxy   [N*m/m]
%     Q  - (2 x nEvalPoints x nEle)  shear force resultants
%            Q(1,:,:) = Qx,   Q(2,:,:) = Qy                      [N/m]
%
%   Theory:  sigma_hat = D_hat * B * a^(e)
%     M = Cb * Bb * De    (bending curvatures -> moments)
%     Q = Cs * Bs * De    (transverse shear strains -> shear forces)
%
%   Shape functions are precomputed at all evaluation points before the
%   element loop (same pattern as getStiffnessMatrix_Plates).

nEle    = size(elements.connectivity, 1);
nDofNod = size(elements.dof, 2) / size(elements.connectivity, 2);  % DOF per node
nNodEle = size(elements.connectivity, 2);
nDofEle = nNodEle * nDofNod;

nEvalPoints = size(evaluationPoints, 1);

% --- Precompute shape functions at all evaluation points ---
dNNodes = shapefunsDer(evaluationPoints, 'Q4');   % (2 x nNodEle x nEvalPoints)
NiNodes = shapefuns(evaluationPoints,    'Q4');   % (1 x nNodEle x nEvalPoints)

M = zeros(3, nEvalPoints, nEle);
Q = zeros(2, nEvalPoints, nEle);

for iEle = 1:nEle

    Cb = properties(elements.property(iEle)).Cb;
    Cs = properties(elements.property(iEle)).Cs;

    eleDofs  = elements.dof(iEle,:)';
    nodesEle = nodes(elements.connectivity(iEle,:), :);
    De       = D(eleDofs);

    for iEval = 1:nEvalPoints

        jac   = dNNodes(:,:,iEval) * nodesEle(:,1:2);   % 2x2
        dNBxy = jac \ dNNodes(:,:,iEval);               % 2 x nNodEle

        Bb = zeros(3, nDofEle);
        Bb(1, 2:3:end) = dNBxy(1,:);
        Bb(2, 3:3:end) = dNBxy(2,:);
        Bb(3, 2:3:end) = dNBxy(2,:);
        Bb(3, 3:3:end) = dNBxy(1,:);

        Bs = zeros(2, nDofEle);
        Bs(1, 1:3:end) =  dNBxy(1,:);
        Bs(2, 1:3:end) =  dNBxy(2,:);
        Bs(1, 2:3:end) = -NiNodes(1,:,iEval);
        Bs(2, 3:3:end) = -NiNodes(1,:,iEval);

        M(:, iEval, iEle) = Cb * Bb * De;
        Q(:, iEval, iEle) = Cs * Bs * De;
    end
end

end
