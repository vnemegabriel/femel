function D = solveSystem(K, R, freeDofs)
% solveSystem  Solve the reduced plate FEM system K*D = R
%
%   D = solveSystem(K, R, freeDofs)
%
%   Inputs:
%     K         - (nDofTot x nDofTot) global stiffness matrix (sparse or dense)
%     R         - (nDofTot x 1) global load vector
%     freeDofs  - indices of unconstrained DOFs (from fixBoundary)
%
%   Output:
%     D  - (nDofTot x 1) full displacement vector
%          D(fixedDofs) = 0  (homogeneous BCs)
%          D(freeDofs)  = K(free,free) \ R(free)
%
%   DOF order per node: [w, phi_x, phi_y]

R = reshape(R, [], 1);

D = zeros(size(R));
D(freeDofs) = K(freeDofs, freeDofs) \ R(freeDofs);

end
