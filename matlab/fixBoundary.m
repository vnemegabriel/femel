function [fixedDofs, freeDofs] = fixBoundary(nodes, coordIdx, coordVal, dofs2fix, nDofNod, tol)
% fixBoundary  Identify constrained DOFs on a boundary edge
%
%   [fixedDofs, freeDofs] = fixBoundary(nodes, coordIdx, coordVal, dofs2fix, nDofNod)
%   [fixedDofs, freeDofs] = fixBoundary(nodes, coordIdx, coordVal, dofs2fix, nDofNod, tol)
%
%   Inputs:
%     nodes      - (nNod x 2) global node coordinates
%     coordIdx   - which coordinate defines the boundary: 1 = x, 2 = y
%     coordVal   - value of that coordinate on the boundary edge
%     dofs2fix   - local DOF indices to constrain at boundary nodes
%                  e.g. 1:3 for clamped (all), [1] for simply supported (w only)
%                  (local DOF order: 1=w, 2=phi_x, 3=phi_y)
%     nDofNod    - number of DOFs per node (3 for Mindlin plates)
%     tol        - coordinate tolerance (default: 1e-6)
%
%   Outputs:
%     fixedDofs  - global DOF indices to be constrained
%     freeDofs   - remaining free global DOF indices
%
%   Example – clamp all four edges of a plate (nodes, nNod x 2):
%     L = 1.0;
%     [fd1, ~] = fixBoundary(nodes, 1, -L, 1:3, 3);
%     [fd2, ~] = fixBoundary(nodes, 1,  L, 1:3, 3);
%     [fd3, ~] = fixBoundary(nodes, 2, -L, 1:3, 3);
%     [fd4, ~] = fixBoundary(nodes, 2,  L, 1:3, 3);
%     fixedDofs = unique([fd1; fd2; fd3; fd4]);
%     freeDofs  = setdiff(1 : nDofNod*size(nodes,1), fixedDofs)';

if nargin < 6
    tol = 1e-6;
end

nNod    = size(nodes, 1);
nDofTot = nDofNod * nNod;

% Find nodes on the boundary
boundaryNodes = find(abs(nodes(:, coordIdx) - coordVal) < tol);

% Map local DOF offsets to global DOF indices
%   node i -> global DOFs: nDofNod*(i-1)+1 ... nDofNod*(i-1)+nDofNod
fixedDofs = zeros(length(boundaryNodes) * length(dofs2fix), 1);
idx = 0;
for k = 1:length(boundaryNodes)
    i = boundaryNodes(k);
    globalBase = nDofNod * (i - 1);
    for d = dofs2fix(:)'
        idx = idx + 1;
        fixedDofs(idx) = globalBase + d;
    end
end
fixedDofs = unique(fixedDofs(1:idx));

freeDofs = setdiff((1:nDofTot)', fixedDofs);

end
