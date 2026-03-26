function K = assembleSystem(elementStiffness, rowIndex, columnIndex)
% assembleSystem  Assemble sparse global stiffness matrix from element contributions
%
%   K = assembleSystem(elementStiffness, rowIndex, columnIndex)
%
%   Inputs:
%     elementStiffness - (nEle x 1) cell array of element stiffness matrices Kele
%     rowIndex         - (nEle x 1) cell array of row    DOF index matrices
%     columnIndex      - (nEle x 1) cell array of column DOF index matrices
%
%   Output:
%     K  - sparse global stiffness matrix
%
%   All three cell arrays are produced by getStiffnessMatrix_Plates.
%   Typical usage in a main script:
%
%     elementStiffness = cell(nEle,1);
%     rowIndex         = cell(nEle,1);
%     columnIndex      = cell(nEle,1);
%
%     for iEle = 1:nEle
%         eleDofs = elements.dof(iEle,:)';
%         nodesEle = nodes(elements.connectivity(iEle,:),:);
%         [elementStiffness{iEle}, ~, ~, rowIndex{iEle}, columnIndex{iEle}] = ...
%             getStiffnessMatrix_Plates(nodesEle, eleDofs, Cb, Cs, integrationMethod);
%     end
%
%     K = assembleSystem(elementStiffness, rowIndex, columnIndex);

I = vertcat(rowIndex{:});
J = vertcat(columnIndex{:});
S = vertcat(elementStiffness{:});

K = sparse(I(:), J(:), S(:));

end
