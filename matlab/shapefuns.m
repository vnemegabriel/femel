function Ni = shapefuns(locationGaussPoints, eleType)
% shapefuns  Shape function values in isoparametric coordinates
%
%   Ni = shapefuns(locationGaussPoints, eleType)
%
%   Input:
%     locationGaussPoints - (nPoints x 2) matrix of [xi, eta] coordinates
%                           (also accepts a single row vector for one point)
%     eleType             - 'Q4' or 'Q8'
%
%   Output:
%     Ni  - (1 x nNodEle x nPoints) shape function array
%             Ni(1,:,p) = [N1, N2, ..., Nn] evaluated at point p
%
%           For a single point (nPoints == 1) the output is (1 x nNodEle)
%           for backward compatibility with existing callers.

nPoints = size(locationGaussPoints, 1);

switch eleType
    case 'Q4';  nNodEle = 4;
    case 'Q8';  nNodEle = 8;
    otherwise;  error('Element type ''%s'' not supported.', eleType);
end

Ni = zeros(1, nNodEle, nPoints);

for iPt = 1:nPoints

    xi  = locationGaussPoints(iPt, 1);
    eta = locationGaussPoints(iPt, 2);

    switch eleType

        case 'Q4'
            N1 = 0.25*(1-xi)*(1-eta);
            N2 = 0.25*(1+xi)*(1-eta);
            N3 = 0.25*(1+xi)*(1+eta);
            N4 = 0.25*(1-xi)*(1+eta);
            Ni(1,:,iPt) = [N1 N2 N3 N4];

        case 'Q8'
            N8 = 0.5*(1-xi)*(1-eta^2);
            N7 = 0.5*(1-xi^2)*(1+eta);
            N6 = 0.5*(1+xi)*(1-eta^2);
            N5 = 0.5*(1-xi^2)*(1-eta);
            N4 = 0.25*(1-xi)*(1+eta) - 0.5*(N7+N8);
            N3 = 0.25*(1+xi)*(1+eta) - 0.5*(N6+N7);
            N2 = 0.25*(1+xi)*(1-eta) - 0.5*(N5+N6);
            N1 = 0.25*(1-xi)*(1-eta) - 0.5*(N8+N5);
            Ni(1,:,iPt) = [N1 N2 N3 N4 N5 N6 N7 N8];
    end
end

% Backward compatibility: single point -> return (1 x nNodEle)
if nPoints == 1
    Ni = Ni(1,:,1);
end

end
