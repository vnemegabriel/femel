function dN = shapefunsDer(locationGaussPoints, eleType)
% shapefunsDer  Shape function derivatives in isoparametric coordinates
%
%   dN = shapefunsDer(locationGaussPoints, eleType)
%
%   Input:
%     locationGaussPoints - (nPoints x 2) matrix of [xi, eta] coordinates
%                           (also accepts a single row vector for one point)
%     eleType             - 'Q4' or 'Q8'
%
%   Output:
%     dN  - (2 x nNodEle x nPoints) shape function derivative array
%             dN(1,:,p) = dN/dxi  at point p
%             dN(2,:,p) = dN/deta at point p
%
%           For a single point (nPoints == 1) the output is squeezed to
%           (2 x nNodEle) for backward compatibility with existing callers.

nPoints = size(locationGaussPoints, 1);

switch eleType
    case 'Q4';  nNodEle = 4;
    case 'Q8';  nNodEle = 8;
    otherwise;  error('Element type ''%s'' not supported.', eleType);
end

dN = zeros(2, nNodEle, nPoints);

for iPt = 1:nPoints

    xi  = locationGaussPoints(iPt, 1);
    eta = locationGaussPoints(iPt, 2);

    switch eleType

        case 'Q4'
            %xi derivatives
            dN1_dxi  = -0.25*(1-eta);
            dN2_dxi  =  0.25*(1-eta);
            dN3_dxi  =  0.25*(1+eta);
            dN4_dxi  = -0.25*(1+eta);

            %eta derivatives
            dN1_deta = -0.25*(1-xi);
            dN2_deta = -0.25*(1+xi);
            dN3_deta =  0.25*(1+xi);
            dN4_deta =  0.25*(1-xi);

            dN_xi   = [dN1_dxi  dN2_dxi  dN3_dxi  dN4_dxi];
            dN_eta  = [dN1_deta dN2_deta dN3_deta dN4_deta];

        case 'Q8'
            %xi derivatives
            dN8_dxi  =  0.5*(-1)*(1 - eta^2);
            dN7_dxi  =  0.5*(-2*xi)*(1 + eta);
            dN6_dxi  =  0.5*(1)*(1 - eta^2);
            dN5_dxi  =  0.5*(-2*xi)*(1 - eta);
            dN4_dxi  =  0.25*(-1)*(1 + eta) - 0.5*(dN7_dxi + dN8_dxi);
            dN3_dxi  =  0.25*(1)*(1 + eta)  - 0.5*(dN6_dxi + dN7_dxi);
            dN2_dxi  =  0.25*(1)*(1 - eta)  - 0.5*(dN5_dxi + dN6_dxi);
            dN1_dxi  =  0.25*(-1)*(1 - eta) - 0.5*(dN5_dxi + dN8_dxi);

            %eta derivatives
            dN8_deta =  0.5*(-2*eta)*(1 - xi);
            dN7_deta =  0.5*(1)*(1 - xi^2);
            dN6_deta =  0.5*(-2*eta)*(1 + xi);
            dN5_deta =  0.5*(-1)*(1 - xi^2);
            dN4_deta =  0.25*(1)*(1 - xi)   - 0.5*(dN7_deta + dN8_deta);
            dN3_deta =  0.25*(1)*(1 + xi)   - 0.5*(dN6_deta + dN7_deta);
            dN2_deta =  0.25*(-1)*(1 + xi)  - 0.5*(dN5_deta + dN6_deta);
            dN1_deta =  0.25*(-1)*(1 - xi)  - 0.5*(dN5_deta + dN8_deta);

            dN_xi  = [dN1_dxi  dN2_dxi  dN3_dxi  dN4_dxi  dN5_dxi  dN6_dxi  dN7_dxi  dN8_dxi];
            dN_eta = [dN1_deta dN2_deta dN3_deta dN4_deta dN5_deta dN6_deta dN7_deta dN8_deta];
    end

    dN(:,:,iPt) = [dN_xi; dN_eta];
end

% Backward compatibility: single point -> return 2D matrix (2 x nNodEle)
if nPoints == 1
    dN = dN(:,:,1);
end

end
