function [w, gp, npg] = gauss(n)
% gauss  Tensor-product 2D Gauss quadrature points and weights
%
%   [w, gp, npg] = gauss(n)
%
%   Input:
%     n   - (2x1) number of Gauss points per direction [n_xi; n_eta]
%           e.g. repmat(2, 2, 1) gives [2; 2] for 2x2 quadrature
%
%   Outputs:
%     w   - (npg x 1) combined weights  w(p) = w_xi(i) * w_eta(j)
%     gp  - (npg x 2) Gauss point coordinates  [xi, eta]
%     npg - total number of Gauss points  (n_xi * n_eta)
%
%   The ordering follows the outer-xi / inner-eta convention, consistent
%   with how shapefunsder and shapefuns index their third dimension.

n_xi  = n(1);
n_eta = n(2);

[w_xi,  gp_xi]  = gauss1D(n_xi);
[w_eta, gp_eta] = gauss1D(n_eta);

npg = n_xi * n_eta;
w   = zeros(npg, 1);
gp  = zeros(npg, 2);

idx = 0;
for i = 1:n_xi
    for j = 1:n_eta
        idx = idx + 1;
        w(idx)    = w_xi(i) * w_eta(j);
        gp(idx,:) = [gp_xi(i), gp_eta(j)];
    end
end

end
