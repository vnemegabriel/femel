function [Cb, Cs] = getElementConstitutive(material, t)
% getElementConstitutive  Generalized constitutive matrices for Mindlin plate
%
%   [Cb, Cs] = getElementConstitutive(material, t)
%
%   Inputs:
%     material - struct with fields:
%                  .E   Young's modulus
%                  .nu  Poisson's ratio
%                  .G   Shear modulus  (G = E / (2*(1+nu)))
%     t        - plate thickness
%
%   Outputs:
%     Cb - (3x3) bending constitutive matrix:  Cb = t^3/12 * Db
%     Cs - (2x2) shear constitutive matrix:    Cs = 5/6 * t * G * I2
%
%   Generalized stresses:  sigma_hat = Cb * eps_b_hat
%     [Mx; My; Mxy] = Cb * [kappa_x; kappa_y; kappa_xy]
%
%   Generalized shear:     sigma_s_hat = Cs * eps_s_hat
%     [Qx; Qy] = Cs * [gamma_xz; gamma_yz]

E  = material.E;
nu = material.nu;
G  = material.G;

% Plane-stress bending modulus
Ep = E / (1 - nu^2);

% Bending constitutive (thickness integrated analytically: int z^2 dz = t^3/12)
Cb = t^3/12 * [Ep      nu*Ep  0;
               nu*Ep   Ep     0;
               0       0      G];

% Shear constitutive (shear correction factor kappa = 5/6)
% D_hat_s = t * D_s = 5/6 * t * G * I2
Cs = 5/6 * t * G * eye(2);

end
