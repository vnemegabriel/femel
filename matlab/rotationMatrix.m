function R = rotationMatrix(beta,nDofNod)

if nDofNod == 2

    R = [cos(beta) sin(beta)
        -sin(beta) cos(beta)];

end