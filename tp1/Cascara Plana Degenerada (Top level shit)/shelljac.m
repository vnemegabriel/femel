function jac = shelljac(N,dN,zeta,xyz,t,v3)

tt = [t; t; t];
v3t = (v3.*tt)';
jac = [ dN*(xyz + zeta*v3t/2)
         N*(v3t)/2 ];

end 

