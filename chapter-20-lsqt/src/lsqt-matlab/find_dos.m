function DOS=find_dos(M,E_max,E_scaled,H_scaled,phi)
C=find_moments(M,H_scaled,phi,phi);
DOS=chebyshev_summation(M,C,E_scaled,E_max);
