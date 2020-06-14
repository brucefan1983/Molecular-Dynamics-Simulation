function [g] = find_rdf(r, L, pbc, Ng, rc)

% determine some parameters
N = size(r, 1);         % number of particles
L_times_pbc = L .* pbc; % deal with boundary conditions
rho = N / prod(L);      % global particle density
dr = rc / Ng;           % bin size

% accumulate
g = zeros(Ng, 1);
for n1 = 1 : (N - 1)                               % sum over the atoms
    for n2 = (n1 + 1) : N                          % skipping half of the pairs
       r12 = r(n2, :) - r(n1, :);                  % position difference vector
       r12 = r12 - round(r12 ./L ) .* L_times_pbc; % minimum image convention
       d12 = sqrt(sum(r12 .* r12));                % distance
       if d12 < rc                                 % there is a cutoff
           index = ceil(d12 / dr);                % bin index
           g(index) = g(index) + 1;                % accumulate
       end
   end
end

% normalize
for n = 1 : Ng
    g(n) = g(n) / N * 2;           % 2 because half the pairs have been skipped 
    dV = 4 * pi * (dr * n)^2 * dr; % volume of a spherical shell
    g(n) = g(n) / dV;              % now g is the local density
    g(n) = g(n) / rho;             % now g is the RDF
end
