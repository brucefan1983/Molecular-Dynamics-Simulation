clear; close all; 
load r.txt; % length in units of Angstrom

% parameters from MD simulation
N = 256;                % number of particles
L = 5.60 * [4, 4, 4]; % box size
pbc = [1, 1, 1];        % boundary conditions

% number of bins (number of data points in the figure below)
Ng = 100;

% parameters determined automatically
rc = min(L) / 2;     % the maximum radius
dr = rc / Ng;        % bin size
Ns = size(r, 1) / N; % number of frames

% do the calculations
g = zeros(Ng, 1); % The RDF to be calculated
for n = 1 : Ns
    r1 = r(((n - 1) * N + 1) : (n * N), :); % positions in one frame
    g = g + find_rdf(r1, L, pbc, Ng, rc);    % sum over frames
end
g = g / Ns;                                 % time average in MD

% plot the data
r = (1 : Ng) * dr / 3.405;
figure;
plot(r, g, 'o-');
xlim([0, 3.5]);
ylim([0, 3.5]);
xlabel('r^{\ast}', 'fontsize', 15)
ylabel('g(r)', 'fontsize', 15)
set(gca, 'fontsize', 15);

