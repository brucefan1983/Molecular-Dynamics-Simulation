clear; close all;
r0 = [0.0, 0.0, 0.5, 0.5; ...
      0.0, 0.5, 0.0, 0.5; ...
      0.0, 0.5, 0.5, 0.0].';
n0 = size(r0, 1);
nxyz = 7 * [1, 1, 1];
N = nxyz(1) * nxyz(2) * nxyz(3) * n0;
a = 5.385 * [1, 1, 1];
box_length = diag(a .* nxyz);

r = zeros(N, 3);
n = 0;
for nx = 0 : nxyz(1) - 1
    for ny = 0 : nxyz(2) - 1
        for nz = 0 : nxyz(3) - 1
            for m = 1 : n0
                n = n + 1;
                r(n, :) = a .* ([nx,ny,nz] + r0(m, :));   
            end
        end
    end
end

fid = fopen('xyz.in', 'w');
fprintf(fid, '%d\n', N);
fprintf(fid, '%g %g %g %g %g %g %g %g %g\n', box_length);
for n =1 : N
    fprintf(fid, 'Ar %g %g %g %g\n', r(n, :), 40);
end
fclose(fid);
