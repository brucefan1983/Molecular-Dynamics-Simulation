function [g]=find_rdf(type1,type2,type,r,box,pbc,Ng,rc)

% determine some parameters
N=size(r,1);              % number of particles
dr=rc/Ng;                 % bin size
N_type1=sum(type==type1); % number of type1 atoms
N_type2=sum(type==type2); % number of type2 atoms
rho=N_type2/det(box);     % particle density for type 2

% accumulate
g=zeros(Ng,1);
for n1=1:N                                    % sum over the atoms
    if type(n1)~=type1                        % type1 is the center atom
        continue;
    end
    for n2=1:N                                % loop over the atoms again
        if type(n2)~=type2 || n1==n2          % type2 is the other atom
            continue;
        end
        r12=r(n2,:)-r(n1,:);                  % position difference vector
        r12=r12.';                            % column vector
        r12=box\r12;                          % transform to cubic box
        r12=r12-pbc.*round(r12);              % mininum image convention
        r12=box*r12;                          % transform back
        d12=sqrt(sum(r12.*r12));              % distance
        if d12<rc                             % there is a cutoff
            index=ceil(d12/dr);               % bin index
            g(index)=g(index)+1;              % accumulate
        end
   end
end

% normalize
for n=1:Ng
    g(n)=g(n)/N_type1;   % average over the center atoms
    dV=4*pi*(dr*n)^2*dr; % volume of a spherical shell
    g(n)=g(n)/dV;        % now g is the local density
    g(n)=g(n)/rho;       % now g is the RDF
end
