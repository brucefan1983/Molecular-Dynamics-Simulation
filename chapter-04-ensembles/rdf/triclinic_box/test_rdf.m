clear; close all; 

% read in position data
load r.xyz; 
N=132; % number of particles
% your atom types
% C H N O Ru Zr
% 0 1 2 3 4  5
type=r(1:N,1); % atom types
r=r(:,2:4);  % coordinates in units of Angstrom

% choose the atom types you want to consider:
type1=0; % C
type2=4; % Ru

% parameters from MD simulation
pbc=[1,1,1].'; % boundary conditions
% a=5.6*4;
% b=5.6*4;
% c=5.6*4;
% alpha=pi/2;
% beta=pi/2;
% gamma=pi/2;
a=14.6599;
b=14.6518;
c=14.6567;
alpha=59.9849/180*pi;
beta=59.9923/180*pi;
gamma=59.9893/180*pi;

% choose the number of bins (number of data points in the figure)
Ng=100;

% determine other parameters automatically
ax=a;
bx=b*cos(gamma);
by=b*sin(gamma);
cx=c*cos(beta);
cy=(b*c*cos(alpha)-bx*cx)/by;
cz=sqrt(c*c-cx*cx-cy*cy);
box=[ax,bx cx;0,by,cy;0,0,cz]; % box matrix
volume=det(box);
h1=volume/(a*b*sin(gamma));
h2=volume/(b*c*sin(alpha));
h3=volume/(c*a*sin(beta));
rc=min([h1,h2,h3])/2; % the maximum cutoff radius that can be considered
dr=rc/Ng;         % bin size
Ns=size(r,1)/N;  % number of frames

% do the calculations
g=zeros(Ng,1); % The RDF to be calculated
for n=1:Ns
    r1=r(((n-1)*N+1):(n*N),:); % positions in one frame
    g=g+find_rdf(type1,type2,type,r1,box,pbc,Ng,rc); % sum over frames
end
g=g/Ns; % time average in MD

% plot the data
r=(1:Ng)*dr;
figure;
plot(r,g,'o-','linewidth',2);
xlabel('r (Angstrom)','fontsize',15)
ylabel('g(r)','fontsize',15)
title('C-Ru');
set(gca,'fontsize',15);

