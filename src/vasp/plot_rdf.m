clear; close all; 

% read in position data
load xf.txt;r=xf(752*64+1:end,1:3); 
N=64; % number of particles
type=zeros(N,1); % atom types
type1=0; % Si
type2=0; % Si

% parameters from MD simulation
pbc=[1,1,1].'; % boundary conditions
a=11;
b=11;
c=11;
alpha=pi/2;
beta=pi/2;
gamma=pi/2;

% choose the number of bins (number of data points in the figure)
Ng=150;

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
Ns=size(r,1)/N  % number of frames

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
plot(r/0.529,g,'o-','linewidth',2);
xlabel('r (Bhor)','fontsize',15)
ylabel('g(r)','fontsize',15)
axis tight;
title('liquid silicon');
set(gca,'fontsize',15,'xtick',0:1:10,'ytick',0:0.2:2.4);
grid on;

