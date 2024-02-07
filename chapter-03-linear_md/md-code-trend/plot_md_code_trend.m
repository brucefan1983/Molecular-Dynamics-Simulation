clear;close all;
gromacs=[1300,1590,1820,2220,2540,2680,3000,3180,3660,3770,4630,5290,6450,7000];
lammps=[566,750,1050,1330,1720,2040,2310,2490,2990,3230,3630,4640,4920,5600];
openmm=[18,26,25,23,33,60,83,145,244,245,452,595,782,880];
hoomd_blue=[11,29,35,48,71,74,73,81,112,134,131,143,144,150];
gpumd=[0,0,0,0,0,0,0,6,10,9,9,26,33,70];
year=2010:2023;
figure;
semilogy(year,gromacs,'d-','linewidth',2); hold on;
semilogy(year,lammps,'s-','linewidth',2); hold on;
semilogy(year,openmm,'s-','linewidth',2); hold on;
semilogy(year,hoomd_blue,'h-','linewidth',2)
plot(year,gpumd,'o-','linewidth',2)
xlabel('year');
ylabel('#citations per year')
legend('GROMACS','LAMMPS','OpenMM','HOOMD-blue','GPUMD')