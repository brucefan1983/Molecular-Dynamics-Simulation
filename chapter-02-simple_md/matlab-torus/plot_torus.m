clear; close all;
[theta,phi] = meshgrid(linspace(0,2*pi,50));
r = 2;
R = 5;
x = (R + r*cos(theta)).*cos(phi);
y = (R + r*cos(theta)).*sin(phi);
z = r*sin(theta);
figure;
surf(x,y,z);
shading interp
colormap summer
axis equal
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
set(gca,'Visible','off')
print -dpng test_print
