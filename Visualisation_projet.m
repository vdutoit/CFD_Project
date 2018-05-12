% NUMERICAL METHODS IN FLUID MECHANICS
% PROJECT - 23/03/2018
% CHAPELLE GREGOIRE & DUTOIT VALENTIN 
close all;
M = 128*2;
N = 1.5*M;
H = 1;
L = 2*H/3;
[X,Y] = meshgrid(linspace(0,L,M),linspace(0,H,N));
dt = 0.01;
t_end = 1000;
nt = t_end/dt;
T = importdata('temperature.txt',' ');

fig = figure;

for i = 0:nt-1
    pcolor(X,Y,T(i*N+1:i*N+N,1:M));
    title(['time :',num2str(dt*i,'%.2f')]);
    caxis([-5e-3 5e-3])
    colorbar;
    shading interp;
    axis equal
    drawnow
end
