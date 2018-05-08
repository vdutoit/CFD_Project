% NUMERICAL METHODS IN FLUID MECHANICS
% PROJECT - 23/03/2018
% CHAPELLE GREGOIRE & DUTOIT VALENTIN 
close all;
h = 0.5;
L = 2;
H = 1.5*L;
M = L/h;
N = H/h;
% [X,Y] = meshgrid(linspace(0,L,M),linspace(0,H,N));
dt = 0.1;
t_end = 30;
nt = t_end/dt;

T = importdata('temperature.txt',' ');

fig = figure;
axis equal;
hold on;
for i = 0:nt-1
    pcolor(T(i*N+1:i*N+N,1:M));
    shading interp
    drawnow
end