%np3217 01333401
%AERO96021 High Performance Computing
clear all; close all; clc;
%Script to read in data from vorticity.txt and streamfunction.txt and plot
%the relevant graphs

%import data - each to its own struct
v100 = importdata('Data/vorticity_100.txt');
s100 = importdata('Data/streamfunction_100.txt');
v400 = importdata('Data/vorticity_400.txt');
s400 = importdata('Data/streamfunction_400.txt');
v1000 = importdata('Data/vorticity_1000.txt');
s1000 = importdata('Data/streamfunction_100.txt');
v3200 = importdata('Data/vorticity_3200.txt');
s3200 = importdata('Data/streamfunction_3200.txt');

s2x1 = importdata('Data/streamfunction_2x1.txt');
%s1x2 = importdata('Data/streamfunction_1x2.txt');

%% vorticity and streamfunction contours for Re = 100, Lx = Ly = 1, 1000 timesteps
figure; contourf(v100.data);
title('Voriticty contour for (Lx,Ly) = (1,1), Re = 100');

figure; contourf(s100.data);
title('Streamfunction contour for (Lx,Ly) = (1,1), Re = 100');

%% streamfunction contours for (Lx,Ly) = (2,1) and (Lx,Ly) = (1,2)
figure; contourf(s2x1.data);
title('Streamfunction contour for (Lx,Ly) = (2,1), Re = 100');

%figure; contourf(s1x2.data);
%title('Streamfunction contour for (Lx,Ly) = (1,2), Re - 100');

%% set x and y positions
for i = 1:39
    x(i) = i;
    y(i) = i;
end

dy = 1/(39-1);
dx = 1/(39-1);

%calculate horizontal (u) and vertical (v) velocities
u100 = s100.data/dy;
v100 = -s100.data/dx;
u400 = s400.data/dy;
v400 = -s400.data/dx;
u1000 = s1000.data/dy;
v1000 = -s1000.data/dx;
u3200 = s3200.data/dy;
v3200 = -s3200.data/dx;

%plot of u along x = 0.5 for all Re
figure; hold on; grid on;
plot(u100(:,20),y);
plot(u400(:,20),y);
plot(u1000(:,20),y);
plot(u3200(:,20),y);
legend('Re = 100','Re = 400', 'Re = 1000','Re = 3200','Location','SouthWest');
title('Horizontal velocity along the line x = 0.5');
xlabel('Velocity'); ylabel('y position')

%plot of v along y = 0.5 for all Re
figure; hold on; grid on;
plot(x,v100(20,:));
plot(x,v400(20,:));
plot(x,v1000(20,:));
plot(x,v3200(20,:));
legend('Re = 100','Re = 400','Re = 1000','Re = 3200','Location','NorthWest');
title('Vertical velocity along the line y = 0.5');
xlabel('x position'); ylabel('Velocity')




