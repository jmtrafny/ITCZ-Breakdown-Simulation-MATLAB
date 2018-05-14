%{
    Adopted from MA360_BVEModel_DonatoRachaelsonTrafny.m to be at Earth-scale
    instead of zoomed in on a 2000 km patch.  Takes much longer to render.
%}
clc
clear
close all

%% Define constants needed for the model
Re = 6.371e06;               %Radius of Earth (m)
Omega =  7.292e-05;          %Angular velocity of Earth (s^-1)
lat = 9;                     %Latitude (degrees N) for determining beta
beta = 2*Omega*cosd(lat)/Re; % beta (df/dy) value at "lat" (1/ms)
tol = 1e-7;                  % Defines tolerance for accuracy for the relaxation method

%% Define grid characteristics.
Lx = 40000*1e03;         % Circumference of the Earth in meters
Ly = 14000*1e03;         % Pole-to-pole in meters

imax = 256;             % Number of grid points in y dir (i.e., rows)
jmax = floor(imax * (Lx/Ly));  % Number of grid points in x dir (i.e., columns)

psi = zeros(imax,jmax); %Sets size of psi grid, uses zeroes as initial guess

x = linspace(0,Lx,jmax);     % Create x grid
y = linspace(0,Ly,imax);     % Create y grid

dx = x(2)-x(1);        % Grid space distance in x (dx should equal dy).
dy = y(2)-y(1);        % Grid space distance in y (not used).

[X,Y] = meshgrid(x/1e6,y/1e6); % Create meshgrid in 1000's of km for plotting

%% Define time parameters.
dt = 120;                    % Time-step in seconds; can be decreased for finer detail
time_hrs = 256;              % Total number of hours of simulation
time_sec = time_hrs*3600;    % Total time of simulation in seconds
nmax = ceil(time_sec/dt);    % Total number of iterations needed

%% Create initial condition
xc = Lx/2;                  % Used to find center of grid
yc = Ly/2;                  % Used to find center of grid

zeta_max = 6e-04;           % Max Vorticity
a = 1;                      % Gaussian fn peak
b = imax/2;                 % Gaussian fn center
c = 1;                      % Gaussian fn width (1.0 ~~ 400km)

zeta_grid = zeros(imax, jmax); %Sets size of zeta grid

rng(1); % Freezes random seed for testing
for ii = 1:imax
    for jj = 1:jmax
        % creates random initial value between 6.18e-4 -> 5.82e-4
        zeta_o = zeta_max*(1+.06*(.5-rand));  
        % distributes that random value amongst a gaussian curve
        zeta_grid(ii,jj) = zeta_o * a * exp(-((ii-b).^2)/2*(c.^2));
    end 
end

zetamean = sum(sum(zeta_grid))/(imax*jmax); %Starts model with zero initial vorticity
zeta_grid = zeta_grid - zetamean; %Removes average vorticity from every grid point

%% Plot Initial Iteration
contourf(X,Y,zeta_grid,5)
colormap(jet)
colorbar
axis equal
header = sprintf('Initial Vorticity Field at 0 Hours');
title(header)
xlabel('Thousands of Kilometers')
ylabel('Thousands of Kilometers')

%% Model Code
Fn = -Jac_PS(psi,zeta_grid,dx)-dfdx(psi,dx)*beta; % BVE
%    Relative Vort. Adv.  Planetary Vort. Adv.

zetanp1 = zeta_grid + dt * Fn;
psinp1 = relaxationcyclic(psi,zeta_grid,dx,tol); % Uses relaxation code to get psi field

zetanm1 = zeta_grid;
zeta_grid = zetanp1;
psi = psinp1;

iterations = (time_hrs*3600)/ dt - 1;
fprintf('Total Iterations Required >> %d\n', iterations);
cit = 1; % Initiates counter for Courant number array
Cmat = zeros(1,nmax); % Array for Courant numbers to be read into

% create video
filepath = strcat('Results', num2str(time_hrs), 'H_', num2str(imax), 'x', num2str(jmax));
v = VideoWriter(filepath);
v.FrameRate = 30;
open(v);

for ii=2:nmax % Calculates
    Fn = -Jac_PS(psi,zeta_grid,dx)-dfdx(psi,dx)*beta; % BVE
    zetanp1 = zetanm1+2*dt*Fn;
    psinp1 = relaxationcyclic(psi,zetanp1,dx,tol);  
    zetanm1 = zeta_grid+0.25*(zetanm1-2*zeta_grid+zetanp1); % Time filter (provided by Dr. Guinn)    
    zeta_grid = zetanp1;
    psi = psinp1;  
    
    % Print iterations remaining
    iterations = iterations - 1;
    clc;
    fprintf('Iterations remaining: %d\n', iterations);
    
    %% Plotting code
    if mod(ii*dt,0.2*3600) == 0
        figure('visible','off') % Although paper says print 3 hours, this supresses figure and prints at .5 for fluid video
        cit = cit+1; % Counter for Courant array
        ug=-dfdy(psi,dx);
        vg=dfdx(psi,dx);
        c=max(max(sqrt(ug.^2+vg.^2))); % Calculates c
        C = c*(dt/dx); % Calculates Courant Number
        Cmat(cit) = C; % Reads Courant Number of every third iteration into array Cmat
        contourf(X,Y,zeta_grid,5) % Plots all iterations after first
        colormap(jet)
        colorbar
        header = sprintf('Vorticity Field at %d hours, C = %.4f', floor(ii*dt/3600), C);
        title(header)
        xlabel('Thousands of Kilometers')
        ylabel('Thousands of Kilometers')
        axis equal
        
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end

close(v);