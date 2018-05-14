%{
The original BVE Model was built by Ann Clay, Alexander Donato, and Patrick Sult, with
a huge amount of help from Dr. Thomas Guinn, for WX 462 Numerical Weather Prediction.
Credit goes to Tyler Geen for building the initial code for "relaxationcyclic.m".
The model has been adapted and modified by Alexander Donato, Samuel Rachelson, and 
James Trafny for MA 360 Mathematical Modelling and Simulation I.
%}

clc
clear
close all

%% Define constants needed for the model
Re = 6.371e06;               % Radius of Earth (m)
Omega = 7.292e-05;           % Angular velocity of Earth (s^-1)
lat = 9;                     % Latitude (degrees N) for determining beta
beta = 2*Omega*cosd(lat)/Re; % beta (df/dy) value at "lat" (1/ms)
tol = 1e-7;                  % Defines tolerance for accuracy for the relaxation method

%% Define grid characteristics
jmax = 256;                    % Number of grid points in x dir (i.e., columns)
imax = 256;                    % Number of grid points in y dir (i.e., rows)

psi = zeros(imax,jmax);        % Sets size of psi grid, uses zeroes as initial guess

Lx = 2000*1e03;                % Physical grid width in m (~ a hemisphere)
Ly = 2000*1e03;                % Physical grid height in m (~ equator to pole)

x = linspace(0,Lx,jmax);       % Create x grid
y = linspace(0,Ly,imax);       % Create y grid

dx = x(2)-x(1);                % Grid space distance in x (dx should equal dy).
dy = y(2)-y(1);                % Grid space distance in y (not used).

[X,Y] = meshgrid(x/1e6,y/1e6); % Create meshgrid in 1000's of km for plotting

%% Define time parameters
dt = 120;                   % Time-step in seconds; can be increased to at least 360
time_hrs = 336;              % Total number of hours of simulation
time_sec = time_hrs*3600;   % Total time of simulation in seconds
nmax = ceil(time_sec/dt);   % Total number of iterations needed

%% Create initial condition - Hermite-like
r = zeros(imax, jmax);

yb = 0.98e6;                                % Begining y-location for Hermite
ye = 0.92e6;                                % End y-loccation Hermite
for ii=1:imax/2
    for jj = 1:jmax
        r(ii,jj) = abs((y(ii)-yb)/(ye-yb)); % Bottom half distance
    end
end

yb = 1.02e6;                                % Begining y-location for Hermite
ye = 1.08e6;                                % End y-loccation Hermite
for ii=imax/2:imax
    for jj = 1:jmax
        r(ii,jj) = abs((y(ii)-yb)/(ye-yb)); % Top half distance
    end
end

zeta_grid = zeros(imax, jmax);              % Sets size of zeta grid
zeta_max = 6e-04;                           % Max value at any point

rng(1); % Keeps seed of random number same each time

for ii = 1:imax
    for jj = 1:jmax
        zeta_o = zeta_max*(1+.06*(.5-rand)); % Random number used to add initial disturbance to model
        
        if y(ii) >= 1.08e6 || y(ii) <= 0.92e6
            zeta_grid(ii,jj) = 0;
        elseif y(ii) >= 1.02e6 || y(ii) <= 0.98e6
            zeta_grid(ii,jj) = zeta_o*(1-(3*r(ii,jj)^2)+(2*r(ii,jj)^3)); % Creates Hermite-like shape
        else
            zeta_grid(ii,jj) = zeta_max*(1+.06*(.5-rand));
        end
    end
end

%% Subtracting zetamean
zetamean = sum(sum(zeta_grid))/(imax*jmax); % Starts model with zero initial vorticity
zeta_grid = zeta_grid - zetamean; % Removes average vorticity from every grid point

%% Plot Initial Iteration
contourf(X,Y,zeta_grid,5)
colormap(jet)
colorbar
axis equal
header = sprintf('Initial Vorticity Field at 0 Hours');
title(header)
xlabel('Thousands of Kilometers')
ylabel('Thousands of Kilometers')

figure
plot(y,zeta_grid(:,jmax/2))
set(gca, 'XDir','reverse')
header = sprintf('Initial Vorticity Field (Cross Section)');
title(header)
xlabel('Thousands of Kilometers')
ylabel('Vorticity')
xlim([0.75e6 1.25e6])
ylim([-2e-4 7e-4])

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
filepath = strcat('Results\\', num2str(time_hrs), 'H_', num2str(imax), 'x', num2str(jmax));
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