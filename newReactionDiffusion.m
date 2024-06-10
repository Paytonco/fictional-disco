clear;clc;

% Define diffusion coefficients for P and Mdm2
Dp = 8172;
Dm = 8172;

% Define a constant s
s = 10;

% Set the maximum time and the time step
tmax = 5/60; % 5 seconds converted to minutes
dt = 0.00001; % time step

% Define the spatial domain for x and y
xspan = -10:0.1:10; % x range from -10 to 10 with step size 0.1
yspan = -10:0.1:10; % y range from -10 to 10 with step size 0.1

% Create a 2D grid of x and y values
[X,Y] = meshgrid(xspan,yspan);

% Initialize time and random initial conditions for P and Mdm2
t = 0;
p = rand(size(X)); % Initial condition for P
m = rand(size(Y)); % Initial condition for Mdm2

% Initialize matrices to store the rate of change of P and Mdm2
dpdtMatrix = zeros(size(p));
dmdtMatrix = zeros(size(m));

% Time-stepping loop
while t < tmax
    % Compute diffusion terms using the Laplace function
    diffusionMatrixP = Dp * laplace(p);
    diffusionMatrixM = Dm * laplace(m);
    
    % Update the rates of change for each point in the grid
    for ix = 1:length(xspan)
        for iy = 1:length(yspan)
            pos = [xspan(ix),yspan(iy)];
            pLoc = p(ix,iy);
            mLoc = m(ix,iy);

            dpdtMatrix(ix,iy) = diffusionMatrixP(ix,iy) + f(pLoc,mLoc,s,pos);
            dmdtMatrix(ix,iy) = diffusionMatrixM(ix,iy) + g(pLoc,mLoc,pos);
        end
    end
    
    % Update P and Mdm2 using Euler's method and ensure non-negative values
    p = max(p + dpdtMatrix * dt,0);
    m = max(m + dmdtMatrix * dt,0);
    t = t + dt;

    % Apply boundary conditions (Dirichlet boundary condition)
    p(1,:) = zeros(size(p(1,:)));
    p(end,:) = zeros(size(p(end,:)));
    p(:,1) = zeros(size(p(:,1)));
    p(:,end) = zeros(size(p(:,end,:)));

    m(1,:) = zeros(size(m(1,:)));
    m(end,:) = zeros(size(m(end,:)));
    m(:,1) = zeros(size(m(:,1)));
    m(:,end) = zeros(size(m(:,end,:)));

    % Plot the current state of P
    figure(1)
    surf(X,Y,p)
    zlabel("P (nM)")
    xlabel("x (\mu m)")
    ylabel("y (\mu m)")
    title("P at t="+3600*t+"s") % Convert t to seconds
    xlim([-10,10])
    ylim([-10,10])

    % Plot the current state of Mdm2
    figure(2)
    surf(X,Y,m)
    title("Mdm2 at t="+3600*t+"s") % Convert t to seconds
    zlabel("M (nM)")
    xlabel("x (\mu m)")
    ylabel("y (\mu m)")
    xlim([-10,10])
    ylim([-10,10])
end

% Function to compute the discrete Laplacian of matrix A
function out = laplace(A)
    kernel = [0,1,0;
              1,-4,1;
              0,1,0];
    out = conv2(A,kernel,"same");
end

% Function to determine if a position is within a specified radius
function out = indicator(pos)
    x = pos(1);
    y = pos(2);
    r = sqrt(x^2 + y^2);

    out = r <= 2.5; % Return 1 if within radius 2.5, else 0
end

% Function f representing the dynamics of P
function out = f(p,m,s,pos)
    k1 = 1000;
    alph = 2;
    kcat = 11;
    km = 10^5;
    q = 0.1;

    I = indicator(pos); % Check if within the specified region

    out = I * k1 * (1 + alph*s^2) / (1 + s^2) - kcat * p * m / (km + p) * (1 - I) - q * p;
end

% Function g representing the dynamics of Mdm2
function out = g(p,m,pos)
    k2 = 0.47;
    alph = 2;
    q = 0.2;

    I = indicator(pos); % Check if within the specified region

    out = I * k2 * (1 + alph * p^2) / (1 + p^2) - q * m;
end

% Uncomment the function below if you need an Euler step function for the ODE
% function u = eulerStep(u0,pos,odeHandle,dt)
%     dxdt = odeHandle(u0,pos);
%     u = u0 + dxdt*dt;
% end
