function reactionDiffusion
    % Parameters
    R = 10;       % Radius of the circular domain
    Nr = 50;      % Number of grid points in radial direction
    Ntheta = 100; % Number of grid points in angular direction
    Du = 0.1;     % Diffusion coefficient for u
    Dv = 0.05;    % Diffusion coefficient for v
    tspan = [0, 10]; % Time span for the simulation

    % Create spatial grid
    r = linspace(0, R, Nr);
    theta = linspace(0, 2*pi, Ntheta);
    [R, Theta] = meshgrid(r, theta);
    dr = r(2) - r(1);
    dtheta = theta(2) - theta(1);

    % Convert to Cartesian for plotting
    X = R .* cos(Theta);
    Y = R .* sin(Theta);

    % Initial conditions
    u0 = exp(-((X - 0).^2 + (Y - 0).^2));
    v0 = exp(-((X - R/2).^2 + (Y - 0).^2));

    % Flatten initial conditions
    U0 = [u0(:); v0(:)];

    % Time integration
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);
    [t, U] = ode15s(@(t, U) pde_rhs(t, U, Nr, Ntheta, dr, dtheta, Du, Dv), tspan, U0, options);

    % Find indices for plotting
    t_initial = 1;
    t_middle = round(length(t)/2);
    t_final = length(t);

    % Plot initial time
    figure;
    u = reshape(U(t_initial, 1:Nr*Ntheta), [Ntheta, Nr]);
    v = reshape(U(t_initial, Nr*Ntheta+1:end), [Ntheta, Nr]);
    subplot(1, 2, 1);
    surf(X, Y, u, 'EdgeColor', 'none');
    title(['u at t = ', num2str(t(t_initial))]);
    zlim([0 1]);
    xlabel('x'); ylabel('y'); zlabel('u');
    colorbar;
    subplot(1, 2, 2);
    surf(X, Y, v, 'EdgeColor', 'none');
    title(['v at t = ', num2str(t(t_initial))]);
    zlim([0 1]);
    xlabel('x'); ylabel('y'); zlabel('v');
    colorbar;

    % Plot middle time
    figure;
    u = reshape(U(t_middle, 1:Nr*Ntheta), [Ntheta, Nr]);
    v = reshape(U(t_middle, Nr*Ntheta+1:end), [Ntheta, Nr]);
    subplot(1, 2, 1);
    surf(X, Y, u, 'EdgeColor', 'none');
    title(['u at t = ', num2str(t(t_middle))]);
    zlim([0 1]);
    xlabel('x'); ylabel('y'); zlabel('u');
    colorbar;
    subplot(1, 2, 2);
    surf(X, Y, v, 'EdgeColor', 'none');
    title(['v at t = ', num2str(t(t_middle))]);
    zlim([0 1]);
    xlabel('x'); ylabel('y'); zlabel('v');
    colorbar;

    % Plot final time
    figure;
    u = reshape(U(t_final, 1:Nr*Ntheta), [Ntheta, Nr]);
    v = reshape(U(t_final, Nr*Ntheta+1:end), [Ntheta, Nr]);
    subplot(1, 2, 1);
    surf(X, Y, u, 'EdgeColor', 'none');
    title(['u at t = ', num2str(t(t_final))]);
    zlim([0 1]);
    xlabel('x'); ylabel('y'); zlabel('u');
    colorbar;
    subplot(1, 2, 2);
    surf(X, Y, v, 'EdgeColor', 'none');
    title(['v at t = ', num2str(t(t_final))]);
    zlim([0 1]);
    xlabel('x'); ylabel('y'); zlabel('v');
    colorbar;

    function dUdt = pde_rhs(~, U, Nr, Ntheta, dr, dtheta, Du, Dv)
        % Extract u and v
        u = reshape(U(1:Nr*Ntheta), [Ntheta, Nr]);
        v = reshape(U(Nr*Ntheta+1:end), [Ntheta, Nr]);

        % Initialize Laplacians
        Lap_u = zeros(size(u));
        Lap_v = zeros(size(v));

        % Radial Laplacian
        for j = 2:Nr-1
            Lap_u(:, j) = (u(:, j+1) - 2*u(:, j) + u(:, j-1))/dr^2 + ...
                          (1./r(j)) .* (u(:, j+1) - u(:, j-1))/(2*dr);
            Lap_v(:, j) = (v(:, j+1) - 2*v(:, j) + v(:, j-1))/dr^2 + ...
                          (1./r(j)) .* (v(:, j+1) - v(:, j-1))/(2*dr);
        end

        % Angular Laplacian
        for i = 1:Ntheta
            ip = mod(i, Ntheta) + 1;
            im = mod(i-2, Ntheta) + 1;
            Lap_u(i, :) = Lap_u(i, :) + (u(ip, :) - 2*u(i, :) + u(im, :))/(dtheta^2 * r.^2);
            Lap_v(i, :) = Lap_v(i, :) + (v(ip, :) - 2*v(i, :) + v(im, :))/(dtheta^2 * r.^2);
        end

        % Reaction terms
        f = u .* (1 - u) - v .* (u - 0.5);
        g = v .* (1 - v) - u .* (v - 0.5);

        % PDE system
        du_dt = Du * Lap_u + f;
        dv_dt = Dv * Lap_v + g;

        % Flatten
        dUdt = [du_dt(:); dv_dt(:)];
    end
end
