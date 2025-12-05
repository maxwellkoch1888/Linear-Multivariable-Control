function InvertedPendulum_part_b()
close all;
    % Define parameters
    P.g = 9.8; % Gravity constant
    P.l = 0.25; % Pendulum length
    P.m = 1/9.8; % Pendulum mass
    P.b = 1.0; % Friction coefficient
    
    % Seup for convergences to pi/4ths
    P = controlParamters_PiFourths(P); % Create control
    u = @(t, x) control_PiFourths(t, x, P); % Set the control handle
    x0 = [pi/4 - 0.15; 0; 0]; % Initial state with initial integrator
    
    % Simulate the state forward in time the state
    dt = 0.01;
    t = [0:dt:10];
    [tmat, xmat] = ode45(@(t,x)f_true(t,x,u(t, x), P), t, x0);
    tmat = tmat';
    xmat = xmat';
    
    % Calculate the energy
    len = length(tmat);
    E = zeros(1,len);
    u_mat = zeros(1,len);
    for k = 1:len
        E(k) = calculateEnergy(xmat(:,k), P);
        u_mat(k) = u(tmat(k), xmat(:,k));
    end
    
    %% Plot the results
    fontsize = 12;
    
    % Plot the states
    figure;
    subplot(3,1,1);
    plot([tmat(1) tmat(end)], [P.xd(1) P.xd(1)], ':r', 'linewidth', 3); hold on;
    plot(tmat, xmat(1,:), 'b', 'linewidth', 2);
    ylabel('$\theta(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);
    
    subplot(3,1,2);
    plot([tmat(1) tmat(end)], [P.xd(2) P.xd(2)], 'r:', 'linewidth', 3); hold on;
    plot(tmat, xmat(2,:), 'b', 'linewidth', 3);
    ylabel('$\dot{\theta}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);
    
    % Plot the energy
    subplot(3,1,3);
    plot([tmat(1) tmat(end)], [P.u_ff P.u_ff], 'r:', 'linewidth', 3); hold on;
    plot(tmat, u_mat, 'b', 'linewidth', 2);
    ylabel('Torque');
    xlabel('time (sec)');
    set(gca, 'fontsize', fontsize);
    
    %% Create a simulation
    plotter = PendulumEnergyPlotter(P.l);
    for k = 1:len
        plotter.plot(xmat(1,k), E(k), tmat(k));
        pause(dt/10);
    end

end

function xdot = f_true(t,x, u, P)
    % Define parameters
    g = P.g;
    l = P.l + 0.1;
    m = P.m + 0.2;
    b = P.b - 0.1;
    
    % Extract state
    theta = x(1);
    thetadot = x(2);
    
    % Saturate the input
    u = max(u,-1);
    u = min(u, 1);
    
    % Calculate dynamics
    xdot = zeros(3,1);
    xdot(1) = thetadot;
    xdot(2) = g/l*sin(theta) - b/(m*l^2)*thetadot + 1/(m*l^2)*u;
    xdot(3) = x(1) - P.xd(1); % This should be your integral state
end

function e = calculateEnergy(x, P)
    % Define parameters
    g = P.g;
    l = P.l;
    m = P.m;
    b = P.b;
    
    % Extract state
    theta = x(1);
    thetadot = x(2);
    
    % Calculate position
    x_pos = l*sin(theta);
    y_pos = l*cos(theta);
    
    % Calculate velocity
    x_pos_dot = l*cos(theta)*thetadot;
    y_pos_dot = l*sin(theta)*thetadot;
    v = norm([x_pos_dot; y_pos_dot]);
    
    % Calculate the potential energy
    h = y_pos + l;
    PE = m*g*h;
    
    % Calculate the kinetic energy
    KE = 0.5*m*v^2;
    
    % Return total energy
    e = PE + KE;
end

function P = controlParamters_PiFourths(P)
    % Define the linearized dynamics
    A = [0, 1, 0; ...
        sqrt(2)/2*P.g/P.l, -P.b/(P.m*P.l^2), 0;...
        1, 0, 0];
    B = [0; 1/(P.m*P.l^2); 0];
    Q = diag([1/(0.15^2), 1/(0.5^2), 1/(0.2^2)]);
    R = 1;
    disp('A')
    disp(A)
    disp('B')
    disp(B)
    disp('Q')
    disp(Q)
    disp('R')
    disp(R)
    
    % Design the control    
    P.K = lqr(A,B,Q,R);
    disp('K')
    disp(P.K)
    
    % Define the desired value
    P.xd = [pi/4; 0; 0];
    
    % Create feedforward term
    A_temp = A * P.xd;
    P.u_ff = -A_temp(2)/ B(2);
    % P.u_ff = -0.138840091817449;
    disp('u_ff')
    disp(P.u_ff)
end

function u = control_verticle(t, x, P)
    u = -P.K*x;
    
    % Saturate the input
    u = max(u,-1);
    u = min(u, 1);
end

function u = control_PiFourths(t, x, P)
    % Create difference in theta
    dtheta = x(1) - pi/4;
    dtheta = atan2(sin(dtheta), cos(dtheta)); % adjust difference to be between -pi and pi
    
    % Calcuate control
    dx = [dtheta; x(2); x(3)];
    u = P.u_ff - P.K*dx;
    
    % Saturate the input
    u = max(u,-1);
    u = min(u, 1);
end
