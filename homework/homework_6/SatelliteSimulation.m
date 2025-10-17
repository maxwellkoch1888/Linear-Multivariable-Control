function SatelliteSimulation()
    close all;
    
    %% Initialize the simulation variables
    % Parameters
    P.mu = 4.302 * 10^(-3); % Graviational constant
    P.R_des = 200; % Desired radius
    
    % Initialize state variables
    r0 = 205;
    rdot0 = 0.0;
    theta0 = 0;
    thetadot0 = 0;
    x0 = [r0; theta0; rdot0; thetadot0];
    
    % Initialize sim variables
    t0 = 0; % initial time
    dt = 0.01; % time step
    tf = 10.0; % final time    
    
    % Calculate the feedback matrix
    K = SatelliteControlDesign(P.R_des);
    % K = zeros(2,4);
    
    % Set the control input function
    u = @(t, x) input(t, x, K, P);   
    
    
    %% Simulate and plot the system using ode
    % Simulate the system
    [tvec, xvec] = matlabOde45(x0, t0, dt, tf, u);
    uvec = getControlVector(tvec, xvec, u);
    
    % Plot the resulting states
    figure;
    plotResults(tvec, xvec, uvec, 'b');
end

function u = input(t, x, K, P)
%input calculates the input for homework 4 satellite problem
%
% Inputs:
%   t: current time
%   x: current state - [r; theta; rdot; thetadot];
%   K: feedback matrix
%   P: struct of useful parameters
%       .mu: Gravitational consant
%       .R_des: Desired radius
%       (anything else you want to have in it)

    % Desired trajectory
    xd = getDesiredState(t, x, P);

    % State deviation
    delta_x = x - xd;

    % Feedback control
    u = -K * delta_x;
end

function xd = getDesiredState(t, x, P)
%getDesiredState calculates the desired state over time
% Inputs:
%   t: current time (not used)
%   x: current state - [r; theta; rdot; thetadot];
%   K: feedback matrix
%   P: struct of useful parameters
%       .mu: Gravitational consant
%       .R_des: Desired radius
%       (anything else you want to have in it)

    w = sqrt(P.mu / P.R_des^3);  % angular velocity for circular orbit
    
    r_des = P.R_des;
    rdot_des = 0;
    theta_des = w * t;   % angle increases linearly with time
    thetadot_des = w;

    xd = [r_des; theta_des; rdot_des; thetadot_des];
end

function u_vec = getControlVector(tvec, xvec, u)
%getControlVector calculate the control vector over the specified time
%interval and state
%
% Inputs:
%   tvec: 1xm vector of time inputs
%   xvec: nxm matrix of states
%   u: function handle that takes time and state as inputs and outputs
%   the control input

    len = size(tvec, 2);
    u_vec = zeros(2, len);
    for k = 1:len
        u_vec(:,k) = u(tvec(k), xvec(:,k));
    end

end

function [tvec, xvec] = matlabOde45(x0, t0, dt, tf, u)
    %MatlabOde45 uses ODE 45 to simulate the state starting at x0 from time
    % t0 to tf
    %
    % Inputs:
    %   x0: nx1 initial state
    %   t0: scalar - initial time
    %   dt: scalar - time increment
    %   tf: scalar - final time
    %   u: function handle that takes time and state as inputs and outputs
    %   the control input
    %
    % Outputs:
    %   tvec: 1xm vector of times associated with the states
    %   xvec: nxm matrix of states where each column is a state at the
    %   associated time in tvec
    
    % Initialize the time
    t = t0:dt:tf;
    
    % Simulate the output
    [tvec xvec] = ode45(@(t,x) f(t,x,u(t,x)), t, x0);
    
    % Transpose the outputs to get in the correct form
    tvec = tvec';
    xvec = xvec';    
end

function [tvec, xvec] = eulerIntegration(x0, t0, dt, tf, u)
    %eulerIntegration uses eulerIntegration to simulate the state starting at x0 from time
    % t0 to tf
    %
    % Inputs:
    %   x0: nx1 initial state
    %   t0: scalar - initial time
    %   dt: scalar - time increment
    %   tf: scalar - final time
    %   u: function handle that takes time and state as inputs and outputs
    %   the control input
    %
    % Outputs:
    %   tvec: 1xm vector of times associated with the states
    %   xvec: nxm matrix of states where each column is a state at the
    %   associated time in tvec
    
    % Initialize values
    tvec = t0:dt:tf;
    len = length(tvec);
    xvec = zeros(size(x0,1), len);
    xvec(:,1) = x0;
    
    % Simulate forward in time
    x = x0;
    for k = 2:len
        % Calculate state dynamics
        t = tvec(k-1);
        xdot = f(t, x, u(t,x));
        
        % Simulate forward in time
        x = x + dt * xdot;
        xvec(:,k) = x;
    end
end

function xdot = f(t, x, u)
    %f calculates the state dynamics using the current time, state, and
    %control input
    %
    % Inputs:
    %   t: current time
    %   x: current state
    %   u: current control input
    %
    % Ouputs:
    %   xdot: time derivative of x(t)
    
    % Dynamic parameters
    mu = 4.302 * 10^(-3);
    
    % Extract states
    r = x(1); % Radius
    %theta = x(2); % Orbit angle
    rdot = x(3); % Time derivative of radius
    thetadot = x(4); % Time derivative of orbit angle
    
    % Extract control inputs
    ar = u(1);
    ai = u(2);
    
    % Calculate the second derivatives
    rddot = thetadot^2*r - mu/(r^2) + ar;
    thetaddot = 1/r*(ai - 2*rdot*thetadot);
    
    % Output the dynamics
    xdot = [rdot; thetadot; rddot; thetaddot];
end

function plotResults(tvec, xvec, uvec, color)

    % Plot variables
    fontsize = 18;
    linewidth = 2;
    
    % Plot the resulting states
    subplot(6,1,1); hold on;
    plot(tvec, xvec(1,:), color, 'linewidth', linewidth);
    ylabel('$r(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(6,1,2); hold on;
    plot(tvec, xvec(2,:), color, 'linewidth', linewidth);
    ylabel('$\theta(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(6,1,3); hold on;
    plot(tvec, xvec(3,:), color, 'linewidth', linewidth);
    ylabel('$\dot{r}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(6,1,4); hold on;
    plot(tvec, xvec(4,:), color, 'linewidth', linewidth);
    ylabel('$\dot{\theta}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(6,1,5); hold on;
    plot(tvec, uvec(1,:), color, 'linewidth', linewidth);
    ylabel('$a_r(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(6,1,6); hold on;
    plot(tvec, uvec(2,:), color, 'linewidth', linewidth);
    ylabel('$a_i(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    xlabel('Time (seconds)', 'fontsize', fontsize);   
    
end

