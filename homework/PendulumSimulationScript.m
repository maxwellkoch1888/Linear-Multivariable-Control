function PendulumSimulationScript()
    %% Import the xddot and thetaddot functions
    addpath("homework\Xddot.m");
    addpath("homework\Thetaddot.m")
    
    %% Initialize the simulation variables
    t0 = 0; % initial time
    dt = 0.1; % time step
    tf = 10.0; % final time
        
    % Set the control input function
    u = @(t, x) sin(t);    
    
    % Set the starting point
    x0 = [0;0;0];
    
    %% Simulate and plot the system using ode
    % Simulate the system
    [tvec, xvec] = matlabOde45(x0, t0, dt, tf, u);
    uvec = getControlVector(tvec, xvec, u);
    
    % Plot the resulting states
    figure;
    plotResults(tvec, xvec, uvec, 'b');
    
    %% Simulate and plot the system using Euler (or other method)
    % Simulate the system
    [tvec, xvec] = eulerIntegration(x0, t0, dt, tf, u);
    uvec = getControlVector(tvec, xvec, u);
    
    % Plot the results
    plotResults(tvec, xvec, uvec, 'r:');
    
    
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
    u_vec = zeros(1, len);
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
    [tvec, xvec] = ode45(@(t,x) f(t,x,u(t,x)), t, x0);
    
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
    % Write the euler simulation code here
    for i = 1 : len - 1
        xvec(:, i+1) = xvec(:,i) + dt * f(tvec(i), xvec(:, i), u(tvec(i), xvec(:, i)));
    end
end

function zdot = f(t, z, u)
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
    x       = z(1);
    xdot    = z(2);
    theta   = z(3);
    thetadot= z(4);
    
    % Use acceleration functions for xddot and thetaddot
    xddot    = Xddot(x, xdot, theta, thetadot, u);
    thetaddot= Thetaddot(x, xdot, theta, thetadot, u);
    
    % Construct zdot
    zdot = [ xdot;
             xddot;
             thetadot;
             thetaddot ];
end

function plotResults(tvec, xvec, uvec, color)

    % Plot variables
    fontsize = 18;
    linewidth = 2;
    
    % Plot the resulting states
    subplot(5,1,1); hold on;
    plot(tvec, xvec(1,:), color, 'linewidth', linewidth);
    ylabel('x(t)', 'fontsize', fontsize);
    
    subplot(5,1,2); hold on;
    plot(tvec, xvec(2,:), color, 'linewidth', linewidth);
    ylabel('xdot(t)', 'fontsize', fontsize);
    
    subplot(5,1,3); hold on;
    plot(tvec, xvec(3,:), color, 'linewidth', linewidth);
    ylabel('theta(t)', 'fontsize', fontsize);
    
    subplot(5,1,4); hold on;
    plot(tvec, uvec, color, 'linewidth', linewidth);
    ylabel('theta_dot(t)', 'fontsize', fontsize);
    xlabel('time (s)', 'fontsize', fontsize);

    subplot(5,1,5); hold on;
    plot(tvec, uvec, color, 'linewidth', linewidth);
    ylabel('u(t)', 'fontsize', fontsize);
    xlabel('time (s)', 'fontsize', fontsize);
end

