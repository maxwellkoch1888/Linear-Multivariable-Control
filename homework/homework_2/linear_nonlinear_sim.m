function linear_nonlinear_sim()
    %% Initialize the simulation variables
    t0 = 0; % initial time
    dt = 0.1; % time step
    tf = 10.0; % final time
        
    % Set the control input functions
    u = @(t, x) 0;    
    
    % Set the starting points
    A = [0 1; sqrt(2)*9.8*(2*0.25) -9.8/(0.25^2)];
    [V, ~] = eig(A);
    x0_1 = V(:,2); % State associated with negative eigenvalue
    x0_2 = V(:,1); % State associated with positive eigenvalue
    x0_3 = V(:,2); % State associated with theta = pi/4 - 0.05, thetad = 0
    x0_list = {x0_1, x0_2, x0_3};
    
    for x0 = x0_list
        %% Simulate and plot the system using ode
        % Simulate the system
        [tvec, xvec] = matlabOde45(x0, t0, dt, tf, u);
        uvec = getControlVector(tvec, xvec, u);
        
        % Plot the resulting states
        figure;
        plotResults(tvec, xvec, uvec, 'b');
    end
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
    
    % Define system parameters
    g = 1.8;
    m = 1/9.8;
    l = 0.25;
    b = 1;

    % Pull out the states
    theta = x(1);
    theta_dot = x(2);

    % Nonlinear system
    theta_ddot = (m*g*l*sin(theta) - b*theta_dot + u) / (m*l^2);

    % Return the state derivative
    xdot = [theta_dot; theta_ddot];
end

function plotResults(tvec, xvec, uvec, color)

    % Plot variables
    fontsize = 18;
    linewidth = 2;
    
    % Plot the resulting states
    subplot(3,1,1); hold on;
    plot(tvec, xvec(1,:), color, 'linewidth', linewidth);
    ylabel('Theta (t)', 'fontsize', fontsize);
    
    subplot(3,1,2); hold on;
    plot(tvec, xvec(2,:), color, 'linewidth', linewidth);
    ylabel('Theta Dot (t)', 'fontsize', fontsize);
    
    subplot(3,1,3); hold on;
    plot(tvec, uvec, color, 'linewidth', linewidth);
    ylabel('u(t)', 'fontsize', fontsize);
    xlabel('Time (s)', 'fontsize', fontsize);
end

