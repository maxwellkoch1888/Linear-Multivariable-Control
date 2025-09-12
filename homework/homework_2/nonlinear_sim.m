function nonlinear_sim()
    %% Initialize the simulation variables
    t0 = 0; % initial time
    dt = 0.1; % time step
    tf = 10.0; % final time
        
    % Set the control input functions
    % From 2.3d, T(t) = b*theta_dot - mgl*sin(t)
    g = 9.8;
    m = 1/9.8;
    l = 0.25;
    b = 1;
    u = @(t, x) b*0 - m*g*l*sin(pi/4);  
    
    % Set the starting points
    A = [0 1; sqrt(2)*9.8*(2*0.25) -9.8/(0.25^2)];
    [V, ~] = eig(A);
    x0_1 = V(:,2); % State associated with negative eigenvalue
    x0_2 = V(:,1); % State associated with positive eigenvalue
    x0_3 = [pi/4-0.05; 0]; % State associated with theta = pi/4 - 0.05, thetad = 0
    x0_list = [x0_1, x0_2, x0_3];
    
    for k = 1:size(x0_list,2)
        % Pull out the state
        x0 = x0_list(:,k);

        %% Simulate and plot the system using ode
        % Simulate the system
        [tvec, xvec] = matlabOde45(x0, t0, dt, tf, u);
        uvec = getControlVector(tvec, xvec, u);
        
        % Define the equilibrium states
        xeq = [pi/4; 0];
        ueq = b*0 - m*g*l*sin(pi/4);

        % Plot the resulting states
        figure;
        plotResults(tvec, xvec, uvec, xeq, ueq, 'b');
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
    
    % Simulate the linear output
    [tvec, xvec] = ode45(@(t,x) nonlinear_xdot(x,u(t,x)), t, x0);
    
    % Transpose the outputs to get in the correct form
    tvec = tvec';
    xvec = xvec';    
end

function xdot = nonlinear_xdot(x, u)
    % Calculates the state dynamics using the 

    % Define the variables
    g = 9.8;
    m = 1/9.8;
    l = 0.25;
    b = 1;   
    
    % Pull out states
    theta = x(1);
    theta_dot = x(2);

    % Calculate theta double dot
    theta_ddot = g/l * sin(theta) - b/(m*l^2) * theta_dot + u/(m*l^2);

    % Build and transpose xdot
    xdot = [theta_dot, theta_ddot]';

end

function plotResults(tvec, xvec, uvec, xeq, ueq, color)
   
    % Build the actual state
    real_state = xvec + xeq;
    real_input = uvec + ueq;
    
    % xvec(:,2)
    % real_state(:,2)

    % uvec(:,2)
    % real_input(:,2)
    
    % Plot variables
    fontsize = 18;
    linewidth = 2;
    
    % Plot the resulting states
    subplot(3,1,1); hold on;
    plot(tvec, real_state(1,:), color, 'linewidth', linewidth);
    ylabel('Theta (t)', 'fontsize', fontsize);
    
    subplot(3,1,2); hold on;
    plot(tvec, real_state(2,:), color, 'linewidth', linewidth);
    ylabel('Theta Dot (t)', 'fontsize', fontsize);
    
    subplot(3,1,3); hold on;
    plot(tvec, real_input, color, 'linewidth', linewidth);
    ylabel('u(t)', 'fontsize', fontsize);
    xlabel('Time (s)', 'fontsize', fontsize);
end

