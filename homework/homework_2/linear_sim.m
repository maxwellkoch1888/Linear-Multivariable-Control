function linear_sim()
    %% Initialize the simulation variables
    t0 = 0; % initial time
    dt = 0.1; % time step
    tf = 10.0; % final time
        
    % Set the control input functions
    g = 9.8;
    m = 1/9.8;
    l = 0.25;
    b = 1;

    % Define the equilibrium states
    xeq = [pi/4; 0];
    ueq = b*0 - m*g*l*sin(pi/4);

    delta_u = @(t, x) 0;  
    
    % Set the starting points
    A = [0 1; sqrt(2)*9.8/(2*0.25) -9.8/(0.25^2)];
    [V, ~] = eig(A);
    delta_x0_1 = V(:,2); % State associated with negative eigenvalue
    delta_x0_2 = V(:,1); % State associated with positive eigenvalue
    delta_x0_3 = [-0.05; 0]; % State associated with theta = pi/4 - 0.05, thetad = 0
    delta_x0_list = [delta_x0_1, delta_x0_2, delta_x0_3];
    
    for k = 1:size(delta_x0_list,2)
        % Pull out the state
        delta_x0 = delta_x0_list(:,k);

        %% Simulate and plot the system using ode
        % Simulate the system
        [tvec, delta_xvec] = matlabOde45(delta_x0, t0, dt, tf, delta_u);
        uvec = getControlVector(tvec, delta_xvec, delta_u);

        % Plot the resulting states
        figure;
        plotResults(tvec, delta_xvec, uvec, xeq, ueq, 'b');
    end
end

function u_vec = getControlVector(tvec, xvec, delta_u)
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
        u_vec(:,k) = delta_u(tvec(k), xvec(:,k));
    end

end

function [tvec, xvec] = matlabOde45(delta_x0, t0, dt, tf, delta_u)
    %MatlabOde45 uses ODE 45 to simulate the state starting at x0 from time
    % t0 to tf
    %
    % Inputs:
    %   delta_x0: nx1 initial state
    %   t0: scalar - initial time
    %   dt: scalar - time increment
    %   tf: scalar - final time
    %   delta_u: function handle that takes time and state as inputs and outputs
    %   the control input
    %
    % Outputs:
    %   tvec: 1xm vector of times associated with the states
    %   xvec: nxm matrix of states where each column is a state at the
    %   associated time in tvec
    
    % Initialize the time
    t = t0:dt:tf;
    
    % Simulate the linear output
    [tvec, xvec] = ode45(@(t,x) f(x,delta_u(t,x)), t, delta_x0);

    % Transpose the outputs to get in the correct form
    tvec = tvec';
    xvec = xvec';    

end

function xdot = f(x, u)
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

    % Define the variables
    g = 9.8;
    m = 1/9.8;
    l = 0.25;
    b = 1;
    
    % Linear system matrices
    A = [0 1; sqrt(2)*g/(2*l) -b/(m*l^2)];
    B = [0; 1/(m*l^2)];
    
    % LTI equation
    xdot = A*x + B*u;
end

function plotResults(tvec, delta_xvec, delta_uvec, xeq, ueq, color)
   
    % Build the actual state
    real_state = delta_xvec + xeq;
    real_input = delta_uvec + ueq;
    
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

