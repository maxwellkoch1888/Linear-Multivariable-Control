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
    xeq = [pi; 0];
    ueq = b*0 - m*g*l*sin(pi);

    delta_u = @(t, x) 0;  
    
    % Set the starting points
    delta_x0 = [-0.05; 0]; % State associated with theta = pi - 0.05, thetad = 0

    %% Simulate and plot the system using ode
    % Simulate the system
    [tvec, delta_xvec] = matlabOde45(delta_x0, t0, dt, tf, delta_u);
    uvec = getControlVector(tvec, delta_xvec, delta_u);

    % Calculate the upper bound using lyapunov equation
    % Calculate the value of P
    A = [0, 1; -g/2, -b/(m*l^2)];
    Q = eye(2);
    P = lyap(A', Q);

    Pvals = eig(P);
    Qvals = eig(Q);
    Pmin = min(Pvals);
    mu = -min(Qvals)/max(Pvals);
    V0 = delta_x0' * P * delta_x0;

    convergence_bound = (1 / Pmin) * exp(mu * (tvec - t0)) * V0; 

    % Plot the resulting states
    figure;
    plotResults(tvec, delta_xvec, uvec, xeq, ueq, convergence_bound, 'b');

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
    A = [0, 1; -g/2, -b/(m*l^2)];
    B = [0; 1/(m*l^2)];
    
    % LTI equation
    xdot = A*x + B*u;
end

function plotResults(tvec, delta_xvec, delta_uvec, xeq, ueq, convergence_bound, color)
   
    % Build the actual state
    real_state = delta_xvec + xeq;
    real_input = delta_uvec + ueq;

    % Norm of deviation
    delta_norm = sqrt(sum(delta_xvec.^2,1));
    bound_on_norm = sqrt(convergence_bound);

    % Plotting variables
    fontsize = 18;
    linewidth = 2;
    
    % Plot the resulting states
    subplot(5,1,1); hold on;
    plot(tvec, real_state(1,:), color, 'linewidth', linewidth);
    ylabel('Theta (t)', 'fontsize', fontsize);
    
    subplot(5,1,2); hold on;
    plot(tvec, real_state(2,:), color, 'linewidth', linewidth);
    ylabel('Theta Dot (t)', 'fontsize', fontsize);
    
    subplot(5,1,3); hold on;
    plot(tvec, real_input, color, 'linewidth', linewidth);
    ylabel('u(t)', 'fontsize', fontsize);
    xlabel('Time (s)', 'fontsize', fontsize);

    % Plot convergence bound vs actual state norm
    subplot(5,1,4); hold on;
    plot(tvec, bound_on_norm, 'b', 'linewidth', linewidth); % Lyapunov bound
    ylabel('Lyap Bound', 'fontsize', fontsize);
    xlabel('Time (s)', 'fontsize', fontsize);

    subplot(5,1,5); hold on;
    plot(tvec, delta_norm, 'b', 'linewidth', linewidth);      % actual deviation norm
    ylabel('||\delta x||', 'fontsize', fontsize);
    xlabel('Time (s)', 'fontsize', fontsize);
end
