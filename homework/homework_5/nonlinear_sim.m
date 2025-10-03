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
    u = @(t, x) b*0 - m*g*l*sin(pi);  

    % Define the equilibrium states
    xeq = [pi; 0];
    ueq = b*0 - m*g*l*sin(pi);

    % Set the starting points
    x0 = [-0.05; 0] + xeq; % State associated with theta = pi - 0.05, thetad = 0

    %% Simulate and plot the system using ode
    % Simulate the system
    [tvec, xvec] = matlabOde45(x0, t0, dt, tf, u);
    uvec = getControlVector(tvec, xvec, u);


    % Calculate the upper bound using lyapunov equation
    % Calculate the value of P
    A = [0, 1; -g/l, -b/(m*l^2)];
    Q = eye(2);
    P = lyap(A', Q);

    Pvals = eig(P);
    Qvals = eig(Q);
    Pmin = min(Pvals);
    mu = -min(Qvals)/max(Pvals);
    V0 = x0' * P * x0;

    convergence_bound = (1 / Pmin) * exp(mu * (tvec - t0)) * V0; 

    % Plot the resulting states
    figure;
    plotResults(xvec, tvec, uvec, convergence_bound, 'b', xeq);
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

function plotResults(xvec, tvec, uvec, convergence_bound, color, xeq)

    % Norm of deviation
    delta_xvec = xvec - xeq;    
    delta_norm = sqrt(sum(delta_xvec.^2,1));
    bound_on_norm = sqrt(convergence_bound);

    % Plotting variables
    fontsize = 18;
    linewidth = 2;
    
    % Plot the resulting states
    subplot(4,1,1); hold on;
    plot(tvec, xvec(1,:), color, 'linewidth', linewidth);
    ylabel('Theta (t)', 'fontsize', fontsize);
    
    subplot(4,1,2); hold on;
    plot(tvec, xvec(2,:), color, 'linewidth', linewidth);
    ylabel('Theta Dot (t)', 'fontsize', fontsize);
    
    subplot(4,1,3); hold on;
    plot(tvec, uvec, color, 'linewidth', linewidth);
    ylabel('u(t)', 'fontsize', fontsize);
    xlabel('Time (s)', 'fontsize', fontsize);

    % Plot convergence bound vs actual state norm
    subplot(4,1,4); hold on;
    plot(tvec, delta_norm, 'k', 'linewidth', linewidth);             % actual deviation norm
    plot(tvec, bound_on_norm, 'r--', 'linewidth', linewidth);        % Lyapunov bound
    ylabel('||\delta x||', 'fontsize', fontsize);
    xlabel('Time (s)', 'fontsize', fontsize);
    legend('Actual norm','Lyapunov bound','Location','northeast');
end

