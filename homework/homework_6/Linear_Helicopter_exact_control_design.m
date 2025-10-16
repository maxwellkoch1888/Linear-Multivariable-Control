function Linear_Helicopter_exact_control_design()
close all;

    %% Create the controller
    % Define open loop control parameters
    x1 = zeros(8,1);
    x2 = 0.5.*ones(8,1);
    t1 = 2.5; % Time to get to zero
    t2 = 2.5; % Time to get to x2
    
    % Calculate the control
    P = getSystemMatrices();
        
    % Set the control input function
    u = @(t,x) zeros(4,1);    
    
    %% Initialize the simulation variables
    % Time variables
    t0 = 0; % initial time
    dt = 0.01; % time step
    tf = t1+t2; % final time
    t = t0:dt:tf;
    
    % set initial conditions:
    x0 = x1;
    %x0 = rand(8,1);
    
    %% Simulate and plot the system using ode
    % Simulate the system          
    [tvec xvec] = ode45(@(t,x) f(t,x, u(t,x), P), t, x0);
    tvec = tvec';
    xvec = xvec';  
    
    uvec = getControlVector(tvec, xvec, u);
    
    % Plot the resulting states
    plotResults(tvec, xvec, uvec, x2);
end

function P = getSystemMatrices()
    %Rationalized Helicopter model
    %Code adapted from: http://folk.ntnu.no/skoge/book/2nd_edition/matlab_m/Sec13_2.m
    
    %% Create system matrices 
    % State matrix
    a01 = [          0                  0                  0   0.99857378005981;
                     0                  0   1.00000000000000  -0.00318221934140;
                     0                  0 -11.57049560546880  -2.54463768005371;
                     0                  0   0.43935656547546  -1.99818229675293;
                     0                  0  -2.04089546203613  -0.45899915695190;
    -32.10360717773440                  0  -0.50335502624512   2.29785919189453;
      0.10216116905212  32.05783081054690  -2.34721755981445  -0.50361156463623;
     -1.91097259521484   1.71382904052734  -0.00400543212891  -0.05741119384766];

    a02 = [0.05338427424431             0                  0                  0;
      0.05952465534210                  0                  0                  0;
     -0.06360262632370   0.10678052902222  -0.09491866827011   0.00710757449269;
                     0   0.01665188372135   0.01846204698086  -0.00118747074157;
     -0.73502779006958   0.01925575733185  -0.00459562242031   0.00212036073208;
                     0  -0.02121581137180  -0.02116791903973   0.01581159234047;
      0.83494758605957   0.02122657001019  -0.03787973523140   0.00035400385968;
                     0   0.01398963481188  -0.00090675335377  -0.29051351547241];

    P.A=[a01 a02];

    % Input matrix
    P.B=[              0                  0                  0                  0;
                      0                  0                  0                  0;
       0.12433505058289   0.08278584480286  -2.75247764587402  -0.01788876950741;
      -0.03635892271996   0.47509527206421   0.01429074257612                  0;
       0.30449151992798   0.01495801657438  -0.49651837348938  -0.20674192905426;
       0.28773546218872  -0.54450607299805  -0.01637935638428                  0;
      -0.01907348632812   0.01636743545532  -0.54453611373901   0.23484230041504;
      -4.82063293457031  -0.00038146972656                  0                 0];
  
    % Output matrix
    P.C = [1 0 0 0 0 0 0 0; ... % Pitch
           0 1 0 0 0 0 0 0; ... % Roll
           0 0 0 0 1 0 0 0; ... % Heading velocity (yaw rate)
           0 0 0 0 0 1 0 0]; ... % Heave velocity (forward velocity)
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
    u_vec = zeros(4, len);
    for k = 1:len
        u_vec(:,k) = u(tvec(k), xvec(:,k));
    end

end

function xdot = f(t, x, u, P)
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
    
    % LTI equation for state
    xdot = P.A*x + P.B*u;
end

function plotResults(tvec, xvec, uvec, xd)

    %% Create figure and parameters
    % Figures: inputs and two figures for states
    fig_input = figure('units','normalized','outerposition',[0 0 1 1]);
    fig_state_2 = figure('units','normalized','outerposition',[0 0 1 1]);
    fig_state_1 = figure('units','normalized','outerposition',[0 0 1 1]);
    
    % Plot variables
    fontsize = 18;
    linewidth = 2;
    
    %% Plot desired value
    for k = 1:4
        set(0,'CurrentFigure', fig_state_1);
        subplot(4,1,k); hold on;
        plot([tvec(1) tvec(end)], [xd(k) xd(k)], 'r:', 'linewidth', linewidth);
        
        set(0,'CurrentFigure', fig_state_2);
        subplot(4,1,k); hold on;
        plot([tvec(1) tvec(end)], [xd(k+4) xd(k+4)], 'r:', 'linewidth', linewidth);
    end
    
    %% Label plots
    linewidth = 3;
    set(0,'CurrentFigure', fig_state_1);
    subplot(4,1,1); hold on;
    ylabel('$\theta$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,2); hold on;
    ylabel('$\phi$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,3); hold on;
    ylabel('$p$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,4); hold on;
    ylabel('$q$', 'fontsize', fontsize, 'Interpreter','latex');
    xlabel('Time (seconds)', 'fontsize', fontsize);
    
    set(0,'CurrentFigure', fig_state_2);
    subplot(4,1,1); hold on;
    ylabel('$\xi$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,2); hold on;
    ylabel('$v_x$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,3); hold on;
    ylabel('$v_y$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,4); hold on;
    ylabel('$v_z$', 'fontsize', fontsize, 'Interpreter','latex');
    xlabel('Time (seconds)', 'fontsize', fontsize);
    
    %% Plot the resulting states
    linewidth = 2;
    color = 'b';
    set(0,'CurrentFigure', fig_state_1);
    subplot(4,1,1); hold on;
    plot(tvec, xvec(1,:), color, 'linewidth', linewidth);
    ylabel('$\theta$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,2); hold on;
    plot(tvec, xvec(2,:), color, 'linewidth', linewidth);
    ylabel('$\phi$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,3); hold on;
    plot(tvec, xvec(3,:), color, 'linewidth', linewidth);
    ylabel('$p$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,4); hold on;
    plot(tvec, xvec(4,:), color, 'linewidth', linewidth);
    ylabel('$q$', 'fontsize', fontsize, 'Interpreter','latex');
    xlabel('Time (seconds)', 'fontsize', fontsize);
    
    set(0,'CurrentFigure', fig_state_2);
    subplot(4,1,1); hold on;
    plot(tvec, xvec(5,:), color, 'linewidth', linewidth);
    ylabel('$\xi$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,2); hold on;
    plot(tvec, xvec(6,:), color, 'linewidth', linewidth);
    ylabel('$v_x$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,3); hold on;
    plot(tvec, xvec(7,:), color, 'linewidth', linewidth);
    ylabel('$v_y$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,4); hold on;
    plot(tvec, xvec(8,:), color, 'linewidth', linewidth);
    ylabel('$v_z$', 'fontsize', fontsize, 'Interpreter','latex');
    xlabel('Time (seconds)', 'fontsize', fontsize);
    
    %% Plot inputs
    color = 'b';
    set(0,'CurrentFigure', fig_input);
    subplot(4,1,1); hold on;
    plot(tvec, uvec(1,:), color, 'linewidth', linewidth);
    ylabel('$u_1$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,2); hold on;
    plot(tvec, uvec(2,:), color, 'linewidth', linewidth);
    ylabel('$u_2$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,3); hold on;
    plot(tvec, uvec(3,:), color, 'linewidth', linewidth);
    ylabel('$u_3$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,4); hold on;
    plot(tvec, uvec(4,:), color, 'linewidth', linewidth);
    ylabel('$u_4$', 'fontsize', fontsize, 'Interpreter','latex');
    xlabel('Time (seconds)', 'fontsize', fontsize);
end

