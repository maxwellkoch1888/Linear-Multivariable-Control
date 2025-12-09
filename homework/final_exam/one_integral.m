function one_integral()
    close all;
    clear all;

    %% Create system
    [A, B, C, E] = get_system();
        
    %% Control Design
    % SOLVE FOR u_ff AND x_d, UNDERDEFINED SYSTEM
    d_nom = 10;

    % MAKE x3 = 20
    % x_d = [x1;x2;20;x4;x5]
    % uff = [uff1;uff2]
    % UNKNOWN VECTOR vars = [x1;x2;x4;x5;u1;u2]
    
    % Build the reduced system: A*x_d + B*uff = b
    A1 = A(:,[1 2 4 5]);   % columns of x except the fixed x3
    A3 = A(:,3)*20;        % contribution of fixed x3
    
    M = [A1 B];            % matrix multiplying [x1;x2;x4;x5;u1;u2]
    rhs = -E*d_nom - A3;
    
    % Minimum-norm solution
    vars = pinv(M)*rhs;
    
    x1 = vars(1);
    x2 = vars(2);
    x4 = vars(3);
    x5 = vars(4);
    uff1 = vars(5);
    uff2 = vars(6);
    
    x_d = [x1; x2; 20; x4; x5];
    u_ff = [uff1; uff2];

    % MAKE SURE uff AND x_d are correct
    % disp(A*x_d + B*u_ff + E*d_nom) % Should be very close to zero

    % BUILD AUGMENTED SYSTEM
    C1 = [0, 0, 1, 0, 0]; % we only care about state 3
    A_aug = [A, zeros(5,1); C1, 0];
    B_aug = [B; zeros(1,2)];

    % CHECK CONTROLLABILITY
    gamma      = ctrb(A,B);
    rank_gamma = rank(gamma);
    % disp('Rank Gamma:')
    % disp(rank_gamma) % rank = 7, completely controllable

    % BUILD Q AND R MATRICES
    Q = diag([0, 0, 1, 0, 0, 1/(20^2)]);

    Cc = [0,0,1,0,0,1/20];
    omega_int = obsv(A_aug,Cc);
    % rank(omega_int)
    
    R = diag([1/(0.5^2), 1/(0.5^2)]);

    % CALCULATE GAINS 
    K_aug = lqr(A_aug, B_aug, Q, R);
    Kx = K_aug(1:2,1:5);
    Ki = K_aug(1:2,6);
  
    % disp('Closed Loop Eigenvalues:')
    % disp(eig(A-B*Kx))

    %% Create the observer (Create the observer)
    omega = obsv(A,C);
    % disp('Rank omega:')
    % disp(rank(omega)) % rank = 5, completely observable

    Q = diag([10,10,10,10,10]);
    R = 1;

    L = lqr(A', C', Q, R)';

    % disp('Observer poles:')
    % disp(eig(A-L*C))
    
    %% Store the values (Feel free to add any additional values here to pass into the control or dynamics functions)
    P.A = A;
    P.B = B;
    P.C = C;
    P.E = E;
    P.n_states = 5;
    P.n_ctrl = 2;
    P.ctrl_max = 100;    

    % ADDED VALUES
    P.uff = u_ff;
    P.Kx = Kx;
    P.Ki = Ki;
    P.x_d = x_d;
    P.d_nom = 10; % from problem statement, assume d=10
    P.L = L;
    

    %% Simulate the system (Don't change anything in this section except x0_ctrl)
    % Create the initial state
    x0_sys = [20.5; 19.; 19.; 11.; 21];
    x0_obs = [15; 15; 15; 15; 15];
    x0_ctrl = [0]; % Any additional states added for control -> initialize each to zero
    x0 = [x0_sys; x0_obs; x0_ctrl];

    % Simulate the system throughout the entire day (86400 seconds)
    t_int = [0 86400];
    [tvec, x_mat] = ode45(@(t,x)dynamics(t, x, P), t_int, x0);
    x_mat = x_mat';

    % Get the control and disturbance over time
    u_mat = get_all_control(x_mat, P);
    d_vec = get_all_disturbances(tvec);
    

    %% Plot the results
    % Plot the states over time
    figure;
    for k = 1:5
        subplot(5,1,k);
        plot([tvec(1) tvec(end)]/3600, [x_d(k) x_d(k)], 'r:', LineWidth=1.5)        
        hold on;
        plot(tvec/3600, x_mat(k+P.n_states, :), 'g', LineWidth=1.5)
        plot(tvec/3600, x_mat(k,:), 'b', LineWidth=1.5);
        ylabel(["x_", num2str(k)])
    end
    xlabel("Time (hr)")
    sgtitle("States vs Time")    

    % Plot other things vs time
    %% Plot state errors over time
    figure;
    for k = 1:5
        subplot(5,1,k);
        e_k = x_mat(k,:) - x_d(k);   % system state error
        plot(tvec/3600, e_k, 'b');
        ylabel(["e_", num2str(k)]);
        grid on;
    end
    xlabel("Time (hr)");
    sgtitle("State Errors vs Time");
end

function [A, B, C, E] = get_system()
    data = load("sys_1.mat");
    A = data.A;
    B = data.B;
    C = data.C;
    E = data.E;
end

function d = disturbance(t)
    % Disturbance as a function of time
    if t < 3600 || t > 82800
        d = 9.0;
    elseif t < 43200
        delta = (t - 3600)/(43200-3600);
        d = 15*(delta) + 9*(1-delta);
    elseif t < 82800
        delta = (t - 43200)/(82800-43200);
        d = 9*(delta) + 15*(1-delta);
    end
    
end

function xdot = dynamics(t, x, P)
    % Defines the dynamics of the system
    % Args:
    %   t: Current time
    %   x: Current augmented state (i.e., actual state, observer state, and
    %   states needed for control)
    %   P: Parameters that are required for simulation
    %
    % Returns:
    %   xdot: The aggregate dynamics

    % Extract the states (don't touch the next three lines)
    x_sys = x(1:P.n_states, :);
    x_obs = x(P.n_states+1:2*P.n_states, :);
    x_ctrl = x(2*P.n_states+1:end, :);

    % Calculate the control (don't touch this line)
    u = control(x_obs, x_ctrl, P);

    % System dynamics (don't touch this line)
    x_sys_dot = P.A*x_sys + P.B*u + P.E*disturbance(t);

    % Observer dynamics (this line should use P.d_nom instead of disturbance(t))
    x_obs_dot = P.A*x_obs + P.B*u + P.E*P.d_nom + P.L*(P.C*x_sys - P.C*x_obs);

    % Control dynamics (definitely fix this line)
    e = x_obs - P.x_d;     % state error

    x_ctrl_dot = [e(3)];     % integrate the error

    xdot = [x_sys_dot; x_obs_dot; x_ctrl_dot];
end

function u = control(x_obs, x_ctrl, P)
    % Defines the control 
    % 
    % Args:
    %   x_obs: The observer state
    %   x_ctrl: The controller state
    %   P: Parameters that are required for control
    %
    % Returns
    %   u: The resulting control

    % Create the feedback control (you'll want to change this)
    e_x = x_obs - P.x_d;   % 5 states
    sigma = x_ctrl;        % integrator state

    u = -P.Kx * e_x - P.Ki * sigma + P.uff;

    % Bound the control (leave the following two lines alone)
    u = max(u, -P.ctrl_max);
    u = min(u, P.ctrl_max);
end

function u_mat = get_all_control(x_mat, P)
    % Loop through and calculate the control given the states
    %
    % Args:
    %   x_mat: matrix of all states over time
    %   P: Necessary parameters for computing control

    % Extract the observer and control states
    x_mat_obs = x_mat(P.n_states+1:2*P.n_states, :);
    x_mat_ctrl = x_mat(2*P.n_states+1:end, :);
    
    % Calculate the control to be returned
    n_times = size(x_mat, 2);
    u_mat = zeros(P.n_ctrl, n_times);
    for k = 1:n_times
        u_mat(:,k) = control(x_mat_obs(:,k), x_mat_ctrl(:,k), P);
    end    
end

function d_vec = get_all_disturbances(tvec)
    % Loop through and calculates the disturbance given the time values
    n_times = length(tvec);
    d_vec = zeros(1,n_times);
    for k = 1:n_times
        d_vec(k) = disturbance(tvec(k));
    end
end
