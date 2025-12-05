%% INITIALIZE THE SYSTEM
close all 
% BUILD SYSTEM MATRICES
A = [0,1;-2,-3];
B = [0;1];
n = 2;

% BUILD Q, R, S MATRICES
Q = [1,0;0,0];
R = 1;
S = [1,0;0,1];

% BUILD THE SIMULATION VARIABLES
t0 = 0;
T  = 10;
dt = 0.1;
x0 = [1;2];

% PLOT P
plot_p(A, B, Q, R, S, T, dt, t0)

% PLOT ERROR
plot_p_error(A, B, Q, R, S, T, dt, t0)

% PLOT DYNAMICS
plot_dyn(A, B, Q, R, S, T, dt, t0, x0)

%% SIMULATE P
% CALCULATE THE NUMERIC SOLUTION TO THE RICATTI EQUATION
function [p_num, t_num] = numeric(A, B, Q, R, S, T, dt, t0)
    % BUILD VALUES OF X AND Y
    X     = eye(2);
    Y     = S;
    
    % CALCULATE THE M MATRIX
    M = [A, -B*inv(R)*B'; -Q, -A'];
    
    % EULER INTEGRATION
    P11_values = [];
    P12_values = [];
    P21_values = [];
    P22_values = [];
    t_num = [];
    t = t0;
    while t < T
        % CALCULATE THE P MATRIX
        result = expm(M*(t-T)) * [X;Y];
        X_t = result(1:2,1:2);
        Y_t = result(3:4,1:2);
        P_t = Y_t * inv(X_t);
    
        % SEPARATE AND SAVE THE P_MATRIX VALUES FOR PLOTTING
        P11_values = [P11_values, P_t(1,1)];
        P12_values = [P12_values, P_t(1,2)];
        P21_values = [P21_values, P_t(2,1)];
        P22_values = [P22_values, P_t(2,2)];   
        t_num = [t_num, t];
        t = t + dt;
    end

    % BUILD p_num AND t_num
    p_num = [P11_values; P12_values; P21_values; P22_values];
end 

% CALCULATE THE ODE45 SOLULTION TO THE RICATTI EQN
function [p_ode, t_ode] = ode(A, B, Q, R, S, T, dt, t0)
    % SIMULATE WITH ODE45
    [t_ode,p_mat] = ode45(@(t,p_vec) DRE(t, p_vec, A, B, Q, R), [T:-dt:t0], S(:) );
    
    % TRANSPOSE P AND FLIP T TO MATCH DIMENSIONS
    p_ode = p_mat';
    % p_ode = flip(p_ode);
    % t_ode = flip(t_ode);
end

% DIFFERENTIAL RICCATI EQN FOR ODE
function p_dot_vec = DRE(t, p_vec, A, B, Q, R)
    % EXTRACT DIMENSION
    n = size(A,1);

    % RESHAPE VECTOR TO MATRIX
    P = reshape(p_vec, n,n);

    % CREATE DYNAMICS
    p_dot = -A'*P - P*A - Q + P*B*inv(R)*B'*P;

    % RESHAPE RESULT INTO A VECTOR
    p_dot_vec = reshape(p_dot, n*n,1);
end 

%% SIMULATE DYNAMICS
% USE EULER INTEGRATION TO SIMULATE THE DYNAMICS
function [t_vec, x_mat] = euler_dynamics(A, B, Q, R, S, T, dt, t0, x0)
    % Get numeric Riccati solution backward in time
    [p_mat, t_vec] = numeric(A, B, Q, R, S, T, dt, t0);

    t_vec = flip(t_vec);

    N = length(t_vec);
    n = length(x0);

    x_mat = zeros(n, N);
    x_mat(:,1) = x0;

    for k = 2:N
        x_prev = x_mat(:, k-1);
        P_prev = reshape( p_mat(:,k-1), n, n );
        u_prev = - (1/R) * B' * P_prev * x_prev;
        x_mat(:,k) = x_prev + dt*(A*x_prev + B*u_prev);
    end
end

% USE ODE45 TO SIMULATE DYNAMICS
function [t_sol, x_sol] = ode45_dynamics(A, B, Q, R, S, T, x0)

    n = size(A,1);

    p0 = S(:); 

    solP = ode45(@(t,p) DRE(t,p,A,B,Q,R), [T 0], p0);

    function dx = closed_loop(t, x)
        p_vec = deval(solP, t); 
        P = reshape(p_vec, n, n);
        u = -(1/R) * B' * P * x; 
        dx = A*x + B*u;
    end

    [t_sol, x_sol] = ode45(@closed_loop, [0 T], x0);

end

%% PLOTTING FUNCTIONS
% PLOT P 
function plot_p(A, B, Q, R, S, T, dt, t0)
    figure;
    
    [p_ode, t_ode] = ode(A, B, Q, R, S, T, dt, t0);
    [p_num, t_num] = numeric(A, B, Q, R, S, T, dt, t0);
    
    subplot(4,1,1)
    plot(t_num, p_num(1,:), 'b', LineWidth=1.5)
    hold on 
    plot(t_ode, p_ode(1,:), 'r--', LineWidth=1.0)
    xlabel('Time (s)')
    ylabel("P11")
    legend('Numeric', 'ode45')
    hold off
    
    subplot(4,1,2)
    plot(t_num, p_num(2,:), 'b', LineWidth=1.5)
    hold on 
    plot(t_ode, p_ode(2,:), 'r--', LineWidth=1.0)
    xlabel('Time (s)')
    ylabel("P12")
    legend('Numeric', 'ode45')
    hold off
    
    subplot(4,1,3)
    plot(t_num, p_num(3,:), 'b', LineWidth=1.5)
    hold on 
    plot(t_ode, p_ode(3,:), 'r--', LineWidth=1.0)
    xlabel('Time (s)')
    ylabel("P21")
    legend('Numeric', 'ode45')
    hold off
    
    subplot(4,1,4)
    plot(t_num, p_num(4,:), 'b', LineWidth=1.5)
    hold on 
    plot(t_ode, p_ode(4,:), 'r--', LineWidth=1.0)
    xlabel('Time (s)')
    ylabel("P22")
    legend('Numeric', 'ode45')
    hold off

    sgtitle('P Solution Using Matrix Exponential and ODE45')
end 

% PLOT DYNAMICS
function plot_dyn(A, B, Q, R, S, T, dt, t0, x0)
    [t_ode, x_ode] = ode45_dynamics(A, B, Q, R, S, T, x0);
    [t_eul, x_eul] = euler_dynamics(A, B, Q, R, S, T, dt, t0, x0);

    figure;

    subplot(2,1,1)
    plot(t_eul, x_eul(1,:), 'b', LineWidth= 1.5)
    hold on 
    plot(t_ode, x_ode(:,1), 'r--', 'LineWidth', 1.5)
    xlabel('Time (s)')
    ylabel("x_1")
    legend('Numeric', 'ode45')
    hold off

    subplot(2,1,2)
    plot(t_eul, x_eul(2,:), 'b', LineWidth= 1.5)
    hold on 
    plot(t_ode, x_ode(:,2), 'r--', LineWidth=1.5)
    xlabel('Time (s)')
    ylabel("x_2")
    legend('Numeric', 'ode45')
    hold off

    sgtitle('Dynamics Using Euler Integration and ODE45')
end 
    
% PLOT ERROR
function plot_p_error(A, B, Q, R, S, T, dt, t0)
    % CALCULATE TWO SOLUTIONS
    [p_num, t_num] = numeric(A, B, Q, R, S, T, dt, t0);
    [p_ode, t_ode] = ode(A, B, Q, R, S, T, dt, t0);

    % INTERPOLATE DIFFERENCE
    p_num_interp = zeros(size(p_ode));
    for i = 1:4
        p_num_interp(i,:) = interp1(t_num, p_num(i,:), t_ode, 'linear');
    end

    error_mat = abs(p_num_interp - p_ode);

    % Plot error
    figure;
    subplot(4,1,1)
    plot(t_ode, error_mat(1,:), 'r', 'LineWidth', 1.5)
    xlabel('Time (s)')
    ylabel('P11 Error')

    subplot(4,1,2)
    plot(t_ode, error_mat(2,:), 'r', 'LineWidth', 1.5)
    xlabel('Time (s)')
    ylabel('P12 Error')

    subplot(4,1,3)
    plot(t_ode, error_mat(3,:), 'r', 'LineWidth', 1.5)
    xlabel('Time (s)')
    ylabel('P21 Error')

    subplot(4,1,4)
    plot(t_ode, error_mat(4,:), 'r', 'LineWidth', 1.5)
    xlabel('Time (s)')
    ylabel('P22 Error')

    sgtitle('Absolute Error between Numeric and ode45 Riccati Solutions')
end


% % use deval to find 
% p_sol = ode45(@(t,p_vec) DRE(t, p_vec, A, V,Q, inv(R), S[t_f:-dt, t0], S(:)))
% pvecim1 = deval(p_soln,tvec(k-1));
% pim1 = reshape(pvecim1,n,n);

