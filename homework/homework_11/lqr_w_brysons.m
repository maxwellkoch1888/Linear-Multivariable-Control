%% DEFINE STATE MATRICES
A = [3,6,4;9,6,10;-7,-7,-9];
B = [-2/3, 1/3; 1/3, -2/3; 1/3, 1/3];

xmax = [1/1, 1/100^2, 1/100^2];
umax = [1/25, 1/100];
xmax_terminal = [1/100, 1/1, 1/4];
Q = diag(xmax);
R = diag(umax);
S = diag(xmax_terminal);

% CALCULATE CONTROLLABILITY
gamma = ctrb(A,B);
disp('Rank of Gamma:')
disp(rank(gamma))

T = [orth(gamma), null(gamma')];

Abar = inv(T) * A * T;
Bbar = inv(T) * B;

% INITIALIZE SIMULATION VARIABLES 
t0 = 0;
tf = 1;
dt = 0.01;

% PLOT P
plot_p(A, B, Q, R, S, T, dt, t0)

% PLOT ERROR
plot_p_error(A, B, Q, R, S, T, dt, t0)

% PLOT DYNAMICS
% plot_dyn(A, B, Q, R, S, T, dt, t0, x0)

%% SIMULATE P
% CALCULATE P NUMERICALLY
function [t_num, p_num] = numeric()
    % BUILD VALUES OF X AND Y
    X     = eye(3);
    Y     = S;
    X_dot = A*X - B*R*B'*Y;
    Y_dot = -Q*X - A'*Y;
    
    % CALCULATE THE M MATRIX
    M = [A, -B*R*B'; -Q, -A'];
    
    % EULER INTEGRATION
    P11_values = [];
    P12_values = [];
    P13_values = [];
    P23_values = []; 
    P22_values = [];
    P33_values = [];
    t_num = [];
    
    t = t0;
    while t < tf
        % CALCULATE THE P MATRIX
        result = expm(M*(t-tf)) * [X;Y];
        X_t = result(1:3,1:3);
        Y_t = result(4:6,1:3);
        P_t = Y_t * inv(X_t);
    
        % SEPARATE AND SAVE THE P_MATRIX VALUES FOR PLOTTING
        P11_values = [P11_values, P_t(1,1)];
        P12_values = [P12_values, P_t(1,2)];
        P13_values = [P13_values, P_t(1,3)];    
        P23_values = [P23_values, P_t(2,3)];
        P22_values = [P22_values, P_t(2,2)];   
        P33_values = [P33_values, P_t(3,3)];   
        
        t_num = [t_num, t];
        t = t + dt;
    end
end 

% SOLVE RICATTI EQN USING ODE45
function [t_ode, p_ode] = ode45(A, B, Q, R, S, T, dt, t0)
    % SIMULATE WITH ODE45
    [t_vec,p_mat] = ode45(@(t,p_vec) DRE(t, p_vec, A, B, Q, R), [T:-dt:t0], S(:) );
    
    % TRANSPOSE P AND FLIP T TO MATCH DIMENSIONS
    p_ode = p_mat';
    p_ode = flip(p_ode);
    t_ode = flip(t_vec);
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

%% PLOT THE SOLUTIONS
% PLOT P 
function plot_p(A, B, Q, R, S, T, dt, t0)
    [t_num, p_num]
    figure;
    
    subplot(6,1,1)
    plot(t_num, p_num, 'b', LineWidth=1.5)
    xlabel('Time (s)')
    ylabel("P11")
    legend('Numeric', 'ode45')
    
    subplot(6,1,2)
    plot(t_num, p_num, 'b', LineWidth=1.5)
    xlabel('Time (s)')
    ylabel("P22")
    legend('Numeric', 'ode45')
    
    subplot(6,1,3)
    plot(t_num, p_num, 'b', LineWidth=1.5)
    xlabel('Time (s)')
    ylabel("P33")
    legend('Numeric', 'ode45')
    
    subplot(6,1,4)
    plot(t_num, p_num, 'b', LineWidth=1.5)
    xlabel('Time (s)')
    ylabel("P12")
    legend('Numeric', 'ode45')
    
    subplot(6,1,5)
    plot(t_num, p_num, 'b', LineWidth=1.5)
    xlabel('Time (s)')
    ylabel("P13")
    legend('Numeric', 'ode45')
    
    subplot(6,1,6)
    plot(t_num, p_num, 'b', LineWidth=1.5)
    xlabel('Time (s)')
    ylabel("P23")
    legend('Numeric', 'ode45')

    subplot(6,1,7)
    plot(t_num, p_num, 'b', LineWidth=1.5)
    xlabel('Time (s)')
    ylabel("P23")
    legend('Numeric', 'ode45')    

    subplot(9,1,8)
    plot(t_num, p_num, 'b', LineWidth=1.5)
    xlabel('Time (s)')
    ylabel("P23")
    legend('Numeric', 'ode45')

    subplot(9,1,9)
    plot(t_num, p_num, 'b', LineWidth=1.5)
    xlabel('Time (s)')
    ylabel("P23")
    legend('Numeric', 'ode45')    
end 

% PLOT DYNAMICS
function plot_dyn(A, B, Q, R, S, T, dt, t0, x0)
    [t_ode, x_ode] = ode45_dynamics(A, B, Q, R, S, T, x0);
    [t_eul, x_eul] = euler_dynamics(A, B, Q, R, S, T, dt, t0, x0);

    figure;

    subplot(3,1,1)
    plot(t_eul, x_eul(1,:), 'b', LineWidth= 1.5)
    hold on 
    plot(t_ode, x_ode(:,1), 'r--', 'LineWidth', 1.5)
    xlabel('Time (s)')
    ylabel("x_1")
    legend('Numeric', 'ode45')
    hold off

    subplot(3,1,2)
    plot(t_eul, x_eul(2,:), 'b', LineWidth= 1.5)
    hold on 
    plot(t_ode, x_ode(:,2), 'r--', LineWidth=1.5)
    xlabel('Time (s)')
    ylabel("x_2")
    legend('Numeric', 'ode45')
    hold off

    subplot(3,1,3)
    plot(t_eul, x_eul(3,:), 'b', LineWidth= 1.5)
    hold on 
    plot(t_ode, x_ode(:,3), 'r--', LineWidth=1.5)
    xlabel('Time (s)')
    ylabel("x_3")
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



