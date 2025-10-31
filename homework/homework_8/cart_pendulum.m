function cart_pendulum()
    create_controller();
end

function create_controller()
    % GET STATE SPACE MATRICES
    [A,B,C,D] = get_dynamics();

    % CALCULATE STABILITY OF EQ POINT
    disp('Stability of theta=pi:')
    disp(eig(A))

    % DESIGN LQR CONTROLLER TO KEEP PENDULUM UPRIGHT
    Q = diag([10, 1, 100, 1]);
    R = 0.1;
    K = lqr(A,B,Q,R);
    disp('LQR gain K = ')
    disp(K)

    % CLOSED LOOP SYSTEM
    A_cl = A - B*K;

    disp('Eigenvalues of (A - BK):')
    disp(eig(A_cl))

    % SIMULATION
    tspan = [0 10];
    x0 = [0; 0; -0.25; 0]; 
    [t, x] = ode45(@(t,x) (A_cl*x), tspan, x0);

    % COMPUTE CONTROL INPUT
    u = -K*x';

    % PLOTS
    figure;
    subplot(5,1,1)
    plot(t, x(:,1))
    xlabel('Time (s)')
    ylabel('x')

    subplot(5,1,2)
    plot(t, x(:,2))
    xlabel('Time (s)')
    ylabel('x dot')

    subplot(5,1,3)
    plot(t, x(:,3))
    xlabel('Time (s)')
    ylabel('theta')

    subplot(5,1,4)
    plot(t, x(:,4))
    xlabel('Time (s)')
    ylabel('theta dot')

    subplot(5,1,5)
    plot(t, u)
    xlabel('Time (s)')
    ylabel('Control Input (u)')

    
end

function [A,B,C,D] = get_dynamics()
% DEFINE PARAMETERS
M = .5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
l = 0.3;

p = I*(M+m)+M*m*l^2;

% BUILD STATE SPACE MODEL 
    A = [0      1              0           0;
         0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
         0      0              0           1;
         0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
    B = [     0;
         (I+m*l^2)/p;
              0;
            m*l/p];
    C = [1 0 0 0;
         0 0 1 0];
    D = [0;
         0];
end