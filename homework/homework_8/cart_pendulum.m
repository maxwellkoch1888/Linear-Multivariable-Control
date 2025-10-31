function cart_pendulum()
    create_controller();
end

function create_controller()
    % GET STATE SPACE MATRICES
    [A,B,C,D] = get_dynamics();

    % POLE PLACEMENT CONTROLLER
    desired_poles = [-10+1i, -10-1i, -1+1i, -1-1i];
    K_place = place(A, B, desired_poles);
    disp('Place gain K = ')
    disp(K_place)

    % LQR CONTROLLER
    Q = diag([10, 1, 100, 1]);
    R = 0.1;
    K_lqr = lqr(A,B,Q,R);
    eig_lqr = eig(A - B*K_lqr);
    disp('LQR gain K = ')
    disp(K_lqr)
    disp('LQR poles:')
    disp(eig_lqr)

    % SIMULATION SETTINGS
    tspan = [0 10];
    x0 = [0; 0; -0.25; 0]; 
    
    % LQR RESPONSE
    [t_lqr, x_lqr] = ode45(@(t,x) (A - B*K_lqr)*x, tspan, x0);
    u_lqr = -K_lqr*x_lqr';
    
    % PLACE RESPONSE
    [t_place, x_place] = ode45(@(t,x) (A - B*K_place)*x, tspan, x0);
    u_place = -K_place*x_place';

    % PLOTS
    figure;
    subplot(2,1,1)
    plot(t_lqr, x_lqr(:,3), 'b', 'LineWidth', 1.5)
    hold on
    plot(t_place, x_place(:,3), 'r--', 'LineWidth', 1.5)
    xlabel('Time (s)')
    ylabel('\theta (rad)')
    legend('LQR', 'Pole Placement')
    title('Pendulum Angle Response')
    
    subplot(2,1,2)
    plot(t_lqr, u_lqr(1,:), 'b', 'LineWidth', 1.5)
    hold on
    plot(t_place, u_place(1,:), 'r--', 'LineWidth', 1.5)
    xlabel('Time (s)')
    ylabel('Control Input u (N)')
    legend('LQR', 'Pole Placement')
    title('Control Input Comparison')
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