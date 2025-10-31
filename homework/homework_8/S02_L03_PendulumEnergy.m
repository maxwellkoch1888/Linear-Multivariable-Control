function S02_L03_PendulumEnergy()
close all;
    % Define parameters
    P.g = 9.8; % Gravity constant
    P.l = 0.25; % Pendulum length
    P.m = 1/9.8; % Pendulum mass
    P.b = 0.01; % Friction coefficient
    P.x1_max = 0.5; 
    P.x2_max = 0.1; 
    P.u1_max = 0.5;

    % Simulate the state forward in time the state
    x0 = [pi-0.1; 0];
    dt = 0.01;
    t = [0:dt:20];
    [tmat, xmat] = ode45(@(t,x)f(t,x,@zeroControl, P), t, x0);
    tmat = tmat';
    xmat = xmat';
    
    % Calculate the energy
    len = length(tmat);
    E = zeros(1,len);
    umat = zeros(1,len);

    for k = 1:len
        E(k) = calculateEnergy(xmat(:,k), P);
        umat(k) = zeroControl(tmat(k),xmat(:,k),P);
    end
    
    %% Plot the results
    fontsize = 12;
    
    % Plot the states
    subplot(3,1,1);
    plot(tmat, xmat(1,:), 'b', 'linewidth', 3);
    ylabel('$\theta(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);
    
    subplot(3,1,2);
    plot(tmat, xmat(2,:), 'b', 'linewidth', 3);
    ylabel('$\dot{\theta}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);
    
    % Plot the torque
    subplot(3,1,3);
    plot(tmat, umat, 'r', 'linewidth', 3);
    ylabel('Torque');
    xlabel('time (sec)');
    set(gca, 'fontsize', fontsize);
    
    %% Create a simulation
    plotter = PendulumEnergyPlotter(P.l);
    for k = 1:len
        plotter.plot(xmat(1,k), E(k), tmat(k));
        pause(dt/10);
    end

end

function xdot = f(t,x, u_function, P)
    % Define parameters
    g = P.g;
    l = P.l;
    m = P.m;
    b = P.b;
    
    % Calculate the input
    u = u_function(t, x, P);
    
    % Threshold the control input
    T = max(-1, u);
    T = min(1, u);
    
    % Extract state
    theta = x(1);
    thetadot = x(2);
    
    % Calculate dynamics
    xdot = zeros(2,1);
    xdot(1) = thetadot;
    xdot(2) = g/l*sin(theta) - b/(m*l^2)*thetadot + 1/(m*l^2)*T;
end

function u = zeroControl(t, x, P)
    [k,uff] = create_controller2();
    dtheta = x(1) - pi/4;
    dtheta = atan2(sin(dtheta), cos(dtheta));
    u = uff - k*([dtheta; x(2)]);
    
    u = max(u,-1);
    u = min(u,1);
end

function e = calculateEnergy(x, P)
    % Define parameters
    g = P.g;
    l = P.l;
    m = P.m;
    b = P.b;
    
    % Extract state
    theta = x(1);
    thetadot = x(2);
    
    % Calculate position
    x_pos = l*sin(theta);
    y_pos = l*cos(theta);
    
    % Calculate velocity
    x_pos_dot = l*cos(theta)*thetadot;
    y_pos_dot = l*sin(theta)*thetadot;
    v = norm([x_pos_dot; y_pos_dot]);
    
    % Calculate the potential energy
    h = y_pos + l;
    PE = m*g*h;
    
    % Calculate the kinetic energy
    KE = 0.5*m*v^2;
    
    % Return total energy
    e = PE + KE;
end

function section_26_control_design()
    % BUILD CONTROLLERS
    create_controller1();
    create_controller2();
    create_controller3();
end

%% QUESTIONS
% where does the feedforward term come from when the systems are already
% linearized?
% How do we access the functions in pendulum_energy?

%% Functions 
function [k,uff] = create_controller1()
    % disp('Problem 2.3a')
    % disp('------------------------------')

    % GET A AND B
    [A,B] = get_23a();
    
    % BUILD GAINS  
    k = build_controller(A,B);
    % disp('k:')
    % disp(k)

    % BUILD CONTROLLER
    uff = 0;
end

function [k,uff] = create_controller2()
    % disp('Problem 2.3c')
    % disp('------------------------------')

    g = 9.8 ;
    m = 1/9.8;
    l = 0.25; 
    b = 0.01;

    % GET A AND B
    [A,B] = get_23c;
    
    % BUILD GAINS 
    k = build_controller(A,B);
    % disp('k:')
    % disp(k)

    % BUILD CONTROLLER
    uff = -m*g*l*sin(pi/4);
end 

function k = build_controller(A,B)
    [n,~] = size(B);

    % BUILD GAMMA AND CHECK CONTROLLABILITY
    gamma = [B, A*B];
    gamma_rank = rank(gamma);
    % if gamma_rank == n 
    %     disp('System is completely controllable, rank of gamma = n.')
    % else 
    %     disp('System is not completely controllable, rank of gamma/= n.')
    % end 

    % BUILD Q AND R
    Q = [1/0.5^2, 0; 0, 1/0.1^2];
    R = [1/0.5^2];

    % CALCULATE CONTROL INPUT
    k = lqr(A,B,Q,R);

end 

function [A,B] = get_23a()
    % DEFINE VALUES
    g = 9.8 ;
    m = 1/9.8;
    l = 0.25; 
    b = 0.01;

    A = [0,1; g/l, -b/(m*l^2)];
    B = [0; 1/(m*l^2)];
end 

function [A,B] = get_23c()
    % DEFINE VALUES
    g = 9.8 ;
    m = 1/9.8;
    l = 0.25; 
    b = 0.01;

    A = [0,1; g*sqrt(2)/(2*l), -b/(m*l^2)];
    B = [0; 1/(m*l^2)];
end 

