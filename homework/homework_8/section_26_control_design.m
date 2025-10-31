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
function create_controller1()
    disp('Problem 2.3a')
    disp('------------------------------')

    % GET A AND B
    [A,B] = get_23a();
    
    % BUILD GAINS  
    k = build_controller(A,B);
    disp('k:')
    disp(k)

    % BUILD CONTROLLER
    uff = 0;
    u = uff + k;

end

function create_controller2()
    disp('Problem 2.3c')
    disp('------------------------------')

    g = 9.8 ;
    m = 1/9.8;
    l = 0.25; 
    b = 1;

    % GET A AND B
    [A,B] = get_23c;
    
    % BUILD GAINS 
    k = build_controller(A,B);
    disp('k:')
    disp(k)

    % BUILD CONTROLLER
    uff = -m*g*l*sin(pi/4);
    u = uff + k;
end 

function create_controller3()
    % GET A AND B
    [A,B] = get_28();

    % BUILD GAINS
    k = build_controller(A,B);
    disp('k:')
    disp(k)

    % BUILD CONTROLLER
    uff = 0;
    u = uff + k;
end 

function k = build_controller(A,B)
    [n,~] = size(B);

    % BUILD GAMMA AND CHECK CONTROLLABILITY
    gamma = [B, A*B];
    gamma_rank = rank(gamma);
    if gamma_rank == n 
        disp('System is completely controllable, rank of gamma = n.')
    else 
        disp('System is not completely controllable, rank of gamma/= n.')
    end 

    % BUILD Q AND R
    Q = [1/0.5^2, 0; 0, 1/0.1^2];
    R = [1/0.5^2];

    % CALCULATE CONTROL INPUT
    k = lqr(A,B,Q,R);
    disp('Evaluate System Stability:')
    disp(eig(A-B*k))
end 

function [A,B] = get_23a()
    % DEFINE VALUES
    g = 9.8 ;
    m = 1/9.8;
    l = 0.25; 
    b = 1;

    A = [0,1; g/l, -b/(m*l^2)];
    B = [0; 1/(m*l^2)];
end 

function [A,B] = get_23c()
    % DEFINE VALUES
    g = 9.8 ;
    m = 1/9.8;
    l = 0.25; 
    b = 1;

    A = [0,1; g*sqrt(2)/(2*l), -b/(m*l^2)];
    B = [0; 1/(m*l^2)];
end 

function[A,B] = get_28()
    l = 1;
    m = 1; 
    b = 0.1;
    g = 9.8;

    A = [0,1;g/l, -b/(m*l^2)];
    B = [0;1/(m*l^2)];
end 
