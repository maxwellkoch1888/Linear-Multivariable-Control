function section_26_control_design()
    % BUILD CONTROLLERS
    create_controller1();
    create_controller2();
end

%% QUESTIONS
% where does the feedforward term come from when the systems are already
% linearized?
% How do we access the functions in pendulum_energy?

%% Functions 

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
    %%DEFINE VALUES
    g = 9.8 ;
    m = 1/9.8;
    l = 0.25; 
    b = 1;

    A = [0,1; g*sqrt(2)/(2*l), -b/(m*l^2)];
    B = [0; 1/(m*l^2)];
end 


function create_controller1()
    disp('Problem 2.3a')
    disp('------------------------------')

    % GET A AND B
    [A,B] = get_23a();
    
    % BUILD CONTROLLER  
    k = build_controller(A,B);
    disp('k:')
    disp(k)
end

function create_controller2()
    disp('Problem 2.3c')
    disp('------------------------------')

    % GET A AND B
    [A,B] = get_23c;
    
    % BUILD CONTROLLER 
    k = build_controller(A,B);
    disp('k:')
    disp(k)
end 

function k = build_controller(A,B)
    [n,m] = size(B);

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