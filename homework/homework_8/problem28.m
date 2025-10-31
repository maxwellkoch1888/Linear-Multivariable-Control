function problem28()
    % BUILD CONTROLLERS
    create_controller1();
    create_controller2();
end

%% QUESTIONS
% where does the feedforward term come from when the systems are already
% linearized?
% How do we access the functions in pendulum_energy?

%% Functions 
function create_controller1()
    disp('Problem 2.8a')
    disp('------------------------------')
    g = 9.8 ;
    m = 1/9.8;
    l = 0.25; 
    b = 1;
    x1 = 0;
    x2 = 0; 

    % GET A AND B
    [A,B] = get_part_a();
    
    % BUILD GAINS  
    k = build_controller(A,B);
    disp('k:')
    disp(k)

    % BUILD CONTROLLER
    uff = -g/l*sin(x1);
    disp('uff:')
    disp(uff)

end

function create_controller2()
    disp('Problem 2.8b')
    disp('------------------------------')
    syms x1 x2

    % GET A AND B
    [A,B] = get_part_b;
    
    % BUILD GAINS 
    k = place(A,B,[-1,-2,-3]);
    disp('k:')
    disp(k)
    disp('Evaluate System Stability:')
    disp(eig(A-B*k))


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

function [A,B] = get_part_a()
    % DEFINE VALUES
    g = 9.8 ;
    m = 1/9.8;
    l = 0.25; 
    b = 1;

    A = [0,1; 0, -b/(m*l^2)];
    B = [0; 1];

    % xdot = [x2; g/l*sin(x1) - b/(m*l^2)];
end 

function [A,B] = get_part_b()
    % DEFINE VALUES
    g = 9.8 ;
    m = 1/9.8;
    l = 0.25; 
    b = 1;

    A = [0,1,0;0,0,1;0,0,-b/(m*l^2)];
    B = [0;0;-g/l];
end 
