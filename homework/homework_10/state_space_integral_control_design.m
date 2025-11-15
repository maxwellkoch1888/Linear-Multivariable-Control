function state_space_integral_control_design()
    % problem_A()
    problem_B()

end 

function problem_A()
    % DEFINE A, B, AND TRUE MATRICES
    A = [0,-2;2,0];
    B = [0;1];

    % DEFINE uff TERM
    uff = 4;

    % BUILD Z DOT
    A = [0, -2, 0; 2, 0, 0; 1, 0, 0];
    B = [0;1;0];
    Atrue = [0, -2.1, 0; 2.01, 0, 0; 1, 0, 0];
    Btrue = [0;0.95;0];

    % DETERMINE CONTROLLABILITY
    gamma = ctrb(A,B);
    % disp('rank of gamma:')
    % disp(rank(gamma))

    % BUILD Q AND R WITH BRYSONS METHOD
    xmax = [1/9, 1/16, 1/100];
    umax = [1/4];
    Q = diag(xmax);
    R = diag(umax);
    
    % FIND GAIN MATRIX
    k = lqr(A, B, Q, R);

    % % FIND STEADY STATE
    % xd = [-2;0;0];
    % xss_true = -inv(Atrue-Btrue*k)*Btrue*uff
    % 
    xss = [-2;0;1.15789473684211];
    xdot = Atrue*xss+Btrue*(uff-k*xss)

% (-0.151637625553134*-2*0.95+3.8 - 4.02) / 0.2
end

function problem_B()
    % DEFINE A, B, AND TRUE MATRICES
    A = [0,-2;2,0];
    B = [0;1];

    % DEFINE uff TERM
    uff = 1/3;

    % BUILD Z DOT
    A = [-1, 0, 0; 0, 1, 0; 0, 1, 0];
    B = [0;1;0];
    Atrue = [-0.95, 0, 0; 0, 1.1, 0; 0, 1, 0];
    Btrue = [0;0.97;0];

    % DETERMINE CONTROLLABILITY
    gamma = ctrb(A,B);
    % disp('rank of gamma:')
    % disp(rank(gamma))

    % BUILD Q AND R WITH BRYSONS METHOD
    xmax = [1/36, 1/1, 1/9];
    umax = 1/16;
    Q = diag(xmax)
    R = diag(umax);
    
    % BUILD CONTROLLABLE DECOMPOSITION
    T = [orth(gamma), null(gamma')];
    Abar = inv(T)*A*T;
    Bbar = inv(T) * B;
    Qbar = T'*Q*T

    Abar11 = Abar(1:2,1:2);
    Bbar1 = Bbar(1:2);
    Qbar11 = Qbar(1:2,1:2);

    % FIND GAIN MATRIX
    kbar1 = lqr(Abar11, Bbar1, Qbar11, R);
    kbar = [kbar1, 0]
    k = kbar*inv(T)



    % % FIND STEADY STATE
    % xd = [-2;0;0];
    % xss_true = -inv(Atrue-Btrue*k)*Btrue*uff
    % 
    xss = [0;3;0];
    xdot = Atrue*xss
    xdot = Btrue * (-1/3)
    xdot = Atrue*xss+Btrue*(uff-k*xss);

% (-0.151637625553134*-2*0.95+3.8 - 4.02) / 0.2
end