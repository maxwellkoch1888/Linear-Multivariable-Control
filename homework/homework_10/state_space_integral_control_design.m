function state_space_integral_control_design()
    problem_A()
    % problem_B()

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
    disp('rank of gamma:')
    disp(rank(gamma))

    % BUILD Q AND R WITH BRYSONS METHOD
    xmax = [1/9, 1/16, 1/100];
    umax = [1/4];
    Q = diag(xmax);
    R = diag(umax);
    
    % FIND GAIN MATRIX
    k = lqr(A, B, Q, R);

    % FIND STEADY STATE
    xd = [-2;0;0];
    xdot = A*xd + B*(-k*xd + uff)
end