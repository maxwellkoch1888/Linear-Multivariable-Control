function decomposition_problem()
    %% Create system
    [A, B, C] = get_system();
    x_d = [-10/3; 3; -3; 43/12; 0];
    disp('---------- System Matrices ----------')
    disp('A:')    
    disp(A)
    disp('B:')   
    disp(B)
    disp('C:')    
    disp(C)

    %% CALCULATE THE uff TERM
    disp('Calculate uff')
    syms uff_sym [3 1]
    uff = solve(B*uff_sym == -A*x_d, uff_sym);
    disp('uff:')
    disp(uff)

    %% CALCULATE THE CONTROLLABILITY AND OBSERVABILITY
    gamma = ctrb(A,B);
    omega = obsv(A,C);
    % disp(rank(gamma))
    % disp(rank(omega))

    %% Design the controller
    % PERFORM CONTROLLABLE DECOMPOSITION
    disp('---------- Design the Controller ----------')
    Tc = [orth(gamma), null(gamma')]
    Acbar = inv(Tc)*A*Tc
    Bcbar = inv(Tc)*B

    % CALCULATE Q AND R
    Qc = diag([1,2,3,4,5])
    Rc = diag([1,2,3]);
    Qcbar = inv(Tc)*Qc*Tc

    % PULL OUT THE CONTROLLABLE PORTION
    Acbar11 = Acbar(1:4, 1:4)
    Bcbar1 = Bcbar(1:4, :)
    Qcbar11 = Qcbar(1:4, 1:4)

    % CALCULATE THE CONTROLLER GAINS
    kcbar11 = lqr(Acbar11, Bcbar1, Qcbar11, Rc)
    kcbar = [kcbar11, [0;0;0]]
    kc = kcbar*inv(Tc);
    disp('k:')
    disp(kc)

    %% Design the observer
    disp('---------- Design the Observer ----------')
    % PERFORM OBSERVABLE DECOMPOSITION
    To = [orth(omega'), null(omega)]
    Aohat = inv(To)*A*To
    Cohat = C*To

    % CALCULATE THE Q MATRIX
    Qo = diag([1,2,3,4,5])
    Ro = diag([1,2,3])
    Qohat = inv(To)*Qo*To

    % PULL OUT THE OBSERVABLE PORTION
    Aohat11 = Aohat(1:4, 1:4)
    Cohat1 = Cohat(:, 1:4)
    Qohat11 = Qohat(1:4,1:4)
    
    % CALCULATE THE OBSERVABLE GAINS
    kohat11 = lqr(Aohat11', Cohat1', Qohat11, Ro)
    kohat = [kohat11, [0;0;0]]
    L = (kohat * inv(To))';
    disp('L:')
    disp(L)
end

function [A, B, C] = get_system()
    A = ...
        [3.0000         0   -4.0000         0         0;
         0    2.2000    0.8000         0    3.6000;
         0         0   -1.0000         0         0;
    1.0000         0    4.0000    4.0000         0;
         0   -1.4000   -0.6000         0   -3.2000];

    B = ...
        [0    2.0000         0;
        0.8000    0.2000    0.6000;
        1.0000    1.0000         0;
             0   -1.0000         0;
       -0.6000   -0.4000   -0.2000];

    C = ...
     [1     0     0     1     0;
     0     0     1     1     0;
     0     2     0     0     1];
end


