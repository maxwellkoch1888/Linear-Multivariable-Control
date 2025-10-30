function control_decomposition()
    control_design_system_1()
    control_design_system_2()
    control_design_system_3()

end

function [Q,R] = get_Q_R(A,B)
    [n,m] = size(B);
    Q = zeros(n);
    R = zeros(m);
    for i = 1:n
        for j = 1:n
            if i == j 
                Q(i,j) = 1/(i^2);
            end 
        end 
    end
    for i = 1:m
        for j = 1:m
            if i==j
                R(i,j) = 1/(i^2);
            end
        end
    end
end

function control_design_system_1()
    % Get system
    [A, B] = get_system1();
    disp('------------------------------------------------------------------')
    display("Create your controller for system 1");  
    [n, m] = size(B);

    % check controllability
    gamma = ctrb(A,B);
    gamma_rank = rank(gamma);
    disp('rank of gamma:')
    disp(gamma_rank)

    % build transformations
    v = orth(gamma);
    w = null(gamma');
    T = [v, w];
    A_hat = inv(T)*A*T;
    B_hat = inv(T)*B;

    % pull out controllable portion
    A11_hat = A_hat(1:gamma_rank, 1:gamma_rank);
    B1_hat = B_hat(1:gamma_rank, :);

    % check stabilizability
    A22_hat = A_hat(gamma_rank+1:n, gamma_rank+1:n);
    disp('A22 hat eigenvalues:')
    disp(eig(A22_hat))

    % calculate control
    [Q,R] = get_Q_R(A,B);
    Q_hat = T'*Q*T;
    Q11_hat = Q_hat(1:gamma_rank, 1:gamma_rank);

    k1_hat = lqr(A11_hat, B1_hat, Q11_hat, R);
    disp('k1 hat:')
    disp(k1_hat)

    % ensure the system is actually stable 
    k_eig = eig(A11_hat - B1_hat * k1_hat);
    disp('Stabilizable eigenvalues:')
    disp(k_eig)
    disp('All eigenvalues')
    % disp(eig(A_hat-B_hat*k))

    % build actual k and check system
    k_hat = [k1_hat, zeros(m, n - gamma_rank)];
    size(inv(T))
    size(k_hat)
    k = k_hat * inv(T); 
    disp('Final eigenvalues')
    disp(eig(A-B*k))

end

function control_design_system_2()
    % Get system
    [A, B] = get_system2();
    disp('------------------------------------------------------------------')    
    display("Create your controller for system 2");  
    [n, m] = size(B);

    % check controllability
    gamma = ctrb(A,B);
    gamma_rank = rank(gamma);
    disp('rank of gamma:')
    disp(gamma_rank)

    % build transformations
    v = orth(gamma);
    w = null(gamma');
    T = [v, w];
    A_hat = inv(T)*A*T;
    B_hat = inv(T)*B;

    % pull out controllable portion
    A11_hat = A_hat(1:gamma_rank, 1:gamma_rank);
    B1_hat = B_hat(1:gamma_rank, :);

    % check stabilizability
    A22_hat = A_hat(gamma_rank+1:n, gamma_rank+1:n);
    disp(A11_hat)
    disp(B1_hat)    
    disp('A22 hat eigenvalues:')
    disp(eig(A22_hat))
    disp('System is not stabilizable, uncontrollable eigenvalues > 0.')

end

function control_design_system_3()
    % Get system
    [A, B] = get_system3();
    disp('------------------------------------------------------------------')    
    display("Create your controller for system 3");    
    [n, m] = size(B);

    % check controllability
    gamma = ctrb(A,B);
    gamma_rank = rank(gamma);
    disp('rank of gamma:')
    disp(gamma_rank)

    % build transformations
    v = orth(gamma);
    w = null(gamma');
    T = [v, w];
    A_hat = inv(T)*A*T;
    B_hat = inv(T)*B;

    % pull out controllable portion
    A11_hat = A_hat(1:gamma_rank, 1:gamma_rank);
    B1_hat = B_hat(1:gamma_rank, :);

    % check stabilizability
    A22_hat = A_hat(gamma_rank+1:n, gamma_rank+1:n);
    disp('A22 hat eigenvalues:')
    disp(eig(A22_hat))

    % calculate control
    [Q,R] = get_Q_R(A,B);
    Q_hat = T'*Q*T;
    Q11_hat = Q_hat(1:gamma_rank, 1:gamma_rank);

    k1_hat = lqr(A11_hat, B1_hat, Q11_hat, R);
    disp('k1 hat:')
    disp(k1_hat)

    % ensure the system is actually stable 
    k_eig = eig(A11_hat - B1_hat * k1_hat);
    disp('Stabilizable eigenvalues:')
    disp(k_eig)


    % build actual k and check system
    k_hat = [k1_hat, zeros(m, n - gamma_rank)];
    k = k_hat * inv(T); 
    disp('Final eigenvalues')
    disp(eig(A-B*k))
end


function [A, B] = get_system1()
    A = ...
       [-2.200         0   -0.4000         0         0;
             0    1.0000         0   -1.0000    1.0000;
        0.6000         0   -0.8000         0         0;
             0    1.0000         0    3.0000         0;
             0    3.0000         0    2.0000    2.0000];

    B = ...
        [0     0;
         1     0;
         0     0;
        -1     1;
         0     0];
end

function [A, B] = get_system2()

    A = ...
        [2.2000         0    0.4000         0         0;
         0    2.0000         0   -1.0000    3.0000;
   -0.6000         0    0.8000         0         0;
         0    2.0000         0    5.0000   -2.0000;
         0    2.0000         0    2.0000    2.0000];

    B = ...
        [0         0   -0.2000;
         0   -1.0000         0;
         0         0    0.6000;
         0    2.0000         0;
    1.0000    1.0000         0];
end

function [A, B] = get_system3()
    A = ...
        [5.0000         0   -2.0000    2.0000         0;
        0.4000    1.6000    0.8000    0.4000    1.8000;
        2.0000         0    2.0000    2.0000         0;
       -1.0000         0    3.0000    2.0000         0;
       -0.8000   -1.2000   -1.6000   -0.8000   -2.6000];

    B = ...
        [0    2.0000         0;
        0.8000    0.2000    0.6000;
        1.0000    1.0000         0;
             0   -1.0000         0;
       -0.6000   -0.4000   -0.2000];
end


