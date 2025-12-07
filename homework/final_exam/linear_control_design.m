% GET SYSTEM MATRICES 
[A, B, C, x_d] = get_system(); 

% DETERMINE THE CONTROLLABILITY OF THE SYSTEM
% disp(size(A)) % 7x7 matrix
gamma = ctrb(A,B); 
% disp(rank(gamma)) % rank = 5, not CC

% DETERMINE THE OBSERVABILITY OF THE SYSTEM 
omega = obsv(A,C); 
% disp(rank(omega)) % rank = 4, not CO

%% CALCULATE A MINIMAL REALIZATION FOR THE STATE SPACE SYSTEM
disp('Minimal Realization')
disp('-------------------------------------------------')

% MINIMALITY WITH MATLAB COMMAND, SECTION 17.5 OF HESPANHA
n = size(C,1);
m = size(B,2);
D = zeros(n,m);

system = ss(A,B,C,D);
[minsys] = minreal(system);
Ahat = minsys.A;
Bhat = minsys.B;
Chat = minsys.C;
Dhat = minsys.D;

% MAKE SURE MINIMAL REPRESENTATION HAS SAME TRANSFER FUNCTION
syms s 
original_sys = simplify(C / (s*eye(7) - A) * B + D);
min_sys = simplify(Chat  / (s*eye(3) - Ahat) * Bhat + Dhat);

% SYMBOLIC TRANSFER FUNCTIONS DON'T LOOK IDENTICAL, CHECK BY NORMALIZING THE
% TRANSFER FUNCTION AT INFINITY
% REFERENCE SECTION 4.2.2-4.2.3 OF LAVRETSKY TO NORMALIZE THE TRANSFER FUNCTION
Omega = logspace(-3,3,200);

% Frequency response of TF difference
G_full = ss(A,B,C,D);
G_min  = ss(Ahat,Bhat,Chat,Dhat);
Gdiff = freqresp(G_full - G_min, Omega);   % size: [ny, nu, N]

% Compute infinity norm: max over frequencies of induced 2-norm
norm_vals = zeros(1, length(Omega));

for k = 1:length(Omega)
    Gk = squeeze(Gdiff(:,:,k));   % convert ny x nu x 1 into ny x nu
    norm_vals(k) = norm(Gk, 2);   % induced 2-norm at this freq
end

norm_diff = max(norm_vals);

disp('Infinity norm of TF difference:');
disp(norm_diff)

% SAVE ANSWERS AS A .mat FILE
save("prob2a.mat", "Ahat", "Bhat", "Chat")

%% CONTROL DESIGN 
disp('Control Design')
disp('-------------------------------------------------')

% z = x - x_d;

% FIND u_ff TERM 
syms u_ff1 u_ff2 u_ff3 u_ff4
u_ff = [u_ff1; u_ff2; u_ff3; u_ff4];

variables = [u_ff1, u_ff2, u_ff3, u_ff4];
eqn = A*x_d + B* u_ff == 0.0;

vars = solve(eqn, variables);
u_ff = [vars.u_ff1; vars.u_ff2; vars.u_ff3; vars.u_ff4];

% WRITE FULL CONTROL
% u = u_hat + u_ff 

% TEST CONTROLLABILITY/STABILIZABILITY
A_list = {A}; 
B_list = {B};
controllable(A_list, B_list); % stabilizable, need to do decomposition

% CONTROLLABLE DECOMPOSITION
T_c = [orth(gamma), null(gamma')];
Ahat_c = inv(T_c) * A * T_c;
Bhat_c = inv(T_c) * B;
Chat_c = C * T_c;

Ahat11_c = Ahat_c(1:5, 1:5);
Bhat1_c = Bhat_c(1:5,1:4);
Chat1_c = Chat_c();

% CALCULATE khat 
% [nhat11,mhat11] = size(Ahat11) % 5x5 matrix
poles = [-1,-2,-3,-4,-5];
Khat11 = place(Ahat11_c, Bhat1_c, poles);
Khat = [Khat11, zeros(4,2)];

% TRANSFORM khat INTO k 
K = Khat * inv(T_c);

% RESULTING CLOSED LOOP EIGENVALUES 
disp('Closed loop eigenvalues:')
disp(eig(A - B * K)) % makes sense, uncontrollable portion has eigenvalues of -2, -2

%% OBSERVER DESIGN
disp('Observer Design')
disp('-------------------------------------------------')

% NOT COMPLETELY OBSERVABLE, PERFORM OBSERVABLE DECOMP
T_o = [orth(omega'), null(omega)];
Ahat_o = inv(T_o) * A * T_o;
Bhat_o = inv(T_o) * B;
Chat_o = C * T_o ;

Ahat11_o = Ahat_o(1:4,1:4);
Ahat22_o = Ahat_o(5:7,5:7);
Chat1_o = Chat_o(1:3,1:4);

disp('Unobservable eigenvalues:')
disp(eig(Ahat22_o)) % eig of unobservable = [1,-1,-2], not detectable. Eig in RHP.
% system fails detectability test, continue observer design as per exam
% instructions.

poles_o = [-1,-2,-3,-4];
Lhat1 = place(Ahat11_o',Chat1_o',poles_o)';
Lhat = [Lhat1; zeros(3,3)];
L = T_o * Lhat;

disp('Resulting eigenvalues of observer error:')
disp(eig(A - L* C)) % resulting eigenvalues 1,-1,-1,-2,-2,-3,-4.
%  Makes sense, unobservable eigenvalues are 1,-1,-2 and commanded
%  eigenvalues are -1,-2,-3,-4. 


%% FUNCTIONS 
% GET SYSTEM MATRICES 
function [A, B, C, x_d] = get_system()
    data = load("sys2.mat");
    A = data.A;
    B = data.B;
    C = data.C;
    x_d = data.x_d;
end

% USE FUNCTIONS FROM MY CH 6 HW TO DETERMINE IF THE SYSTEM IS
% CONTROLLABLE/STABILIZABLE

% DETERMINE IF A SET OF PROBLEMS IS CONTROLLABLE
function controllable(list_of_A, list_of_B)
    for i = 1:length(list_of_A)
        % MAKE SURE THERE ARE AN EQUAL NUMBER OF A AND B VALUES
        if length(list_of_A) ~= length(list_of_B)
            error("Different number of A and B matricies in lists.")
        end

        A = list_of_A{i};
        B = list_of_B{i};
        % disp("System " + i)
        % disp("A")
        % disp(A)
        % disp("B")
        % disp(B)
       
        % CALCULATE GAMMA
        [n, gamma_rank] = matrix_gamma_rank(A, B);
        disp("Matrix n: " + n)
        disp("Gamma Rank: " + gamma_rank)

        % DISPLAY IF THE SYSTEM IN COMPLETELY CONTROLLABLE
        if n == gamma_rank
            disp("System " + i + " is completely controllable.")
        else 
            disp("System " + i + " is NOT completely controllable.")
        end

        % TEST IF THE SYSTEM IS STABILIZABLE
        stabilizable(A,B)
    end
end

% DETERMINE IF A SYSTEM IS STABILIZABLE
function stabilizable(A,B)
    % FIND THE EIGENVALUES OF THE SYSTEM
    val = eig(A);
    % disp(val)

    % FIND THE DIMENSIONS OF THE SYSTEM
    [n,~] = size(A);

    % TEST IF THE SYSTEM IS STABILIZABLE
    j = 1;
    ustabilizable = [];
    for j = 1:length(val)
       lambda = val(j);
       if lambda > 0
           lambda_rank = rank([lambda*eye(n) - A, B]);
           if lambda_rank ~= n
                ustabilizable = [ustabilizable, lambda_rank];
           end
       end
    end
    
    % PRINT STATEMENT
    if (size(ustabilizable) == 0)
        disp("System is stabilizable.")
    else 
        disp("System is NOT stabilizable.")
    end
end

% FIND N AND RANK OF GAMMMA
function [n, gamma_rank] = matrix_gamma_rank(A, B)
    % FIND RANK OF ORIGINAL MATRIX
    [n, ~] = size(A);

    % CALCULATE GAMMA
    gamma = [];
    j = 0;
    while j <= n-1
        gamma = [gamma, A^j * B];
        j = j + 1;
    end 
    
    % DETERMINE THE RANK OF GAMMA
    gamma_rank = rank(gamma);
    % disp("gamma: ")
    % disp(gamma)
end 