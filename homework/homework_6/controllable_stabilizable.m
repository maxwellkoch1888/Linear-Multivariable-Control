 %% PROBLEM 2.3 A, B, C
 g = 0.8;
 m = 1/9.8;
 l = 0.25;
 b = 1.0;

% DEFINE THE MATRICIES
A1 = [0, 1; g/l, -b/(m*l^2)];
A2 = A1;
A3 = [0, 1; g*sqrt(2)/(2*l), -b/(m*l^2)];

B = [0; 1/(m*l^2)];

list_of_A = {A1, A2, A3};
list_of_B = {B, B, B};

% DETERMINE IF THE SYSTEMS ARE CONTROLLABLE
disp("-------------------------------------------------")
disp("-------------------------------------------------")
disp("Problem 2.3 a, b, c")
controllable(list_of_A, list_of_B)

%% PROBLEM 2.4 C
k = 1.0;

A = [0, 1; -k, 0];
B = [0;1];

list_of_A = {A};
list_of_B = {B};

% DETERMINE IF THE SYSTEMS ARE CONTROLLABLE
disp("-------------------------------------------------")
disp("Problem 2.4c")
controllable(list_of_A, list_of_B)

%% PROBLEM 2.6 B
g = 0.8;
m = 1/9.8;
I = 1.0;
b = 2.0;

A = [0,1; g*m/I, -b/I];
B = [0; 1/I];

list_of_A = {A};
list_of_B = {B};

% DETERMINE IF THE SYSTEMS ARE CONTROLLABLE
disp("-------------------------------------------------")
disp("Problem 2.6b")
controllable(list_of_A, list_of_B)

%% PROBLEM 2.7 B AND D
A1 = [0,0,0; 0,0,0; 0,0,0];
B1 = [1,0; 0,0; 0,1];

A2 = [0,1,0; -1,0,0; 0,0,0];
B2 = [1,-1; 0,0; 0,1];

list_of_A = {A1, A2};
list_of_B = {B1, B2};

% DETERMINE IF THE SYSTEMS ARE CONTROLLABLE
disp("-------------------------------------------------")
disp("Problem 2.7 b and d")
controllable(list_of_A, list_of_B)

%% Simple System Design Problems

AS1 = [1,2;3,4];
BS1 = [0;1];

AS2 = [2,1;-4,-3];
BS2 = [1,1;-1,-1];

AS3 = [1,2,3; 4,5,6; 7,8,9];
BS3 = [1,0; 0,1; 1,1];

AS4 = [-1,-2; 6,7];
BS4 = [-0.5;1];

AS5 = [3,0,0; 5,-1,1; -5,0,2];
BS5 = [1;1;-1];

AS6 = [2,1; 0,1];
BS6 = [1;-1];

list_of_A = {AS1, AS2, AS3, AS4, AS5, AS6};
list_of_B = {BS1, BS2, BS3, BS4, BS5, BS6};

% DETERMINE IF THE SYSTEMS ARE CONTROLLABLE
disp("-------------------------------------------------")
disp("Simple System Design")
controllable(list_of_A, list_of_B)

%% SATELLITE PROBLEM
mu = 4.306*10^(-3);
A = [0, 0, 1, 0;
     0, 0, 0, 1; 
     3*mu/200, 0, 0, 2*200*sqrt(mu/(200^3));
     0, 0, -2/200*sqrt(mu/(200^3)), 0];
B = [0,0;
     0,0;
     1,0
     0,1/200];

list_of_A = {A};
list_of_B = {B};

% DETERMINE IF THE SYSTEMS ARE CONTROLLABLE
disp("-------------------------------------------------")
disp("Satellite Problem")
controllable(list_of_A, list_of_B)

% PLACE POLES FOR SATELLITE SYSTEM
desiredPoles = [-1, -1.5, -2, -2.5]; 
K = place(A, B, desiredPoles);
disp("Gain matrix K:")
disp(K)

%% FUNCTIONS

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
    disp("gamma: ")
    disp(gamma)
end 

% BUILD THE CHARACTERISTIC EQUATION FOR THE CONTROLLED SYSTEM
function abar_char(A, B)
    % DETERMINE THE DIMENSIONS OF THE MATRICES
    [n,m] = size(B);

    % BUILD THE K MATRIX FOR UP TO n = m = 4
    syms k1 k2 k3 k4 k5 k6 k7 k8
    k_init = [k1, k2, k3, k4; k5, k6, k7, k8];

    % make k mxn matrix
    k = k_init(1:m, 1:n);
    

    % CALCULATE ABAR
    Abar = A - B*k;

    disp("Abar:")
    disp(Abar)
    % disp("k:")
    % disp(k)
    syms s
    Abar_char = det(s*eye(n) - Abar);

    disp("Symbolic Characteristic eqn:")
    disp(simplify(Abar_char))
end

% DETERMINE IF A SYSTEM IS STABILIZABLE
function stabilizable(A,B)
    % FIND THE EIGENVALUES OF THE SYSTEM
    val = eig(A);
    disp(val)

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

% DETERMINE IF A SET OF PROBLEMS IS CONTROLLABLE
function controllable(list_of_A, list_of_B)
    for i = 1:length(list_of_A)
        % MAKE SURE THERE ARE AN EQUAL NUMBER OF A AND B VALUES
        if length(list_of_A) ~= length(list_of_B)
            error("Different number of A and B matricies in lists.")
        end

        A = list_of_A{i};
        B = list_of_B{i};
        disp("System " + i)
        disp("A")
        disp(A)
        disp("B")
        disp(B)
       
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

        % CALCULATE THE ABAR MATRIX AND CHARACTERISTIC EQN
        abar_char(A, B)
        disp("-------------------------------------------------")
    end
end
