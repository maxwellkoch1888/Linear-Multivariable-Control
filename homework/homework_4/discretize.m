% INITIALIZE VARIABLES
initial_time = 0 ;
timestep = 0.1;

% BUILD THE FOUR A MATRICES
A1 = [0, -2; 2, 0];
A2 = [-1, 0; 0, -1];
A3 = [1, -2; 1, 0];
A4 = [-100, 0; 0, 100];
A_matrices = {A1, A2, A3, A4};

% BUILD THE FOUR B MATRICES 
B1 = [0; 1];
B2 = [0; 1];
B3 = [1; 0];
B4 = [0; 0];
B_matrices = {B1, B2, B3, B4};

cases = length(A_matrices);

% BUILD THE SIMULATION 
for k = 1:cases
    A = A_matrices{k};
    B = B_matrices{k};

    % CALCULATE ABAR
    Abar = expm(A*timestep);
    
    % CALCULATE BBAR
    Bbar = computeBbar_ode45(A, B, timestep);

    % PRINT RESULTS
    disp('---------------------')
    disp(['System ', num2str(k)]);
    disp('-----------')
    disp('Abar ='); disp(Abar);
    disp('Bbar ='); disp(Bbar);
    disp('Eigenvalues: '); disp(eig(Abar))
    disp('Eigenvalue mag: '); disp(abs(eig(Abar)))
    disp('---------------------')
    disp('')
end

% CALCUALTE THE VALUE OF BBAR WITH ODE45
function Bbar = computeBbar_ode45(A, B, T)

    n = size(A,1);
    m = size(B,2);

    % MAKE THE INITIAL STATE A SINGLE COLUMN
    x0 = zeros(n*m,1);

    % DEFINE FUNCTION FOR OD45: expm(A*(dt - tau)) * B
    f = @(tau,x) reshape(expm(A*(T - tau))*B, n*m, 1);

    % Integrate from tau=0 to tau=dt
    [~, xvec] = ode45(f, [0 T], x0);

    % Take final value
    Bbar_final = xvec(end,:)';

    % Reshape into n x m
    Bbar = reshape(Bbar_final, n, m);
end