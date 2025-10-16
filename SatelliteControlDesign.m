function K = SatelliteControlDesign(R_des)
% Inputs: R_des - desired radius of the orbit
%
% Outputs: K - Feedback control matrix

%% Create the dynamics using symbolic variables
% Symbolic variables for the states
theta = sym('theta', 'real'); % Orbit angle
thetadot = sym('thetadot', 'real'); % Orbit angle time derivative
r = sym('r', 'real'); % Radius of orbit
rdot = sym('rdot', 'real'); % time derivative of radius

% Symbolic variables for the inputs
ar = sym('ar', 'real'); % Radial component of acceleration
ai = sym('ai', 'real'); % In track component of acceleration

% Symbolic constants
mu = sym('mu', 'real'); % Gravity constant

% Create dynamics
rddot = thetadot^2*r - mu/(r^2) + ar; % Second time-derivative of the radius
thetaddot = 1/r*(ai - 2*rdot*thetadot); % Second time-derivative of the orbit angle

%% Create the dynamic equations
% State and input vectors (Define the 4x1 state x and the 2x1 input u as
% vectors of previously defined variables)
x = [r; theta; ]; %... Need to completed this vector
u = []; % Need to complete this vector

% State dynamics (create a vector of the dynamics, i.e., f = [rdot; thetadot; rddot; thetaddot] )


%% Find trajectory to linearize about
R = sym('R', 'real'); % Desired radius of orbit

% Create a trajectory and solve for the equilibrium
% rdot and rddot must be zero (i.e. f(1) and f(3) must be zero) - use the
% subs command
f_sol = subs(f, r, R); % This substitues the big R variable for the little r
f_sol = subs(f_sol, ar, 0); % This substitutes in 0 for ar
f_sol = []; % Substitute 0 in for ai

% Solve for the solution of thetadot (you can use Matlab's "solve" command)
thetadot_sol = solve([],[]); % You can do it using two inputs, the third element of f_sol and the symbolic variable for \dot{\theta}

%% Linear system about the trajectory
% Create the jacobians
df_dx = jacobian(f, x); % Creates a jacobian by evaluating f wrt x
df_du = []; % Create a jacobian with respect to u

% Create the solution trajectory by defining all of the variables in the
% trajectory - replace all "[]" values with the correct value
w = sym('w', 'real'); % Solution for theta_dot
t = sym('t', 'real'); % Time
r = [];
rdot = [];
theta = [];
thetadot = [];
ar = [];
ai = [];

% Evaluate linearization at the trajectory solution
df_dx_sol = simplify(subs(df_dx)); % The "subs" command will substitute all of the redefined variables into df_dx
df_du_sol = []; % Do a very similar thing to get the solution for df_du

%% Evaluate the linearized system
mu = 4.302 * 10^(-3);
R = R_des; % Note, this is passed into the function as a parameter
w =[]; % What is the solution

% Get the A and B matrices (need to convert the symbolic matrices to double
% matrices)
A = double(subs(df_dx_sol)); % Takes the solution for df_dx and creates a matrix of doubles
B = []; % Need to do the same thing for B

% Evaluate eigenvalues of A matrix
eig_A = eig(A)

% Evaluate the controllability
Gamma = []; % create the controllability matrix
rank_Gamma = rank(Gamma)

%% Develop the control law
% Create control (formulate a feedback control matrix K)

% Evaluate equilibrium (check the eigenvalues of the feedback matrix
A_feedback = []; % This is a function of A, B, and K
eig_feedback = []; % The eigenvalues of the feedback matrix

