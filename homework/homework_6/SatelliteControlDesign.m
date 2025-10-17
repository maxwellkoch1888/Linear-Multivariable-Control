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
x = [r; theta; rdot; thetadot];
u = [ar; ai];

% State dynamics (create a vector of the dynamics, i.e., f = [rdot; thetadot; rddot; thetaddot] )
f = [x(3);
     x(4);
     -mu/x(1)^2 + u(1) + x(4)^2*x(1);
     (u(2) - 2*x(3)*x(4))/x(1)];

%% Find trajectory to linearize about
R = sym('R', 'real'); % Desired radius of orbit

% Create a trajectory and solve for the equilibrium
% rdot and rddot must be zero (i.e. f(1) and f(3) must be zero) - use the
% subs command
f_sol = subs(f, r, R); % This substitues the big R variable for the little r
f_sol = subs(f_sol, ar, 0); % This substitutes in 0 for ar
f_sol = subs(f_sol, ai, 0); % Substitute 0 in for ai
f_sol = subs(f_sol, f(1), 0); % Substitute 0 for rdot
f_sol = subs(f_sol, f(3), 0); % Substitute 0 for rddot

% Solve for the solution of thetadot (you can use Matlab's "solve" command)
thetadot_sol = solve(f_sol(3)==0, thetadot); % You can do it using two inputs, the third element of f_sol and the symbolic variable for \dot{\theta}
thetadot_sol = thetadot_sol(1)

%% Linear system about the trajectory
% Create the jacobians
df_dx = jacobian(f, x); % Creates a jacobian by evaluating f wrt x
df_du = jacobian(f, u); % Create a jacobian with respect to u

% Create the solution trajectory by defining all of the variables in the
% trajectory - replace all "[]" values with the correct value
w = sym('w', 'real'); % Solution for theta_dot
t = sym('t', 'real'); % Time
r_val = R;
rdot_val = 0;
theta_val = w * t; % System does not depend on theta, still LTI
thetadot_val = w;
ar_val = 0;
ai_val = 0;

% Evaluate linearization at the trajectory solution
df_dx_sol = simplify(subs(df_dx, [r, theta, rdot, thetadot, ar, ai, mu], ...
            [r_val, theta_val, rdot_val, thetadot_val, ar_val, ai_val, mu]));
df_du_sol = simplify(subs(df_du, [r, theta, rdot, thetadot, ar, ai, mu], ...
            [r_val, theta_val, rdot_val, thetadot_val, ar_val, ai_val, mu]));
    
%% Evaluate the linearized system
mu_val = 4.302e-3;
R_val = R_des; % Note, this is passed into the function as a parameter
w_val = double(subs(thetadot_sol, [R, mu], [R_val, mu_val])); % What is the solution

% Get the A and B matrices (need to convert the symbolic matrices to double
% precision matrices)
A = double(subs(df_dx_sol, [R, w, mu], [R_val, w_val, mu_val])) % Takes the solution for df_dx and creates a matrix of doubles
B = double(subs(df_du_sol, [R, w, mu], [R_val, w_val, mu_val])) % Need to do the same thing for B

% Evaluate eigenvalues of A matrix
eig_A = eig(A)

% Evaluate the controllability
Gamma = [B, A*B, A^2*B, A^3*B]; % create the controllability matrix
rank_Gamma = rank(Gamma)

%% Develop the control law
% Create control (formulate a feedback control matrix K)
desired_poles = [-1, -2, -3, -4];
K = place(A, B, desired_poles);

% Evaluate equilibrium (check the eigenvalues of the feedback matrix
A_feedback = A - B*K; % This is a function of A, B, and K
eig_feedback = eig(A_feedback); % The eigenvalues of the feedback matrix

