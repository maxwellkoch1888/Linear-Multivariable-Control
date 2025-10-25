function K = CreateSegwayControl

%% Define Variables
% Model Parameters
% mc = sym('mc', 'positive');             % mass of wheel
% ms = sym('ms', 'positive');             % mass of body
% d = sym('d', 'positive');               % distance from C to Cg
% L = sym('L', 'positive');               % half distance between wheels
% R = sym('R', 'positive');               % Radius of wheel
% I2 = sym('I2', 'real');             % intertia about n2
% I3 = sym('I3', 'real');             % inertia about n3
% g = sym('g', 'positive');               % Gravity = 9.8 m/s^2

mc    = .503;         % mass of wheel (kg)
ms    = 4.315;        % mass of body (kg)
d     = .1;           % distance from C to Cg (m)
L     = .1;           % half distance between wheels (m)
R     = .073;         % Radius of wheel (m)
I2    = 3.679*10^-3;  % intertia about n2 (kg m^2)
I3    = 28.07*10^-3;  % inertia about n3 (kg m^2)
g     = 9.8;          % Gravity  (m/s^2)

% States
x = sym('x', 'real'); % x-position of the robot
y = sym('y', 'real'); % y-position of the robot
v = sym('v', 'real');         % velocity in forward direction (i.e. v)
vdot = sym('vdot', 'real');       
phi = sym('phi', 'real');           % angle of shaft from verticle
phidot = sym('phidot', 'real');
phiddot = sym('phiddot', 'real');
psi = sym('psi', 'real');           % Orientation
omega = sym('omega', 'real');
omegadot = sym('omegadot', 'real');

% Controls
u1 = sym('u1', 'real');         % u1 = \alpha + \Beta (torque summation)
u2 = sym('u2', 'real');         % u2 = \alpha - \Beta (torque difference)

%% Solve Equations (each fi = 0)
% Setup Equations
p = sym('p', 'real'); % to substitute in for 1/(2R^2)
f1 = 3*(mc+ms)*vdot - ms*d*cos(phi)*phiddot + ms*d*sin(phi)*(phidot^2+omega^2) + u1/R;
f2 = ( (3*L^2 + 1/(2*R^2))*mc + ms*d^2*sin(phi)^2 + I2 )*omegadot + ms*d^2*sin(phi)*cos(phi)*omega*phidot - L/R*u2;
f3 = ms*d*cos(phi)*vdot + (-ms*d^2-I3)*phiddot + ms*d^2*sin(phi)*cos(phi)*phidot^2 + ms*g*d*sin(phi) - u1;

% Solve Equations
soln = solve(f1, f2, f3, vdot, phiddot, omegadot);
phiddot = soln.phiddot;
omegadot = soln.omegadot;
vdot = soln.vdot;

%% Create Function Calls
% Dynamics
% matlabFunction(phiddot, 'file', 'Phiddot.m');
% matlabFunction(omegadot, 'file', 'Omegadot.m');
% matlabFunction(vdot, 'file', 'Vdot.m');

%% Create the generalized linear equations
% State equations
X = [x; y; psi; omega; v; phi; phidot];
U = [u1; u2];
f = [v*cos(psi); v*sin(psi); omega; omegadot; vdot; phidot; phiddot];

% Calculate the jacobians
df_dx = jacobian(f, X);
df_du = jacobian(f, U);

%% Linearize Equations about x = 0
% Set the states to zero and evaluate the dynamics at the zero state and control
zero_state_f = subs(f, [x, y, psi, omega, v, phi, phidot, u1, u2], [0, 0, 0, 0, 0, 0, 0, 0, 0]);

if zero_state_f == 0
    disp('Zero-state and zero input is an equilibrium point.')
else 
    disp('Zero-state and zero input is not an equilibrium point.')
end
disp(' ')

% Create the state matrices (A,B) from df_dx and df_du
A_sym = df_dx;
B_sym = df_du;
A = subs(df_dx, [x, y, psi, omega, v, phi, phidot, u1, u2], [0, 0, 0, 0, 0, 0, 0, 0, 0]); 
B = subs(df_du, [x, y, psi, omega, v, phi, phidot, u1, u2], [0, 0, 0, 0, 0, 0, 0, 0, 0]);

% Calculate controllability
gamma = ctrb(A, B);
gamma_rank = rank(gamma);
[n, ~] = size(A); 
if gamma_rank == n
    disp('The system is completely controllable.');
else
    disp('The system is not completely controllable.');
end
disp(' ')

%% Linearize Equations about z = 0 (Note, you will want to comment out the previous section so that you 
%%                                  are still working with symbolic variables at this point)
% State equations of reduced state (i.e., define z, u, and dynamics of z state -fz)
z = [omega; v; phi; phidot];
fz = [omegadot; vdot; phidot; phiddot];

% Calculate the jacobians (use the dynamics of the z state to calculate dfz/dz and dfz/du)
dfz_dx = jacobian(fz, z); 
dfz_du = jacobian(fz, U);

% Set the states to zero
% Show that this is an equilibrium
zero_state_fz = subs(fz, [omega, v, phi, phidot, u1, u2], [0, 0, 0, 0, 0, 0]);
if zero_state_fz == 0
    disp('Reduced state at zero-state and zero input is an eq point.')
else 
    disp('Reduced state at zero-state and zero input is not an eq point.')
end 
disp(' ')

% Create the state matrices (A,B) from df_dx and df_du

Az = subs(dfz_dx, [omega, v, phi, phidot, u1, u2], [0, 0, 0, 0, 0, 0]);
Bz = subs(dfz_du, [omega, v, phi, phidot, u1, u2], [0, 0, 0, 0, 0, 0]);

% Evaluate controllability
gamma_reduced_state = ctrb(Az, Bz);
gamma_reduced_rank = rank(gamma_reduced_state);
[n,~] = size(Az);
if gamma_reduced_rank == n 
    disp('The reduced state is completely controllable.')
else
    disp('The reduced state is not completely controllable.')    
end 
disp(' ')

% Explain why z is controllable in isolation from the other states
% Build a T matrix for controllability decomposition
v_T = orth(gamma);
w_T = null(gamma');
T = [v_T, w_T];
T_inv = T^-1;

% Calculate transform
A_hat = T*A*T
B_hat = T*B



%% Create a stabilizing control

