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

% %% Linearize Equations about x = 0 (comment out these lines for z-based control design)
% % Set the states to zero
% x = 0;
% y = 0;
% v = 0;
% phi = 0;
% phidot = 0;
% psi = 0;
% omega = 0;
% u1 = 0;
% u2 = 0;
% 
% % Evaluate the dynamics at the zero state and control
% f_zero = subs(f)
% 
% % State matrices
% A_mat = double(subs(df_dx));
% B_mat = double(subs(df_du)); 
% 
% % Calculate controllability
% Gamma = ctrb(A_mat, B_mat);
% rank_Gamma = rank(Gamma);
% 
% % Create a simple transformation matrix
% I = eye(7);
% T = I(:,[4,5,6,7,1,2,3]);
% T_inv = inv(T);
% 
% % Create similarity transform
% A_bar = T_inv * A_mat * T
% B_bar = T_inv * B_mat
% 
% % Look at controllability
% G_cc = ctrb(A_bar(1:4, 1:4), B_bar(1:4, :));
% rank_Gcc = rank(G_cc)

%% Linearize Equations about z = 0 (Comment out above section for control design)
% State equations
Z = [omega; v; phi; phidot];
U = [u1; u2];
fz = f(4:end);

% Calculate the jacobians
dfz_dz = jacobian(fz, Z);
dfz_du = jacobian(fz, U);

% Set the states to zero
x = 0;
y = 0;
v = 0;
phi = 0;
phidot = 0;
psi = 0;
omega = 0;
u1 = 0;
u2 = 0;

% Show that this is an equilibrium
f_eq = subs(fz)

% State matrix
A_mat = double(subs(dfz_dz));
B_mat = double(subs(dfz_du)); 

% Calculate controllability
Gamma = ctrb(A_mat, B_mat);
rank_Gamma = rank(Gamma);

%% Create a stabilizing control
p = [-1,-2,-3,-4];
K = place(A_mat, B_mat, p);
