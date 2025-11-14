clear;
clc;

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
v = sym('v', 'real');         % velocity in forward direction (i.e. v)
vdot = sym('vdot', 'real');       
phi = sym('phi', 'real');           % angle of shaft from verticle
phidot = sym('phidot', 'real');
phiddot = sym('phiddot', 'real');
psi = sym('psi', 'real');           % Orientation
omega = sym('omega', 'real');
omegadot = sym('omegadot', 'real');

% Controls
% A = sym('A', 'real');         % torque on wheel 1 (alpha)
% B = sym('B', 'real');         % torque on wheel 2 (Beta)
u1 = sym('u1', 'real');         % u1 = A+B
u2 = sym('u2', 'real');         % u2 = A-B

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

%% Linearize Equations
% State equations
Z = [omega; v; phi; phidot];
U = [u1; u2];
f = [omegadot; vdot; phidot; phiddot];

% State matrix
A_mat = jacobian(f, Z);                     % Linearize
A_mat = subs(A_mat, Z, zeros(4, 1));     % Around X = 0    
A_mat = subs(A_mat, U, zeros(2, 1))      % Around U = 0

% Input matrix
B_mat = jacobian(f, U);                     % Linearize
B_mat = subs(B_mat, Z, zeros(4, 1));     % Around X = 0
B_mat = subs(B_mat, U, zeros(2, 1))      % Around U = 0

%% Create Function Calls
% Dynamics
matlabFunction(phiddot, 'file', 'Phiddot.m');
matlabFunction(omegadot, 'file', 'Omegadot.m');
matlabFunction(vdot, 'file', 'Vdot.m');

% Linear Matrices
matlabFunction(A_mat, 'file', 'AMat_zero.m');
matlabFunction(B_mat, 'file', 'BMat_zero.m');