%% Inverted pendulum variables
% These variables are described at http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling
M = 0.5; % Mass of cart
m = 0.2; % mass of pendulum
b = 0.1; % coefficient of friction for cart
l = 0.3; % length to pendulum center of mass
I = 0.006; % mass moment of inertia of the pendulum
g = 9.8; % Gravity constant

% Symbolic variables for solving for dynamics
syms u real % Force applied to the cart
syms theta real % pendulum angle from vertical
syms thetadot real % time derivative of theta
syms thetaddot real % time derivative of thetadot
syms x real % car position coordinate
syms xdot real % car velocity
syms xddot real % car acceleration

%% Inverted pendulm equations of motion
% Solve for the equations
eqn1 = (M+m)*xddot + b*xdot + m*l*thetaddot*cos(theta) - m*l*thetadot^2*sin(theta) - u;
eqn2 = (I+m*l^2)*thetaddot + m*g*l*sin(theta) + m*l*xddot*cos(theta);
soln = solve(eqn1, eqn2, xddot, thetaddot);
thetaddot = simplify(soln.thetaddot);
xddot = simplify(soln.xddot);

% Create the matlab functions
matlabFunction(thetaddot, 'file', 'Thetaddot.m');
matlabFunction(xddot, 'file', 'Xddot.m');

% Check simplifications
xddot_check = -(100*u - 10*xdot + 147*cos(theta)*sin(theta) + 6 * thetadot^2*sin(theta))/(15*cos(theta)^2 - 70);
xddot_err = simplify(xddot - xddot_check)

s = sin(theta);
c = cos(theta);
thetaddot_check = (343*s + 50*c*u - 5*c*xdot + 3*c*s*thetadot^2)/(3*c^2-14);
thetaddot_err = simplify(thetaddot_check - thetaddot)

