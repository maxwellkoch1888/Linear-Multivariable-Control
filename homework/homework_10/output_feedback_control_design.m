% BUILD A, B, C, Q, AND R FROM PROBLEM STATEMENT
A = [-5,0,0;0,3,0;0,2,1];
B = [0,0;0,1;1,0];
C = [0,0,1;-1,0,0];

xmax = [1/25,1/9,1];
umax = [1/100, 1];
Q = diag(xmax);
R = diag(umax);

% DETERMINE CONTROLLABILITY
gamma = ctrb(A,B);
disp(rank(gamma))

% PERFORM CONTROLLABLE DECOMPOSITION
T = [orth(gamma), null(gamma')];

Abar = inv(T)*A*T;
Bbar = inv(T)*B;
Qbar = T'*Q*T;

% USE CONTROLLABLE DECOMPOSITION TO FIND k MATRIX
Abar11 = Abar(1:2,1:2);
Bbar1 = Bbar(1:2,1:2);
Qbar11 = Qbar(1:2,1:2);

kbar1 = lqr(Abar11, Bbar1, Qbar11, R);
kbar = [kbar1,[0;0]];
k = kbar*inv(T);

% FIND THE EIGENVALUES OF THE CONTROLLER TO DETERMINE THE PLACEMENT OF
% THE OBSERVER EIGENVALUES
disp(eig(A-B*k))

obs_eig = [-110, -35, -55];
L = place(A', C',  obs_eig);
format long g
disp(L)