%% DEFINE THE STATE MATRICES
A = [-0.038,  18.984, 0,     -32.174; ...
     -0.001, -0.632,  1,      0; ...
      0,     -0.759, -0.518,  0; ...
      0,      0,      1,      0];
B = [10.1,   0; ...
     0,     -0.0086; ... 
     0.025, -0.011; ...
     0,      0];
x0 = [10; 0.1; 0.1; 0.0];

% DEFINE Q AND R MATRICES
Q = diag([1, 100, 500, 500]);
R = diag([0.1, 1.0]);

% CHECK CONTROLLABILITY
gamma = ctrb(A,B);
% disp(rank(gamma))

% BUILD GAIN MATRIX
K = lqr(A,B,Q,R);

% BUILD CLOSED LOOP SYSTEM
Acl = A - B*K;

% SIMULATION PARAMETERS
tspan = [0 10];  

% DEFINE ODE FUNCTION
f = @(t,x) Acl * x;

% SIMULATE FORWARD IN TIME
[t, x] = ode45(f, tspan, x0);

% PLOT STATES
figure;
subplot(4,1,1)
plot(t, x(:,1)), ylabel('x_1'), grid on
subplot(4,1,2)
plot(t, x(:,2)), ylabel('x_2'), grid on
subplot(4,1,3)
plot(t, x(:,3)), ylabel('x_3'), grid on
subplot(4,1,4)
plot(t, x(:,4)), ylabel('x_4'), grid on
xlabel('Time (s)')
sgtitle('Closed Loop State Response')