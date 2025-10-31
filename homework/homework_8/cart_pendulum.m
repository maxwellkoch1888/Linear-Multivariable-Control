function cart_pendulum()
    create_controller();
end

function create_controller()
    % GET STATE SPACE MATRICES
    [A,B,C,D] = get_dynamics();

    
    
end

function [A,B,C,D] = get_dynamics()
% DEFINE PARAMETERS
M = .5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
l = 0.3;

p = I*(M+m)+M*m*l^2;

% BUILD STATE SPACE MODEL 
    A = [0      1              0           0;
         0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
         0      0              0           1;
         0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
    B = [     0;
         (I+m*l^2)/p;
              0;
            m*l/p];
    C = [1 0 0 0;
         0 0 1 0];
    D = [0;
         0];
end