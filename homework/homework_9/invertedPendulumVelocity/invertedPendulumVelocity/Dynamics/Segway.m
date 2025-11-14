classdef Segway < handle
    
    
    properties
        % State indices
        n = 7
        ind_x1 = 1;
        ind_x2 = 2;
        ind_psi = 3;
        ind_omega = 4;
        ind_v = 5;
        ind_phi = 6;
        ind_phidot = 7;
        
        % Feedback indices
        z_ind = 4:7; % Indices of z within x
        z_ind_omega = 1;
        z_ind_v = 2;
        z_ind_phi = 3;
        z_ind_phidot = 4;
        
        % Linearized matrices (\dot{z} = Az + Bu)
        A % State matrix
        B % Input matrix
        
        % Output matrix (y = [omega, v, phi])
        C = [1 0 0 0; 0 1 0 0; 0 0 1 0];
        
        % Feedback matrices
        K % Full state control feedback matrix
        L % Observer gain matrix
        
        % Velocity control variables
        v_d = 0
        omega_d = 0
    end
    
    methods
        function obj = Segway(v_d, omega_d)
            % Store the selected control and associated variables
            obj.v_d = v_d;
            obj.omega_d = omega_d;                 
                        
            % Initialize linearized matrices
            obj.A = AMat_zero();
            obj.B = BMat_zero();
            
            % Create control feedback matrix
            obj.K = obj.createFeedbackControl();
            
            % Create observer matrix
            obj.L = obj.createObserverGains();
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Control Design -- Functions to be implemented
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods         
        % Functions to be implemented
        function K = createFeedbackControl(obj)
            %createFeedbackControl creates the feedback control matrix
            %assuming that the states is given by [omega; v; phi; phidot]
            %with the linearized system matrices (A,B)

            poles = [-5,-2,-3,-4];
            
            % Place holder: zero control
            K = place(obj.A,obj.B, poles);
        end
        
        function [u1, u2] = calculateFeedbackControl(obj, t, x)
            %calculateFeedbackControl Calculates the feedback needed to
            %control the vehicle to the desired velocities
            
            % Extract the control state
            z = x(obj.z_ind);
            
            % Find uff
            uff = [0, 0];

            u_hat = obj.K * z;

            % Placeholder: zero control
            u1 = 0;
            u2 = 0;
            u = u_hat + uff;
            u1 = u(2);
            u2 = u(2);
        end        
    end
    
    %%% Observer functions
    methods
        function L = createObserverGains(obj)
            %createFeedbackControl creates the feedback control matrix
            %assuming that the states is given by [omega; v; phi; phidot]
            %with the linearized system matrices (A,B)            
            
            % Place holder: zero update
            L = zeros(3,4);
        end
        
        function xhat_dot = observerDynamics(obj, xhat, u, y)
            %observerDynamics Calculates the dynamic update for
            %the state defined by [omega; v; phi; phidot] given the output,
            %the current estimate, and the input
            
            % Place holder: zero update
            xhat_dot = zeros(4,1);
        end
        
    end
    
    %%% Dynamics functions  %%%
    methods
        function xdot = dynamics(obj,t, x, u1, u2)
            % Extract needed states
            psi = x(obj.ind_psi);
            omega = x(obj.ind_omega);
            v = x(obj.ind_v);
            phi = x(obj.ind_phi);
            phidot = x(obj.ind_phidot);
            
            % Get unicycle dynamics for updating state
            x1dot = v * cos(psi);
            x2dot = v * sin(psi);
                        
            % Calculate dynamics for \dot{omega}, \ddot{phi}, and \dot{v}
            omegadot = Omegadot(omega, phi, phidot, u2);
            phiddot = Phiddot(omega, phi, phidot, u1);
            vdot = Vdot(omega, phi, phidot, u1);
            
            % Output the aggregate dynamics
            xdot = zeros(obj.n, 1);
            xdot(obj.ind_x1) = x1dot;
            xdot(obj.ind_x2) = x2dot;
            xdot(obj.ind_psi) = omega;
            xdot(obj.ind_omega) = omegadot;
            xdot(obj.ind_v) = vdot;
            xdot(obj.ind_phi) = phidot;
            xdot(obj.ind_phidot) = phiddot;            
        end
        
        function xdot = dynamicsWithObserver(obj, t, x)
            % Extract states
            x_state = x(1:obj.n);
            x_hat = x(obj.n+1:end);
            
            % Calculate control
            [u1, u2] = obj.calculateFeedbackControl(t, x_hat);
            
            % Calculate state dynamics
            x_dot_state = obj.dynamics(t, x_state, u1, u2);
            
            % Calculate the measurement
            y = obj.C * x_state(obj.z_ind); % i.e. C*z
            
            % Calculate z_hat dynamics
            x_hat_dot = obj.observerDynamics(x_hat(obj.z_ind), [u1; u2], y);
            
            % Output the aggregate dynamics (note that the zeros are place
            % holders since we are not estimating the position and
            % orientation)
            xdot = [x_dot_state; zeros(3,1); x_hat_dot];
        end
        
        function xdot = dynamicsWithoutObserver(obj, t, x)            
            % Calculate control
            [u1, u2] = obj.calculateFeedbackControl(t, x);
            
            % Calculate system dynamics
            xdot = obj.dynamics(t, x, u1, u2);            
        end
    end
    
    %%% Plotting functions %%%
    methods
        function plotTiltAndVelocities(obj, t, x, varargin)
           % Extract the state estimate
           if nargin > 3
               xhat = varargin{1};
           else
               xhat = x;
           end
            
           % Plot phi
           subplot(3,1,1);
           plot(t, xhat(obj.ind_phi, :), 'k--', 'linewidth', 3); hold on;
           plot(t, x(obj.ind_phi, :), 'linewidth', 2); 
           plot([t(1), t(end)], [0 0], 'r:', 'linewidth', 2);
           xlabel('Time (s)');
           ylabel('Tilt angle (rad)');
           
           % Plot translational velocity
           subplot(3,1,2);
           plot(t, xhat(obj.ind_v, :), 'k--', 'linewidth', 3); hold on;
           plot(t, x(obj.ind_v, :), 'linewidth', 2); 
           plot([t(1), t(end)], [obj.v_d obj.v_d], 'r:', 'linewidth', 2);
           xlabel('Time (s)');
           ylabel('Trans. Vel (m/s)');
           
           % Plot rotational velocity
           subplot(3,1,3);
           plot(t, xhat(obj.ind_omega, :), 'k--', 'linewidth', 3); hold on;
           plot(t, x(obj.ind_omega, :), 'linewidth', 2); 
           plot([t(1), t(end)], [obj.omega_d obj.omega_d], 'r:', 'linewidth', 2);
           xlabel('Time (s)');
           ylabel('Rot. Vel (rad/s)');
        end
        
        function plot3D(obj, tvec, xvec)
            % Plot the trajectory
            figure;
            plot(xvec(obj.ind_x1, :), xvec(obj.ind_x2, :));
            xlabel('x position', 'fontsize', 14);
            ylabel('y position', 'fontsize', 14);
            title('Path traversed by robot', 'fontsize', 18);

            % Plot the movement of the vehicle
            seg_plot = SegwayPlot();
            for k = 1:10:length(tvec)
                % Extract configuration variables
                q1 = xvec(obj.ind_x1, k); % position variables
                q2 = xvec(obj.ind_x2, k);
                psi = xvec(obj.ind_psi, k); % orientation
                phi = xvec(obj.ind_phi, k); % tilt angle

                % Calculate the transform to the new position and orientation
                T = ObjectPlotter.getTransform(0, 0, psi, [q1; q2; 0]);

                % Plot the segway
                seg_plot.plotSegway(T, phi);
                pause(.1);
            end
        end
    end
end

