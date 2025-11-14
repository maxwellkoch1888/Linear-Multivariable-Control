classdef LuenbergerObserver < handle
    %LuenbergerObserver A simple observer for a linear system
    
    properties
        % State information
        xhat % state estimate
        pend % instance of the SimplePendulum class
        
        % Plant model
        A % Plant matrix
        B % Input matrix
        C % Output matrix
        
        % Observer design
        poles % Vector of desired poles
        L % Observer gain matrix
    end
    
    methods
        function obj = SimpleObserver(poles, xhat, A, B, C)
            %SimpleObserver Construct an instance of this class
            
            % Store variables
            obj.poles = poles;
            obj.xhat = xhat;
            obj.A = A;
            obj.B = B;
            obj.C = C;
            
            % Create observer
            obj.calculateObserverGains();
        end
        
        function calculateObserverGains(obj)
            % Check for observability
            observability_mat = obsv(obj.A, obj.C);
            if rank(observability_mat) == size(obj.A, 1)
                disp('System is observable');
            else
                disp('System is not observable');                
            end  
            
            % Place poles - create observer gain
            obj.L = ( place(obj.A', obj.C', obj.poles) )';            
        end
        
        function xhat_dot = dynamics(obj, xhat, u, y)
            xhat_dot = obj.A*xhat + obj.B*u - obj.L*(obj.C*xhat - y);
        end
    end
end

