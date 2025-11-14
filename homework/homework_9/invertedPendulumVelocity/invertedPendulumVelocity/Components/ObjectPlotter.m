classdef ObjectPlotter < handle
    %ObjectPlotter Plots a 3D object with two 2D views
    
    properties
        % Defining the shape
        V % Verticies
        F % Faces
        C % Colors of faces
        
        % Define the plot handles
        handle_3d = [];
        
    end
    
    methods
        function obj = ObjectPlotter(params)
            %% Store parameter variables
            % Defining the shape
            obj.V = params.V;
            obj.F = params.F;
            obj.C = params.C;
            
            %% Plot the object
            % Plot 3d
            obj.handle_3d = patch(params.ax3, 'Vertices', obj.V, 'Faces', obj.F, ...
                'FaceVertexCData', obj.C, 'FaceColor', 'flat');
            title(params.ax3, '3D Plot');
            xlabel(params.ax3, 'x axis');
            ylabel(params.ax3, 'y axis');
            zlabel(params.ax3, 'z axis');
        end
        
        function plotAll(obj, g)
            % Get Euler angles
            [roll, pitch, yaw] = obj.getEuler(g);
            
            % Plot 3d
            obj.draw3d(g);
        end
        
        function draw3d(obj, g)
            %drawbody transforms the verticies of the aircraft by g and
            %replots the aircraft
            
            % Create the Augmented points
            V_aug = [obj.V, ones(size(obj.V, 1), 1)]';
            
            % Transform the points and extract the original points
            Vg = g*V_aug;
            V_new = Vg(1:3,:)';
            
            % Plot the new points
            set(obj.handle_3d, 'Vertices', V_new);
        end
        
        
    end
    
    methods(Static)
        
        function g = getTransform(roll, pitch, yaw, p)
            R_roll = [...
                1, 0, 0; ...
                0, cos(roll), -sin(roll); ...
                0, sin(roll), cos(roll)];
            R_pitch = [...
                cos(pitch), 0, sin(pitch);...
                0, 1, 0;...
                -sin(pitch), 0, cos(pitch)];
            R_yaw = [...
                cos(yaw), -sin(yaw), 0;...
                sin(yaw), cos(yaw), 0;...
                0, 0, 1];

            R = R_yaw*R_pitch*R_roll;
            g = [R, p; zeros(1,3), 1];
        end
        
        function [roll, pitch, yaw] = getEuler(g)
            %getEuler computes the X-Y-Z Euler coordinates (roll,pitch,yaw)
            %from a given transformation matrix
            %
            % The derivation was taken from chapter 2 of "Introduction to
            % Robotics" by John Craig
            R = g(1:3, 1:3);
            B = atan2(-R(3,1), sqrt(R(1,1)^2 + R(2,1)^2 ));
            cB = cos(B);
            alpha = atan2(R(2,1)/cB, R(1,1)/cB);
            gamma = atan2(R(3,2)/cB, R(3,3)/cB);
            
            roll = gamma;
            pitch = B;
            yaw = alpha;
        end
        
        function V = get2DCircleVerticies(r, nop)
            THETA=linspace(0,2*pi,nop);
            RHO=ones(1,nop)*r;
            [X,Y] = pol2cart(THETA,RHO);

            % Create the points
            V = [X', Y'];
        end
    end
end

