classdef SegwayPlot < handle
    
    properties
        % Components
        left_wheel = [] % Holds a class for the left wheel
        right_wheel = [] % Holds a class for the right wheel
        axle = [] % Horizontal shaft connecting the two wheels
        vert_shaft = [] % Vertical shaft  
        
        % Transforms
        T_left_wheel % transform to the left wheel when pendulum in nominal position
        T_right_wheel % transform to the right wheel when pendulum in nominal position
        T_axle % Transform to the axel when the pendulum is in the nominal position
        T_shaft % Transform to the vertical shaft when the pendulum is in the nominal position
        
        % Axes for plotting
        ax3 = [] % 3D axis        
    end
    
    methods
        function obj = SegwayPlot()
            % Create figure for plotting
            figure('units','inches', 'position', [0.5, 0.5, 5, 5], 'color', [250, 250, 250]./255 ) 
            obj.ax3 = gca;
            axis equal;
            hold on;

            % Create the components that define the inverted pendulum
            obj.left_wheel = Wheel(1.0, obj.ax3);
            obj.right_wheel = Wheel(1.0, obj.ax3);
            obj.axle = Shaft(1.0, obj.ax3);
            obj.vert_shaft = Shaft(1.0, obj.ax3);
            
            % Create the nominal transforms that define the inverted pendulum
            obj.T_left_wheel = ObjectPlotter.getTransform(pi/2,0,0,[00; 0.25; 0]);
            obj.T_right_wheel = ObjectPlotter.getTransform(-pi/2,0,0,[0; -0.25; 0]);
            obj.T_axle = ObjectPlotter.getTransform(pi/2,0,0,[0; -0.25; 0]);
            obj.T_shaft = ObjectPlotter.getTransform(0,pi,0,[0; 0; 0]);
            
            % Plot the inverted pendulum at the nominal position
            obj.left_wheel.plotAll(obj.T_left_wheel);
            obj.right_wheel.plotAll(obj.T_right_wheel);
            obj.axle.plotAll(obj.T_axle);
            obj.vert_shaft.plotAll(obj.T_shaft);
            
            % Adjust plot
            light  % Turn on a light source
            material shiny % Make the inverted pendulum interact with the light - be shiny
            view(obj.ax3, -14, 32); % Set orientation of plot
        end

        function plotSegway(obj, T, tilt)
            % Create a tilt transform
            T_tilt = ObjectPlotter.getTransform(0, tilt, 0, zeros(3,1));
            
            % Update transforms
            T_lw = T * obj.T_left_wheel;
            T_rw = T * obj.T_right_wheel;
            T_ax = T * obj.T_axle;
            T_sh = T * T_tilt * obj.T_shaft;
            
            % Update plots
            obj.left_wheel.plotAll(T_lw);
            obj.right_wheel.plotAll(T_rw);
            obj.axle.plotAll(T_ax);
            obj.vert_shaft.plotAll(T_sh);
        end
        
    end
    
    
    
end