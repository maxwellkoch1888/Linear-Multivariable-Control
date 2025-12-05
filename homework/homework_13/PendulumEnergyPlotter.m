classdef PendulumEnergyPlotter < handle
    %PendulumEnergyPlotter Plots a pendulum and its corresponding energy
    %side-by-side
    
    properties
        l % Length of the pendulum  
        h_shaft = [] % Handle to the pendulum shaft
        h_mass = [] % Handle to the pendulum mass
        h_energy = [] % Handle to the energy plot
    end
    
    methods
        function obj = PendulumEnergyPlotter(pendulum_length)
            %PendulumEnergyPlotter Construct an instance of this class
            %   pendulum_length: length of the pendulum from (0,0)            
            obj.l = pendulum_length;
            
            % Initialize the plot
            obj.initializePlot();            
        end
        
        function initializePlot(obj)
            %% Pendulum figure
            % Plot the shaft
            figure('units', 'inches', 'position', [1 1 12 5]);
            subplot(1,2,1); 
            p1 = [0;0]; 
            p2 = [0; obj.l]; 
            obj.h_shaft = plot([p1(1) p2(1)], [p1(2) p2(2)], 'k', 'linewidth', 2);
            hold on;
            
            % Plot the mass
            obj.h_mass = plot(p2(1), p2(2), 'ro', 'linewidth', 10);
            
            % Set the axis for the figure
            axis equal
            set(gca, 'xlim', [-obj.l, obj.l]);
            set(gca, 'ylim', [-obj.l, obj.l]);
            
            %% Energy figure
            subplot(1,2,2);
            obj.h_energy = plot(0, 0, 'r', 'linewidth', 3);
            xlabel('Time (sec)');
            ylabel('Energy');
            set(gca, 'fontsize', 18);
            
        end
        
        function plot(obj, theta, energy, time)
            %plot: Plots the pendulum and energy given a pendulum angle (theta)
            % the energy (energy) and the time for the given angle and
            % energy
            
            % Calculate position
            x_pos = obj.l*sin(theta);
            y_pos = obj.l*cos(theta);
            
            % Plot the shaft
            p1 = [0;0]; 
            p2 = [x_pos, y_pos];
            set(obj.h_shaft, 'xdata', [p1(1) p2(1)], 'ydata', [p1(2) p2(2)]);
            
            % Plot the mass
            set(obj.h_mass, 'xdata', p2(1), 'ydata', p2(2));
            
            % Plot the energy
            ydata = [get(obj.h_energy, 'ydata') energy];
            xdata = [get(obj.h_energy, 'xdata') time];
            set(obj.h_energy, 'xdata', xdata, 'ydata', ydata);
        end
    end
end

