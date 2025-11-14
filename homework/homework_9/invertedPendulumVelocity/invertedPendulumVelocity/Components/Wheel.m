classdef Wheel < Cylinder
    methods
        function obj = Wheel(scale, ax3)
            % Define the wheel color (black)
            c_wheel = [105, 105, 105]./255;
            
            % Define dimensions
            r = 0.2;
            h = .1; 
            nop = 20;
            
            % Create the buoy
            obj@Cylinder(scale, ax3, ...
                c_wheel, c_wheel, c_wheel, r, h, nop);
        end
    end
end

