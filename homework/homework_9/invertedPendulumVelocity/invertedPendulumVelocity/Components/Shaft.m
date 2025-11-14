classdef Shaft < Cylinder
    methods
        function obj = Shaft(scale, ax3)
            % Define shaft color
            c_shaft = [188, 198, 204] ./255;   % Dark orange
            
            % Define dimensions
            r = .1;
            h = 0.5; 
            nop = 20;
            
            % Create the buoy
            obj@Cylinder(scale, ax3, ...
                c_shaft, c_shaft, c_shaft, r, h, nop);
        end
    end
end

