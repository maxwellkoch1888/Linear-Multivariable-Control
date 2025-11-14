classdef Cylinder < ObjectPlotter
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        r % radius 
        h % height  
        nop = 20 % number of points
    end
    
    methods
        function obj = Cylinder(scale, ax3, c_top, c_bottom, c_side, r, h, nop)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Inputs:
            %   - scale: defines the size of the UAV (i.e. 1 is actual, > 1
            %   magifies, and < 1 shrinks)
            %   - ax3: axis of 3D plot
            %   - axtop: axis of top down view
            %   - axside: axis of side view
            %   - c_top: color of top
            %   - c_bottom: color of bottom of cylinder
            %   - c_side: color of side of cylinder
            %   - r: radius of cylinder
            %   - h: height of cylinder
            %   - nop: number of points used to create one of the circules
            %   of the cylinder
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% Define the top down 2D shape
            params = ObjectParam;
            V_top = ObjectPlotter.get2DCircleVerticies(r, nop) .* scale;
            V_top = [V_top zeros(nop, 1)];
            
            %% Define the 3D shape
            % Create the Vertices
            params.V = zeros(nop*2, 3);
            params.V(1:nop, :) = V_top;
            params.V(nop+1:end, :) = V_top;
            params.V(nop+1:end, 3) = -h.*ones(nop, 1);
            
            % Create the faces and colors
            pad = NaN.*zeros(1,nop-4); % Used to pad face definitions
            params.F = [...
                1:nop;...       % f1: top 
                nop+1:(2*nop)]; % f2: bottom
            params.C = [...
                c_top;...       % Top color
                c_bottom];      % Bottom color
            for k = 1:nop-1     
                % Create a side face
                params.F = [params.F; ...
                    k, k+1, k+nop+1, k+nop, pad];
                
                % Add a side color
                params.C = [params.C; c_side];
            end
            
            %% Create the object
            % Call superclass
            params.ax3 = ax3;
            obj@ObjectPlotter(params);
            
            % Store variables
            obj.r = r;
            obj.h = h;
            obj.nop = nop;
            
        end
    
        
    end
end

