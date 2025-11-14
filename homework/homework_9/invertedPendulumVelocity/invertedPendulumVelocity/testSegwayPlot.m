function testSegwayPlot()
    close all;
    
    % Adds the components folder which holds the plotting objects
    addpath Components
    
    % Initialize the plotting
    seg = SegwayPlot();
        
    % Initialize variables for plotting a circle
    n = 1000;   % Number of points to plot
    V = ObjectPlotter.get2DCircleVerticies(10, n)';  % Verticies of the points
    t = 0; % used for time indexing
    
    % Move the segway around in a circle
    for k = 1:n-1
       % Calculate orientation by looking at the next point
       q_c = V(:,k);
       q_n = V(:,k+1);
       q_diff = q_n-q_c;
       yaw = atan2(q_diff(2), q_diff(1));
       
       % Calculate the tilt angle as a sinusoid
       dt = 0.1;
       tilt = sin(t);
       t = t + dt;
       
       
       % Calculate a transform to the new position and orientation
       T = ObjectPlotter.getTransform(0,0,yaw, [q_c;0]);
       
       % Plot the segway
       seg.plotSegway(T, tilt);
       
       pause(0.1);
    end
    
    

end

