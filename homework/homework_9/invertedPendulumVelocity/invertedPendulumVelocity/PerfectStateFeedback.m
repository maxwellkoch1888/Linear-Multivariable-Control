function PerfectStateFeedback()
    clear;
    clc;
    close all;
    
    % Setup correct path
    addpath Components % Used for plotting
    addpath Dynamics % Used for evaluating dynamics
    
    %% Simulate control
    % Create initial conditions
    v0 = 0.5;
    phi0 = pi/4;
    x0 = [0; 0; 0; 0; v0; phi0; 0];
    
    % Create the segway object
    vd = 1.0;
    omegad = 0.25;
    seg = Segway(vd, omegad);
    
    % Simulate control
    t = 0:.01:20;
    [tvec, xvec] = ode45(@(t, x)seg.dynamicsWithoutObserver(t, x), t, x0);
    xvec = xvec'; % reshape to have each state be a column
    
    % Plot results
    figure;
    seg.plotTiltAndVelocities(t, xvec);
    
    % Create a 3D depiction
    seg.plot3D(tvec, xvec);
end

