function StateSpaceIntegralControl()
    % Run
    close all
    
    figure;
    % System A
    x0 = [0;0];
    [K_A_int, K_A_no_int, u_ff_A, x_d_A, A_true_A, B_true_A, G_A] = system_A_gain();
    simulate_system(A_true_A, B_true_A, G_A, K_A_int,   u_ff_A, x_d_A, x0, true);
    simulate_system(A_true_A, B_true_A, G_A, K_A_no_int,u_ff_A, x_d_A, x0, false);

    figure;
    % System B
    x0 = [1;0];
    [K_B_int, K_B_no_int, u_ff_B, x_d_B, A_true_B, B_true_B, G_B] = system_B_gains();
    simulate_system(A_true_B, B_true_B, G_B, K_B_int,   u_ff_B, x_d_B, x0, true);
    simulate_system(A_true_B, B_true_B, G_B, K_B_no_int,u_ff_B, x_d_B, x0, false);

    %% Functions
    % System A
    function [k_int, k_no_int, u_ff, x_d, A_true, B_true, G] = system_A_gain()
        disp('System A')
    
        % Define system
        A = [0, -2; 2, 0];
        B = [0; 1];
    
        % Evaluate controllability of the augmented system
        G = [1 0];
        Abar = [A, zeros(2,1); G, 0];
        Bbar = [B; 0];
        Gamma = ctrb(Abar,Bbar);
        rank_Gamma = rank(Gamma); % Rank = 3 => CC
    
        % Calculate the control
        Q = diag([1/3^2, 1/4^2, 1/10^2]);
        R = 1/2^2;
        k_int = lqr(Abar, Bbar, Q, R);

        % Calculate the steady state value using the true dynamics
        A_true = [0, -2.1; 2.01, 0];
        B_true = [0; 0.95];
        A_aug = [A_true, zeros(2,1); G 0];
        B_aug = [B_true; 0];
        x_d = [-2; 0];
        u_ff = 4;
        % x_ss = -inv(A_aug - B_aug*K)*(B_aug*K*[x_d; 0] + B_aug*u_ff - [zeros(2,1); G*x_d]);
    
        % Calculate the steady state without the integral
        k_no_int = lqr(A, B, Q(1:2, 1:2), R);
        % x_ss_err = -inv(A_true - B_true*K)*(B_true*K*x_d + B*u_ff);
    end 

    % System B
    function [k_int, k_no_int, u_ff, x_d, A_true, B_true, G] = system_B_gains()
        disp('System B')
    
        % Define system
        A = [-1, 0; 0, 1];
        B = [0; 1];
    
        % Evaluate controllability of the augmented system
        G = [0 1];
        Abar = [A, zeros(2,1); G, 0];
        Bbar = [B; 0];
        Gamma = ctrb(Abar,Bbar);
        rank_Gamma = rank(Gamma); % Rank = 2 => not CC
    
        % Test for stabilizability -- Because all rhp eigenvalues are
        % controllable, then the system is stabilizable
        eig_A = eig(Abar); % Eigenvalues are [-1, 0, 1]
        rank_0 = rank([-Abar, Bbar]); % rank = 3
        rank_1 = rank([eye(3)-Abar, Bbar]); % rank = 3
    
        % Calculate the control - wonderful cheat for those that look at the
        % solutions - matlab's lqr function works on stabilizable systems!
        Q = diag([1/6^2, 1/1^2, 1/3^2]);
        R = 1/4^2;
        k_int = lqr(Abar, Bbar, Q, R)

        % Calculate the steady state value using the true dynamics
        A_true = [-0.95, 0; 0, 1.1];
        B_true = [0; 0.97];
        A_aug = [A_true, zeros(2,1); G 0];
        B_aug = [B_true; 0];
        x_d = [0; 3];
        u_ff = -3;
        % x_ss = -inv(A_aug - B_aug*K)*(B_aug*K*[x_d; 0] + B_aug*u_ff - [zeros(2,1); G*x_d]);
    
        % Calculate the steady state without the integral
        disp('A')        
        disp(A)
        disp('B)')
        disp(B)
        disp('Q')
        disp(Q(1:2, 1:2))
        disp('R')
        disp(R)
        k_no_int = lqr(A, B, Q(1:2, 1:2), R)
        % x_ss_err = -inv(A_true - B_true*K)*(B_true*K*x_d + B*u_ff);
    end

    % Simulate a system
    function simulate_system(A_true, B_true, G, K, u_ff, x_d, x0, use_integrator)
        T = 20;
        i0 = 0;   % integrator state
    
        if use_integrator
            x0_aug = [x0; i0];
            f = @(t,x) dyn_with_integrator(t,x,A_true,B_true,G,K,u_ff,x_d);
            n = 3;
        else
            x0_aug = x0;
            disp('A_true')
            disp(A_true)
            disp('B_true')
            disp(B_true)
            disp('k')
            disp(K)
            disp('uff')
            disp(u_ff)
            disp('xd')
            disp(x_d)

            f = @(t,x) dyn_without_integrator(t,x,A_true,B_true,K,u_ff,x_d);
            n = 2;
        end
    
        [t, X] = ode45(f, [0 T], x0_aug);
    
        % Compute input u(t)
        U = zeros(length(t),1);
        for k=1:length(t)
            if use_integrator
                x = X(k,1:2)';
                i = X(k,3);
                U(k) = -K * [x; i] + u_ff;
            else
                x = X(k,1:2)';
                U(k) = -K * x + u_ff;
            end
        end
    
        % Plot
        sgtitle('States and Control Input')
            
        subplot(4,1,1)
        plot(t, X(:,1), LineWidth=1.5); 
        ylabel('x1')
        xlabel('Time (s)')
        hold on 
        
        subplot(4,1,2)
        plot(t, X(:,2), LineWidth=1.5); 
        ylabel('x2')
        xlabel('Time (s)')
        hold on 

        if (use_integrator == true)
            subplot(4,1,3)
            plot(t,X(:,3), LineWidth=1.5)
            ylabel('x3')
            xlabel('Time (s)')
            hold on
        end 
    
        subplot(4,1,4)
        plot(t, U,'LineWidth',1.6)
        ylabel('Input u')
        hold on 
    end
    
    % Build dynamics if there is an integrator
    function dx = dyn_with_integrator(~,x,A,B,G,K,u_ff,x_d)
        x2 = x(1:2);
        i  = x(3);
        u  = -K*[x2; i] + u_ff;
        dx = [A*x2 + B*u;
              G*(x2 - x_d)];
    end

    % Build dynamics if there is no integrator
    function dx = dyn_without_integrator(~,x,A,B,K,u_ff,x_d)
        u = -K *(x - x_d) + u_ff;
        dx = A*x + B*u;
    end
end

