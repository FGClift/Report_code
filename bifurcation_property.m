%% 
% Parameters from equation (27)
alpha_11 = -0.0135;
alpha_12 = 0.000929;
alpha_22 = 0.00743;
sigma_11 = -2.17;
sigma_12 = 0.872;
sigma_22 = -0.0397;
delta = 19.8;
T = 24;  % Total cycle duration
tau = 24; % Circadian oscillation phase

% Circadian parameters
omega = 2*pi/tau;
theta = 12.7;
a = 0.12;
k_circ = 5.86;
mu = 0.472;




% Number of days to simulate
num_days = 4;

    
    % Initial conditions for day 1
    p_wake_onset = 4.49;
    u_wake_onset = 29.9;

% Create figure
figure;
hold on


    W = 24;

    
    % Time vectors
    dt = 0.01;  % Small for accurate integration
  
    t_asymptote= (0:dt:88);
      time_wake = (0:dt:W);


    %% Asymptote 

    upper_asymptote= u_wake_onset*exp(alpha_22*t_asymptote);

    u0=p_star*ones(size(t_asymptote)); % p_star is the value obtained from eqs_17 with W=24
    
 %% Model simulation

    % Circadian process (same for each day due to 24h periodicity)
    C_wake = a*(0.97*sin(omega*(time_wake-theta)) + ...
               0.22*sin(2*omega*(time_wake-theta)) + ...
               0.007*sin(3*omega*(time_wake-theta)) + ...
               0.03*sin(4*omega*(time_wake-theta)) + ...
               0.001*sin(5*omega*(time_wake-theta)));
    

    
    beta_wake = k_circ*C_wake + mu;

    
    % Construct fundamental solutions (equation 14)
    k11 = 1; k12 = 0;
    k21 = alpha_12/(alpha_22 - alpha_11); k22 = 1;
    
    psi = @(t) [k11*exp(alpha_11*t), k21*exp(alpha_22*t); 
                k12*exp(alpha_11*t), k22*exp(alpha_22*t)];
    psi_inverse= @(t)1/exp((alpha_11+alpha_22)*t)*[k22*exp(alpha_22*t), -k21*exp(alpha_22*t); 
                                                  -k12*exp(alpha_11*t), k11*exp(alpha_11*t)];
    
    k31 = 1; k32 = 0;
    k41 = sigma_12/(sigma_22 - sigma_11); k42 = 1;
    
    phi = @(t) [k31*exp(sigma_11*t), k41*exp(sigma_22*t); 
                k32*exp(sigma_11*t), k42*exp(sigma_22*t)];
    phi_inverse = @(t) 1/exp((sigma_11+sigma_22)*t)*[k42*exp(sigma_22*t), -k41*exp(sigma_22*t); 
                                                       -k32*exp(sigma_11*t), k31*exp(sigma_11*t)];

 
    % Loop through days
    % Preallocate storage for diamond markers
    diamond_times_all = zeros(1, num_days+1);
    diamond_perfs_all = zeros(1, num_days+1);


figure('Position',[200 400 1400 1000]);
    for day = 1:num_days+1

        % Wake period
        p_wake = zeros(size(time_wake));
        u_wake = zeros(size(time_wake));

        for j = 1:length(time_wake)
            t = time_wake(j);

            state_hom = psi(t)*psi_inverse(0) * [p_wake_onset; u_wake_onset];

            state_int = zeros(2,1);
            for k = 1:j-1
                s = time_wake(k);
                integrand = psi(t)*psi_inverse(s) * [beta_wake(k); 0];
                state_int = state_int + integrand * dt;
            end
            state_total = state_hom + state_int;
            p_wake(j) = state_total(1);
            u_wake(j) = state_total(2);
        end

        % Performance at beginning of this day (time = 0 of time_wake)
        diamond_performance = p_wake(1);
        square_performance = u_wake(1);
        diamond_time = (day-1)*24;

        % store this to diamond_times_all every day
        diamond_times_all(day) = diamond_time;
        diamond_perfs_all(day) = diamond_performance;

        % Plot
        day_offset = (day-1)*24;

    
        plot(day_offset + time_wake, p_wake, 'k', 'LineWidth', 2)
        hold on
        % Plot diamond at beginning of each day
        plot(diamond_time, diamond_performance, 'd', 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', 'k', 'MarkerSize', 8)
        hold on
        plot(diamond_time, square_performance, 's', 'MarkerEdgeColor', [0.5 0.5 0.5], ...
            'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
        hold on
        % Connect diamond points up to current day
        idx = 1:day;
        plot(diamond_times_all(idx), diamond_perfs_all(idx), '--k', 'LineWidth', 1.5);

        % Update for next day
        if day <= num_days
            p_wake_onset = p_wake(end);
            u_wake_onset = u_wake(end);
        end

    end

    hold on
    plot(t_asymptote,upper_asymptote,'Color', [0.5 0.5 0.5], 'LineWidth', 3)
    hold on 
    plot(t_asymptote, u0,'Color', [0.5 0.5 0.5],'Linewidth', 1.5 ,'LineStyle','--') 
    
   
x_labels = 0:1:num_days;
set(gca, 'XTick', x_labels * 24, 'FontSize', 40);
set(gca, 'XTickLabel', x_labels, 'FontSize', 40);
hold off
xlabel('Time (days)', 'FontSize', 40)
ylabel('Performance', 'FontSize', 40, 'Position', [-8, 30])
annotation('arrow', [0.04, 0.04], [0.76, 0.83], 'LineWidth', 3); 
annotation('arrow', [0.04, 0.04], [0.3, 0.23], 'LineWidth', 3);
text(-10, 45, 'Worse', 'FontSize', 25, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'Rotation', 0)
text(-10, 15,  'Better', 'FontSize', 25, 'VerticalAlignment', 'top',    'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'Rotation', 0)
ylim([0 60])
yticks(0:10:60)
xlim([0 88])
grid on