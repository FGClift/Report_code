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

% Wake durations to test
W_values = [24, 20, 18, 16];
colors = {'k',[0.8, 0, 0], [0.95, 0.5, 0], [0.1, 0.4, 0]};  

% Number of days to simulate
num_days = 14;

% Pre-initialise results struct for each wake duration group
results = struct();
for w_idx = 1:length(W_values)
    W = W_values(w_idx);
    group_name = sprintf('W%d', W);
    results.(group_name).W = W;
    results.(group_name).time = [];        % Combined time vector
    results.(group_name).performance = []; % Combined p_wake + q_sleep
end

% Create figure
figure;
hold on

for w_idx = 1:length(W_values)
    W = W_values(w_idx);
    group_name = sprintf('W%d', W);
    
    % Time vectors
    dt = 0.01;  % Small for accurate integration
   
      time_wake = (0:dt:W);
    time_sleep = (W:dt:T);
   

    % Circadian process (same for each day due to 24h periodicity)
    C_wake = a*(0.97*sin(omega*(time_wake-theta)) + ...
               0.22*sin(2*omega*(time_wake-theta)) + ...
               0.007*sin(3*omega*(time_wake-theta)) + ...
               0.03*sin(4*omega*(time_wake-theta)) + ...
               0.001*sin(5*omega*(time_wake-theta)));
    
    C_sleep = a*(0.97*sin(omega*(time_sleep-theta)) + ...
                0.22*sin(2*omega*(time_sleep-theta)) + ...
                0.007*sin(3*omega*(time_sleep-theta)) + ...
                0.03*sin(4*omega*(time_sleep-theta)) + ...
                0.001*sin(5*omega*(time_sleep-theta)));
    
    beta_wake = k_circ*C_wake + mu;
    beta_sleep = k_circ*C_sleep + mu;
    
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
    
    % Initial conditions for day 1
    p_wake_onset = 4.49;
    u_wake_onset = 29.9;

    
    % Loop through days
    for day = 1:num_days
        
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
        
        % Sleep period
        v_sleep_onset = u_wake(end) - delta;
        q_sleep_onset = p_wake(end);
        
        q_sleep = zeros(size(time_sleep));
        v_sleep = zeros(size(time_sleep));
        
        for j = 1:length(time_sleep)
            t_rel = time_sleep(j) - W;
            
            state_hom = phi(t_rel)*phi_inverse(0) * [q_sleep_onset; v_sleep_onset];
            
            state_int = zeros(2,1);
            for k = 1:j-1
                s_rel = time_sleep(k) - W;
                integrand = phi(t_rel) * phi_inverse(s_rel) * [beta_sleep(k); 0];
                state_int = state_int + integrand * dt;
            end
            
            state_total = state_hom + state_int;
            q_sleep(j) = state_total(1);
            v_sleep(j) = state_total(2);
        end
       % Plot
day_offset = (day-1)*24;

if W ~= 24 
plot(day_offset + time_wake, p_wake, 'color', colors{w_idx}, 'LineWidth', 2.5) 
plot(day_offset + time_sleep, q_sleep, 'color', colors{w_idx}, 'LineWidth', 2.5)
else 
    % plot wake only up to where day_offset+time_wake < 88
    valid_idx = (day_offset + time_wake) < 88;
    if any(valid_idx)
        plot(day_offset + time_wake(valid_idx), p_wake(valid_idx), 'color', colors{w_idx}, 'LineWidth', 2.5)
    end
end 

      
     
        % Update for next day
        if day < num_days
            p_wake_onset = q_sleep(end);
            u_wake_onset = v_sleep(end) + delta;
        end
    end

    % Plot experimental data 
M = load('vandongenbu.mat');

% Determine fields and plot each as points overlaid on the PVT curves
flds = fieldnames(M);
 data = M.(flds{w_idx});
t = data(:,1);
 val = data(:,2);   


            t_plot = t / 24;

   markers = {'d', 'o', 's', '^'}; % diamond, circle, square, triangle
marker = markers{w_idx};
   
    scatter(t, val, 60, 'Marker',marker, 'MarkerEdgeColor',colors{w_idx}, 'MarkerFaceColor', colors{w_idx},'MarkerFaceAlpha',0.5)
 
end
hold off
x_labels = 0:1:num_days;
set(gca, 'XTick', x_labels * 24, 'FontSize', 40);
set(gca, 'XTickLabel', x_labels, 'FontSize', 40);
hold off
xlabel('Time (days)', 'FontSize', 40) 
ylabel({'Performance'; '(PVT lapses)'}, 'FontSize', 40, 'Position', [-28, 15])
annotation('arrow', [0.045, 0.045], [0.85, 0.92], 'LineWidth', 3); 
annotation('arrow', [0.045, 0.045], [0.24, 0.17], 'LineWidth', 3);
text(-37, 25, 'Worse', 'FontSize', 30, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'Rotation', 0)
text(-37, 5,  'Better', 'FontSize', 30, 'VerticalAlignment', 'top',    'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'Rotation', 0)
ylim([0 30])
yticks(0:5:30)
xlim([0 14*24])
grid on