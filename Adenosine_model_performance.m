
% Parameters
mu_s = 596.4;
mu_w = 869.5;
chi_s = 4.2;
chi_w = 18.18;
Kd1 = 1; % nM
Kd2 = 100; % nM
gamma = 0.9677;
lambda = 291;
R2u = 300; % nM (approximately R2,tot since Kd2 is large)
beta = R2u/(R2u + Kd2);
a = 3.25;
omega = 2*pi/24;
phi = 7.95;
p_max = 60;
D_mid = 583.2;
D_s = 5.872;

wake_time = 8; % All conditions wake up at 8am
T_baseline = 7; % Baseline sleep condition
sleep_start = wake_time - T_baseline; % 1am for baseline
t0 = wake_time;
dt = 0.1; % small timestep for initialization
t = 0:dt:24*13; % Day 0 starts at wake time, simulate 14 days

%% Proper initialization to equilibrium 
T_prestudy = 7.8; % hours sleep prior to inpatient schedule
sleep_start_prestudy = wake_time - T_prestudy; % ensures wake at 8am

% Function to simulate one 24-hour cycle
simulate_one_day = @(Atot_init, R1tot_init, T_sleep, sleep_start_time) ...
    simulate_day_helper(Atot_init, R1tot_init, T_sleep, sleep_start_time, ...
    mu_s, mu_w, chi_s, chi_w, Kd1, beta, gamma, lambda, dt);

%% Initial conditions 
Atot_eq = mu_s + 0.6237*(mu_w - mu_s);
Atot_mean = mu_s + 0.302*(mu_w - mu_s);
R1tot_eq = Atot_mean/gamma - Kd1/((1-gamma)*(1-beta));

fprintf('Initial guesses (from paper formula):\n');
fprintf('Atot_eq = %.2f\n', Atot_eq);
fprintf('Atot_mean = %.2f\n', Atot_mean);
fprintf('R1tot_eq = %.2f\n\n', R1tot_eq);

%% Iterate to equilibrium with M-factor (like paper)
mfactor = 20;
lambda_fast = lambda/mfactor;

tolerance = 1e-6;
max_iterations = 100;

fprintf('Phase 1: Fast convergence with m-factor = %d\n', mfactor);
for iter = 1:max_iterations
    [Atot_new, R1tot_new] = simulate_one_day_mfactor(Atot_eq, R1tot_eq, T_prestudy, ...
        sleep_start_prestudy, mu_s, mu_w, chi_s, chi_w, Kd1, beta, gamma, lambda_fast, dt);
    
    if abs(Atot_new - Atot_eq) < tolerance && abs(R1tot_new - R1tot_eq) < tolerance
        fprintf('Converged after %d iterations (with m-factor)\n', iter);
        break;
    end
    
    Atot_eq = Atot_new;
    R1tot_eq = R1tot_new;
    
    if mod(iter, 20) == 0
        fprintf('  Iteration %d: Atot=%.2f, R1tot=%.2f\n', iter, Atot_eq, R1tot_eq);
    end
end

fprintf('After fast phase: Atot=%.2f, R1tot=%.2f\n\n', Atot_eq, R1tot_eq);

%% Phase 2: Fine convergence without m-factor
fprintf('Phase 2: Fine convergence with normal lambda\n');
for iter = 1:max_iterations
    [Atot_new, R1tot_new] = simulate_one_day(Atot_eq, R1tot_eq, T_prestudy, sleep_start_prestudy);
    
    if abs(Atot_new - Atot_eq) < tolerance && abs(R1tot_new - R1tot_eq) < tolerance
        fprintf('Converged after %d iterations (normal lambda)\n', iter);
        break;
    end
    
    Atot_eq = Atot_new;
    R1tot_eq = R1tot_new;
    
    if mod(iter, 20) == 0
        fprintf('  Iteration %d: Atot=%.2f, R1tot=%.2f\n', iter, Atot_eq, R1tot_eq);
    end
end

fprintf('After fine phase: Atot=%.2f, R1tot=%.2f\n\n', Atot_eq, R1tot_eq);

%% Step 3: Simulate 3 baseline nights with 7.0h sleep
n_baseline_nights = 3;

for night = 1:n_baseline_nights
    [Atot_eq, R1tot_eq] = simulate_one_day(Atot_eq, R1tot_eq, T_baseline, sleep_start);
end

%% Step 4: Get state at wake-up time
[initial_Atot, initial_R1tot] = simulate_to_wake_helper(Atot_eq, R1tot_eq, T_baseline, sleep_start, ...
    mu_s, mu_w, chi_s, chi_w, Kd1, beta, gamma, lambda, dt);

fprintf('Completed 3 baseline nights (7.0h sleep)\n');
fprintf('Final initial conditions: Atot=%.2f, R1tot=%.2f\n\n', initial_Atot, initial_R1tot);

%% Simulate for different TIB conditions
T_conditions = [0, 4, 6, 8]; % hours of time in bed
p_chronic_all = cell(1, 4);
actual_wake_all = cell(1, 4);

for cond = 1:length(T_conditions)
    T_current = T_conditions(cond);
    
    if T_current == 0
        sleep_start_chronic = 0;
    else
        sleep_start_chronic = wake_time - T_current;
        if sleep_start_chronic < 0
            sleep_start_chronic = sleep_start_chronic + 24;
        end
    end
    sleep_end_chronic = wake_time;
    
    if T_current ~= 0
        fprintf('Condition %dh TIB: Sleep %02d:00-%02d:00\n', T_current, ...
            mod(sleep_start_chronic, 24), sleep_end_chronic);
    end
    
    % Initialize arrays
    Atot_chronic = zeros(size(t));
    R1tot_chronic = zeros(size(t));
    R1b_chronic = zeros(size(t));
    D_chronic = zeros(size(t));
    p_chronic = zeros(size(t));
    wake = zeros(size(t));
    
    R1tot_chronic(1) = initial_R1tot;
    Atot_chronic(1) = initial_Atot;
    
    for i = 1:length(t)-1
        day_time = mod(t(i)+t0, 24);
        
        if T_current == 0
            wake(i) = 1;
        elseif sleep_start_chronic < wake_time
            wake(i) = ~(day_time >= sleep_start_chronic && day_time < wake_time);
        else
            wake(i) = ~(day_time >= sleep_start_chronic || day_time < wake_time);
        end
        
        if wake(i) == 0
            Atot_chronic(i+1) = mu_s + (Atot_chronic(i) - mu_s)*exp(-dt/chi_s);
        else
            Atot_chronic(i+1) = mu_w + (Atot_chronic(i) - mu_w)*exp(-dt/chi_w);
        end
        
        R1b_chronic(i) = 0.5*(Atot_chronic(i) + R1tot_chronic(i) + Kd1/(1-beta) - ...
            sqrt((Atot_chronic(i) + R1tot_chronic(i) + Kd1/(1-beta))^2 - 4*Atot_chronic(i)*R1tot_chronic(i)));
        
        dR1tot = (R1b_chronic(i) - gamma*R1tot_chronic(i))/lambda;
        R1tot_chronic(i+1) = R1tot_chronic(i) + dR1tot*dt;
        
        D_chronic(i) = R1b_chronic(i) + a*cos(omega*(t(i)+t0-phi));
        p_chronic(i) = p_max/(1+exp((D_mid-D_chronic(i))/D_s));
    end
    
    % Final point
    i = length(t);
    R1b_chronic(i) = 0.5*(Atot_chronic(i) + R1tot_chronic(i) + Kd1/(1-beta) - ...
        sqrt((Atot_chronic(i) + R1tot_chronic(i) + Kd1/(1-beta))^2 - 4*Atot_chronic(i)*R1tot_chronic(i)));
    D_chronic(i) = R1b_chronic(i) + a*cos(omega*(t(i)+t0-phi));
    p_chronic(i) = p_max/(1+exp((D_mid-D_chronic(i))/D_s));
    
    p_chronic_all{cond} = p_chronic;
    actual_wake_all{cond} = wake;
    
  
end

%% Plotting with 2-hour sleep inertia filter

dt_plot = 2; 
t_plot_indices = 1:round(dt_plot/dt):length(t);
t_plot_days = t(t_plot_indices)/24;
sleep_inertia_hours = 2; 

figure('Units','normalized','Position',[0 0 1 1]); 

% Reduce padding so subplots occupy more of the figure
pad = 0.06;           % smaller padding -> larger subplot area
w = (1 - 3*pad) / 2;  % adjusted to match padding layout
h = (1 - 5*pad) / 2; 

annotation_labels = {'A','B','C','D'};

for cond = 1:4
    % Determine subplot position
    if cond == 1; pos = [pad, pad + h + pad, w, h];
    elseif cond == 2; pos = [pad + w + pad, pad + h + pad, w, h];
    elseif cond == 3; pos = [pad, pad, w, h];
    else; pos = [pad + w + pad, pad, w, h];
    end
    
    subplot('Position', pos)
    hold on
    
    % Plot sleep period
    wake_vec = actual_wake_all{cond};
    % Find transitions: 1 to 0 (sleep start) and 0 to 1 (sleep end)
    sleep_starts = find(diff(wake_vec) == -1);
    sleep_ends = find(diff(wake_vec) == 1);
    
    % Handle edge cases (if simulation starts/ends in sleep)
    if wake_vec(1) == 0; sleep_starts = [1, sleep_starts]; end
    if wake_vec(end) == 0; sleep_ends = [sleep_ends, length(wake_vec)]; end
    
    for s = 1:length(sleep_starts)
        t_start = t(sleep_starts(s))/24;
        t_end = t(sleep_ends(s))/24;
       
        patch([t_start t_end t_end t_start], [0 0 40 40], [0.8 0.8 0.8], ...
              'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

     % Plot experimental data 
M = load('vandongenbu.mat');

% Determine fields and plot each as points overlaid on the PVT curves
flds = fieldnames(M);
 data = M.(flds{cond});
t_data = data(:,1);
 val = data(:,2);   
   
    scatter(t_data/24, val, 'Marker', '^','MarkerEdgeColor','k', 'MarkerFaceColor', 'none')
    % --------------------------------
    
    % Plot ODE results (assuming 'results' exists in your workspace)
    % Map indices to your results structure
    res_fields = {'W24', 'W20', 'W18', 'W16'};
    res_t = results.(res_fields{cond}).time/24;
    res_p = results.(res_fields{cond}).performance;

    % Exclude any data points where time/24 == 13 or time/24 == 0 (avoid plotting those)
    exclude_idx = (res_t == 13) | (res_t == 0);
    if any(exclude_idx)
        res_t(exclude_idx) = [];
        res_p(exclude_idx) = [];
    end
     
  
    plot(res_t, res_p, '.', 'Color', 'b', 'MarkerSize', 20)
 

    % Filter data for plotting (2-hour inertia)
    valid_idx = false(size(t_plot_indices));
    for i = 1:length(t_plot_indices)
        time_since_wake = mod(t(t_plot_indices(i)), 24);
        valid_idx(i) = actual_wake_all{cond}(t_plot_indices(i)) == 1 && ...
                       time_since_wake >= sleep_inertia_hours;
    end
    
    % Plot simulated data
    plot(t_plot_days(valid_idx), p_chronic_all{cond}(t_plot_indices(valid_idx)), ...
        '.', 'Color','r', 'MarkerSize', 20)
    
    
    % Formatting
    set(gca, 'Layer', 'top'); % Ensures grid/axes stay on top of the grey bars
    set(gca, 'FontSize', 43); 
    if cond>2
    xlabel('Time (Days)', 'FontSize', 50)
    end 
    if mod(cond, 2) ~= 0 % Only left subplots get Y-label
        ylabel({'Performance'; '(PVT lapses)'}, 'FontSize', 50)
    end
    ylim([0 40])
    yticks(0:10:40)
    grid on
    
    if cond == 1
        xlim([0 4]); xticks(0:1:4);
        legend('Data (group averages)','McCauley et al.','Phillips et al.', 'Fontsize',32)
    else
        xlim([0 13]); xticks(0:1:13);
    end

    % Add annotation letter in top-left corner of the subplot
    % Compute normalized position for annotation relative to the figure
    % Use a small offset inside the subplot
    x_ann = pos(1) + 0.01*pos(3);
    y_ann = pos(2) + pos(4) - 0.15*pos(4);
    annotation('textbox', [x_ann, y_ann, 0.05, 0.05], 'String', annotation_labels{cond}, ...
        'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 45);
end

%% Helper functions
function [Atot_end, R1tot_end] = simulate_day_helper(Atot_init, R1tot_init, T_sleep, sleep_start_time, ...
    mu_s, mu_w, chi_s, chi_w, Kd1, beta, gamma, lambda, dt)
    
    t_day = 0:dt:24;
    Atot = Atot_init;
    R1tot = R1tot_init;
    
    for i = 1:length(t_day)-1
        day_time = mod(t_day(i), 24);
        
        if sleep_start_time < 0
            sleep_start_time = sleep_start_time + 24;
        end
        
        if sleep_start_time + T_sleep <= 24
            is_awake = ~(day_time >= sleep_start_time && day_time < sleep_start_time + T_sleep);
        else
            is_awake = ~(day_time >= sleep_start_time || day_time < mod(sleep_start_time + T_sleep, 24));
        end
        
        if is_awake
            Atot = mu_w + (Atot - mu_w)*exp(-dt/chi_w);
        else
            Atot = mu_s + (Atot - mu_s)*exp(-dt/chi_s);
        end
        
        R1b = 0.5*(Atot + R1tot + Kd1/(1-beta) - ...
            sqrt((Atot + R1tot + Kd1/(1-beta))^2 - 4*Atot*R1tot));
        
        dR1tot = (R1b - gamma*R1tot)/lambda;
        R1tot = R1tot + dR1tot*dt;
    end
    
    Atot_end = Atot;
    R1tot_end = R1tot;
end

function [Atot_end, R1tot_end] = simulate_one_day_mfactor(Atot_init, R1tot_init, T_sleep, sleep_start_time, ...
    mu_s, mu_w, chi_s, chi_w, Kd1, beta, gamma, lambda_fast, dt)
    
    t_day = 0:dt:24;
    Atot = Atot_init;
    R1tot = R1tot_init;
    
    for i = 1:length(t_day)-1
        day_time = mod(t_day(i), 24);
        
        if sleep_start_time < 0
            sleep_start_time = sleep_start_time + 24;
        end
        
        if sleep_start_time + T_sleep <= 24
            is_awake = ~(day_time >= sleep_start_time && day_time < sleep_start_time + T_sleep);
        else
            is_awake = ~(day_time >= sleep_start_time || day_time < mod(sleep_start_time + T_sleep, 24));
        end
        
        if is_awake
            Atot = mu_w + (Atot - mu_w)*exp(-dt/chi_w);
        else
            Atot = mu_s + (Atot - mu_s)*exp(-dt/chi_s);
        end
        
        R1b = 0.5*(Atot + R1tot + Kd1/(1-beta) - ...
            sqrt((Atot + R1tot + Kd1/(1-beta))^2 - 4*Atot*R1tot));
        
        dR1tot = (R1b - gamma*R1tot)/lambda_fast;
        R1tot = R1tot + dR1tot*dt;
    end
    
    Atot_end = Atot;
    R1tot_end = R1tot;
end

function [Atot_wake, R1tot_wake] = simulate_to_wake_helper(Atot_init, R1tot_init, T_sleep, sleep_start_time, ...
    mu_s, ~, chi_s, ~, Kd1, beta, gamma, lambda, dt)
    
    if sleep_start_time < 0
        sleep_start_time = sleep_start_time + 24;
    end
    
    wake_time = mod(sleep_start_time + T_sleep, 24);
    t_sim = 0:dt:T_sleep;
    
    Atot = Atot_init;
    R1tot = R1tot_init;
    
    for i = 1:length(t_sim)-1
        Atot = mu_s + (Atot - mu_s)*exp(-dt/chi_s);
        
        R1b = 0.5*(Atot + R1tot + Kd1/(1-beta) - ...
            sqrt((Atot + R1tot + Kd1/(1-beta))^2 - 4*Atot*R1tot));
        
        dR1tot = (R1b - gamma*R1tot)/lambda;
        R1tot = R1tot + dR1tot*dt;
    end
    
    Atot_wake = Atot;
    R1tot_wake = R1tot;
end