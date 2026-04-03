% Parameters
mu_s = 596.4;
mu_w = 869.5;
chi_s = 4.2;
chi_w = 18.18;
Kd1 = 1; % nM
Kd2 = 100; % nM
gamma = 0.9677;
lambda = 291; % hours
R2tot= 300;
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

% Caffeine parameters
Kcaf2 = 50000;            % nM, caffeine dissociation constant at A2A
k_cl  = 2/7*log(2);      % per hour, clearance (half-life ~3.5 hours)

% Dosing: two doses on day 2 and day 3
caffeine_dose = 50000; % nM= 100uM which is about 1/5th of the content from a cup of coffee                       % nM equivalent
dose_times    = [t0 + 24, t0 + 48];        % hours


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

%% Caffeine concentration C_u (instantaneous absorption, analytical superposition)
Cu = zeros(size(t));
for j = 1:length(dose_times)
    td = dose_times(j);
    mask = t >= td;
    Cu(mask) = Cu(mask) + caffeine_dose * exp(-k_cl * (t(mask) - td));
end

%% Simulate for three different TIB conditions: 4h, 6h, 8h
T_conditions = [0, 4, 6, 8]; % hours of time in bed
colors = {'k',[0.8, 0, 0], [0.95, 0.65, 0], [0.1, 0.4, 0]};  


p_chronic_all = cell(1, 4);
R2b_chronic_all = cell(1,4);
R2c_chronic_all = cell(1,4);
wake_all = cell(1,4);

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

    % Initialize arrays for this condition
    Atot_chronic = zeros(size(t));
    R1tot_chronic = zeros(size(t));
    R1b_chronic = zeros(size(t));
    D_chronic = zeros(size(t));
    p_chronic = zeros(size(t));
    wake = zeros(size(t));
    actual_wake = zeros(size(t));
    R2b_acute    = zeros(size(t));
R2c_acute    = zeros(size(t));
R2b_chronic   = zeros(size(t));
R2c_chronic   = zeros(size(t));
Au_acute     = zeros(size(t));
Au_chronic    = zeros(size(t));
R2u_chronic  = zeros(size(t));
    
    R1tot_chronic(1) = initial_R1tot;
    Atot_chronic(1) = initial_Atot;
    
    for i = 1:length(t)-1
        day_time = mod(t(i)+t0, 24); % Actual clock time
        
        if T_current == 0
            wake(i) = 1;
        elseif sleep_start_chronic < wake_time
            wake(i) = ~(day_time >= sleep_start_chronic && day_time < wake_time);
        else
            wake(i) = ~(day_time >= sleep_start_chronic || day_time < wake_time);
        end
        
   
        
        % Calculate Atot 
        if wake(i) == 0
            % Sleep
            Atot_chronic(i+1) = mu_s + (Atot_chronic(i) - mu_s)*exp(-dt/chi_s);
        else
            % Wake
            Atot_chronic(i+1) = mu_w + (Atot_chronic(i) - mu_w)*exp(-dt/chi_w);
        end
        
        % Calculate R1b
        R1b_chronic(i) = 0.5*(Atot_chronic(i) + R1tot_chronic(i) + Kd1/(1-beta) - ...
            sqrt((Atot_chronic(i) + R1tot_chronic(i) + Kd1/(1-beta))^2 - 4*Atot_chronic(i)*R1tot_chronic(i)));
        
        % Update R1tot
        dR1tot = (R1b_chronic(i) - gamma*R1tot_chronic(i))/lambda;
        R1tot_chronic(i+1) = R1tot_chronic(i) + dR1tot*dt;



 % Caffeine 
    Cu_val = Cu(i);
        phi_val = 1 + Cu_val / Kcaf2;

    Q_chronic  = Atot_chronic(i) - R1b_chronic(i);
    b_c = phi_val*Kd2 + R2tot - Q_chronic;
    c_c = -phi_val * Kd2 * Q_chronic;
    Au_chronic(i) = 0.5*(-b_c + sqrt(b_c^2 - 4*c_c));

    % R2b and R2c from competitive quasi-static binding
    %   R2b = (Au/Kd2)  / (1 + Au/Kd2 + Cu/Kcaf2) * R2tot
    %   R2c = (Cu/Kcaf2)/ (1 + Au/Kd2 + Cu/Kcaf2) * R2tot

    denom_chronic = 1 + Au_chronic(i)/Kd2 + Cu_val/Kcaf2;


    R2b_chronic(i) = (Au_chronic(i)/Kd2)  / denom_chronic * R2tot;
    R2c_chronic(i) = (Cu_val/Kcaf2)       / denom_chronic * R2tot;
R2u_chronic(i)= R2tot/ denom_chronic;
        
        % Calculate sleep drive D
        D_chronic(i) = R1b_chronic(i)-0.1*R2c_chronic(i)+ a*cos(omega*(t(i)+t0-phi));
        
        % Calculate performance
        p_chronic(i) = p_max/(1+exp((D_mid-D_chronic(i))/D_s));
    end
    
    % Final point
    i = length(t);
    day_time = mod(t(i)+t0, 24); % Actual clock time
     if T_current == 0
            wake(i) = 1;
        elseif sleep_start_chronic < wake_time
            wake(i) = ~(day_time >= sleep_start_chronic && day_time < wake_time);
        else
            wake(i) = ~(day_time >= sleep_start_chronic || day_time < wake_time);
     end
 

    
    
    % Store results
    p_chronic_all{cond} = p_chronic;
    R2b_chronic_all{cond} = R2b_chronic;
R2c_chronic_all{cond} = R2c_chronic;
R2u_chronic_all{cond} = R2u_chronic;
wake_all{cond} = wake;
  
end

% Create a vector (same size as t(t>=t0)) filled with the value of R2b_chronic_all{1} at t==88
idx88 = find(abs(t-107) < eps, 1);
R2b_acute_end = R2b_chronic_all{1}(idx88) * ones(size(t(t>=t0)));
R2c_acute_end =  R2c_chronic_all{1}(idx88) * ones(size(t(t>=t0)));

%% Plotting (downsample for plotting)
dt_plot = 2;
t_plot_indices = 1:round(dt_plot/dt):length(t)-1;
t_plot_days = t(t_plot_indices)/24;
sleep_inertia_hours = 2;


annotation_labels = {'A','B','C','D'};

% Performance figures - combined into one figure with 4 subplots
%figure('Units','normalized','Position',[0 0 1 1]);
pad = 0.06;           % smaller padding -> larger subplot area
w = (1 - 3*pad) / 2;  % adjusted to match padding layout
h = (1 - 5*pad) / 2; 

for cond = 1:4

     if cond == 1; pos = [pad, pad + h + pad, w, h];
    elseif cond == 2; pos = [pad + w + pad, pad + h + pad, w, h];
    elseif cond == 3; pos = [pad, pad, w, h];
    else; pos = [pad + w + pad, pad, w, h];
     end

    subplot('position',pos);
    
    wake_vec = wake_all{cond};
    % Find transitions: 1 to 0 (sleep start) and 0 to 1 (sleep end)
    sleep_starts = find(diff(wake_vec) == -1);
    sleep_ends = find(diff(wake_vec) == 1);
    
    % Handle edge cases (if simulation starts/ends in sleep)
    if wake_vec(1) == 0, sleep_starts = [1, sleep_starts]; end
    if wake_vec(end) == 0, sleep_ends = [sleep_ends, length(wake_vec)]; end
    
    for s = 1:length(sleep_starts)
        t_start = t(sleep_starts(s))/24;
        t_end = t(sleep_ends(s))/24;
        patch([t_start t_end t_end t_start], [0 0 40 40], [0.8 0.8 0.8], ...
              'EdgeColor', 'none', 'HandleVisibility', 'off');
        hold on
    end

    % Only plot after sleep inertia hours and when awake
    valid_idx = false(size(t_plot_indices));
    for i = 1:length(t_plot_indices)
        time_since_wake = mod(t(t_plot_indices(i)), 24);
        valid_idx(i) = wake_all{cond}(t_plot_indices(i)) == 1 && ...
                       time_since_wake >= sleep_inertia_hours;
    end

    plot((t_plot_days(valid_idx)), p_chronic_all{cond}(t_plot_indices(valid_idx)), '.', 'Color', colors{cond}, 'MarkerSize', 20)
    hold on
    for j = 1:length(dose_times)
        xline(dose_times(j)/24, 'b--', 'LineWidth', 4, 'Alpha', 0.8);
    end

    % Formatting
    set(gca, 'Layer', 'top');
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
      
    else
        xlim([0 13]); xticks(0:1:13);
    end

    x_ann = pos(1) + 0.01*pos(3);
    y_ann = pos(2) + pos(4) - 0.15*pos(4);
    annotation('textbox', [x_ann, y_ann, 0.05, 0.05], 'String', annotation_labels{cond}, ...
        'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 45);

    hold off
end

% Dynamics figures
figure;
subplot(2,1,1)
plot((t(t<=96+t0 & t>=t0)-t0)/24, R2b_chronic_all{1}(t<=96+t0 & t>=t0), 'r', 'LineWidth', 3)
hold on
plot((t(t>=t0)-t0)/24,R2b_chronic_all{2}(t>=t0), 'k', 'LineWidth', 3)
hold on
plot((t(t>=t0)-t0)/24,R2b_acute_end, 'r --', 'LineWidth', 2)
hold on 
for j = 1:length(dose_times)
    xline((dose_times(j)-t0)/24,'b--', 'LineWidth', 3, 'Alpha', 0.8);

end
set(gca, 'FontSize', 35); % Set font size for x and y ticks
ylabel('R_{2,b} (nM)', 'FontSize', 35);
xlim([0 8])
ylim([50 200])
yticks(50:50:200)
xticks(0:1:8)
text(0.1, 185, 'A', 'FontSize', 35, 'FontWeight', 'bold', 'Color', 'k');

subplot(2,1,2)
plot((t(t<=96+t0 & t>=t0)-t0)/24, R2c_chronic_all{1}(t<=96+t0 & t>=t0),'r','LineWidth',3)
hold on
plot((t(t>=t0)-t0)/24,R2c_chronic_all{2}(t>=t0), 'k', 'LineWidth', 3)
hold on
for j = 1:length(dose_times)
 xline((dose_times(j)-t0)/24,'b--', 'LineWidth', 3, 'Alpha', 0.8);
end
xlabel('Time (days)', 'FontSize', 35);
ylabel('R_{2,c} (nM)', 'FontSize', 35);
set(gca, 'FontSize', 35); % Set font size for x and y ticks
xlabel('Time (days)', 'FontSize', 35);
ylabel('R_{2,c} (nM)', 'FontSize', 35);
xlim([0 8])
xticks(0:1:8)
ylim([0 150])
yticks(0:50:150)
legend('0 hours time in bed ','4 hours time in bed','Location','northeast','Fontsize',22)
text(0.1, 135, 'B', 'FontSize', 35, 'FontWeight', 'bold', 'Color', 'k');







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