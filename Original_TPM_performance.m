%% Parameters 
% Time constants for process S
Tr = 18.2; % hours (rising time constant during wake)
Td = 4.2;  % hours (declining time constant during sleep)


% Circadian process C parameters
A = 0.12;   % amplitude
t0 = 8.6;  % phase (hours)
yk = [0.97, 0.22, 0.07, 0.03, 0.001]; % harmonic weights 


%% Simulation Setup 
dt = 0.1; % time step in hours (6 minutes)
total_hours = 24*14; % simulate 16 days
wake_up=7; %7 wake up as specified in  paper
time = wake_up:dt:total_hours+wake_up; 


conditions = struct();
% Condition 1: Total Sleep Deprivation (0h TIB for 3 days)
conditions(1).name = '0h TIB';
conditions(1).TIB = 0;
conditions(1).n_restriction_days = 3;
conditions(1).color = 'k';

% Condition 2: 4h TIB for 14 days
conditions(2).name = '4h TIB';
conditions(2).TIB = 4;
conditions(2).n_restriction_days = 14;
conditions(2).color = [0.8, 0, 0]; 

% Condition 3: 6h TIB for 14 days
conditions(3).name = '6h TIB';
conditions(3).TIB = 6;
conditions(3).n_restriction_days = 14;
conditions(3).color = [0.95, 0.5, 0];
% 8h TIB condition
conditions(4).name = '8h TIB';
conditions(4).TIB = 8;
conditions(4).n_restriction_days = 14;
conditions(4).color = [0.1, 0.4, 0];

%% Loop through each condition

for cond_idx = 1:length(conditions)
    cond = conditions(cond_idx);

    
    % Initialize Variables
    wake_status = zeros(size(time));
    S = zeros(size(time)); % Homeostatic process
    U = zeros(size(time)); % Process U (slow modulating process)
    L = zeros(size(time)); % Lower asymptote for U
    C = zeros(size(time)); % Circadian process
    PVT_lapses = zeros(size(time)); % Predicted PVT lapses

    for i=1:length(time) % Define sleep/wake schedule for each time point
        hour_of_day = mod(time(i), 24);
        if cond.TIB == 0
            % Total sleep deprivation - always awake
            wake_status(i) = 1;
        else 
            if cond.TIB == 4
                sleep_start = 3; %3am
                sleep_end = wake_up;
            elseif cond.TIB==6
                sleep_start = 1; %1am
                sleep_end = wake_up;
            elseif cond.TIB==8
                sleep_start = 23; %11pm
                sleep_end = wake_up;
            end   

            % Determine if awake
            if sleep_start < sleep_end
                % Sleep doesn't cross midnight
                if hour_of_day >= sleep_end 
                    wake_status(i) = 1;
                elseif hour_of_day >= sleep_start && hour_of_day < sleep_end
                    wake_status(i) = 0;
                else
                    wake_status(i) = 1;
                end
            elseif sleep_start > sleep_end % Sleep crosses midnight (e.g., 23:00 to 07:00)
                if hour_of_day >= sleep_end && hour_of_day < sleep_start
                    wake_status(i) = 1;
                else
                    wake_status(i) = 0;
                end
            end
        end 
    end 



    %% Circadian Process C
    for i = 1:length(time)
        t = time(i);
        C(i) = A * sum(yk .* sin(2*pi*(1:5)*(t - t0)/24));
    end
        %% Initial Conditions (assuming 8h sleep baseline)
    % After 8h sleep and 16h wake at baseline
    U(1) = 1 ;
 S(1) = U(1) - (1 - exp(-8/Td)) / (1 - exp(-(16/Tr + 8/Td)));
    L(1) = U(1) - 1;

    %% Simulate Homeostatic Process S and Process U
    for i = 2:length(time)
        U(i)=1;
        L(i)=0;
        if wake_status(i) == 1 % During wake
         
     
            S(i) = U(i) + (S(i-1) - U(i-1)) * exp(-dt/Tr);
 
        else % During sleep
        
            S(i) = L(i) + (S(i-1) - L(i-1)) * exp(-dt/Td);
        end
    end 

    %% Predict PVT Performance
    for i = 1:length(time)
        combined= S(i) - C(i);
         a = 13; 
         b=4.49-a*(S(1)-C(1));
        % Linear regression to predict lapses
        PVT_lapses(i) = a * combined+b ;
        % Ensure non-negative lapses
        PVT_lapses(i) = max(0, PVT_lapses(i));
    end
    
    % Plot only the 0h TIB condition for 88h
    if cond.TIB == 0
        plot((time(1:889)-wake_up)/24, PVT_lapses(1:889),'Color', cond.color, 'LineWidth', 2.5, 'DisplayName', cond.name);
    else
        plot((time-wake_up)/24, PVT_lapses,'color', cond.color, 'LineWidth', 2.5, 'DisplayName', cond.name);
    end
    hold on

end
xlabel('Time (days)','FontSize',40)
ylabel({'Performance','(PVT lapses)'}, 'FontSize', 40, 'HorizontalAlignment', 'center','position',[-1, 15]);
set(gca, 'FontSize', 40); 
annotation('arrow', [0.04, 0.04], [0.81, 0.88]); 
annotation('arrow', [0.04, 0.04], [0.31, 0.24]); 
text(-1.6, 28, 'Worse', 'FontSize', 30, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'Rotation', 0)
text(-1.6, 2, 'Better', 'FontSize', 30, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'Rotation', 0)
legend('Location','northeast','FontSize',24);
text(0.2,28,'A','FontSize',35,'FontWeight','bold')
grid on;
xlim([0 14]);
ylim([0 30]);
xticks(0:1:14)
yticks(0:5:30)
figure;
for cond_idx = 1:length(conditions)
    cond = conditions(cond_idx);

    % Plot experimental data 
    M = load('vandongenbu.mat');
    flds = fieldnames(M);
    data = M.(flds{cond_idx});
    t   = data(:,1);
    val = data(:,2);

    t_plot = t / 24;

markers = {'d', 'o', 's', '^'}; % diamond, circle, square, triangle
marker = markers{cond_idx};
scatter(t_plot(t_plot<=14), val(t_plot<=14), 50, 'Marker', marker, ...
    'MarkerEdgeColor', cond.color, 'MarkerFaceColor', cond.color, 'LineWidth', 0.8);

    hold on
    % Draw lines connecting the points in sorted time order
    [t_sorted, idx] = sort(t_plot);
    val_sorted = val(idx);
     if cond_idx == 1
        plot(t_sorted, val_sorted, '-', 'Color', cond.color, 'LineWidth', 1.5);
    else
        plot(t_sorted(t_sorted<=14), val_sorted(t_sorted<=14), '-', 'Color', cond.color, 'LineWidth', 1.5);
    end
    
    hold on
end
xlabel('Time (days)','FontSize',40)
ylabel({'Performance','(PVT lapses)'}, 'FontSize', 40, 'HorizontalAlignment', 'center','position',[-1, 15]);
set(gca, 'FontSize', 40); 
annotation('arrow', [0.04, 0.04], [0.81, 0.88]); 
annotation('arrow', [0.04, 0.04], [0.31, 0.24]); 
text(-1.6, 28, 'Worse', 'FontSize', 30, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'Rotation', 0)
text(-1.6, 2, 'Better', 'FontSize', 30, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'Rotation', 0)
text(0.2,28,'B','FontSize',35,'FontWeight','bold')
grid on;
box on;
xlim([0 14]);
ylim([0 30]);
xticks(0:1:14)
yticks(0:5:30)
    