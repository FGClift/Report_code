% Parameters
mu_s = 596.4;
mu_w = 869.5;
chi_s = 4.2;
chi_w = 18.18;
Kd1 = 1;  % nM
Kd2 = 100;  % nM
gamma = 0.9677;
lambda = 291;  % hours
R2u = 300;  % nM (approximately R2,tot since Kd2 is large)
T = 4; % time asleep per day
t0 = 11; % baseline sleep ends at 11am (4am + 7h)
dt = 0.5;
time = t0:dt:(24*8)+t0;

% Initial condition, A_tot,w(0) given by (24)
initial_value_Atot_t0 = (mu_s*(1-exp(-T/chi_s)) + mu_w*exp(-T/chi_s)*(1-exp((T-24)/chi_w))) / ...
                        (1-exp((T-24)/chi_w-(T/chi_s)));


%% Figure 2B: Calculate Atot 

% Acute sleep restriction condition (red)
Atot_acute = zeros(size(time));
initial_value_Atot_w1= mu_w + (initial_value_Atot_t0 - mu_w)*exp(-t0/chi_w); % The simulation starts at T+7 hours of sleep=11
for i = 1:length(time)
    if time(i) <= 107  % 4 days of continuous wake
        Atot_acute(i) = mu_w + (initial_value_Atot_t0 - mu_w)*exp(-time(i)/chi_w);
    end
end

% Chronic sleep restriction condition (black)
Atot_chronic = zeros(size(time));

for i = 1:length(time)
    day_time = mod(time(i), 24);
    
    % Sleep period: 16:00 to 16+T (4pm to 8pm in terms of day_time)
    if day_time >= 16 && day_time < 16 + T
        % We're in a sleep period
        if i == 1
            % First time point
            initial_s = initial_value_Atot_w1;
        else
       
          initial_s = Atot_chronic(i-1);

        end
        
        % Time elapsed in current sleep period
        time_in_sleep = day_time - 16;
        
        % Calculate Atot using previous value and time step
        if i > 1
            % Here we use forward euler integration to calulate next step 
            % A(t+dt) ≈ A(t) + dt*(mu - A(t))/chi = mu + (A(t) - mu)*(1 - dt/chi)
            Atot_chronic(i) = mu_s + (Atot_chronic(i-1) - mu_s)*(1 - dt/chi_s);
        else
            Atot_chronic(i) = mu_s + (initial_s - mu_s)*exp(-time_in_sleep/chi_s);
        end
        
    else
        % We're in a wake period
        if i == 1
            % First time point (awake)
            Atot_chronic(i) = initial_value_Atot_w1;
        else
                initial_w = Atot_chronic(i-1);
          
            % Use Euler integration to ensure continuity
            Atot_chronic(i) = mu_w + (Atot_chronic(i-1) - mu_w)*(1 - dt/chi_w);
        end
    end
end

Atot_acute_end = Atot_acute(time == 107) * ones(size(time));

%% Figure 2A: Calculate Au
beta = R2u/(R2u + Kd2);
R1tot_approx = 600;  % Approximate for Au calculation

Au_acute = zeros(size(Atot_acute));
Au_chronic = zeros(size(Atot_chronic));

for i = 1:length(time)
    % For acute
    if Atot_acute(i) > 0
        % Equation 11 to calculate R1b
        A1b_temp = 0.5*(Atot_acute(i) + R1tot_approx + Kd1/(1-beta) - ...
            sqrt((Atot_acute(i) + R1tot_approx + Kd1/(1-beta))^2 - 4*Atot_acute(i)*R1tot_approx)); 
        % Au=(1-beta)*(Atot-A1b)
        Au_acute(i) = (1-beta)*(Atot_acute(i) - A1b_temp); % Can use a constant value of R1tot here since occuring on fast timescale where R1tot willl be approx constant
    end
    
    % For chronic
    if Atot_chronic(i) > 0
        A1b_temp = 0.5*(Atot_chronic(i) + R1tot_approx + Kd1/(1-beta) - ...
            sqrt((Atot_chronic(i) + R1tot_approx + Kd1/(1-beta))^2 - 4*Atot_chronic(i)*R1tot_approx));
        Au_chronic(i) = (1-beta)*(Atot_chronic(i) - A1b_temp);
    end
end

Au_acute_end = Au_acute(time-11 == 96) * ones(size(time));


%% Figure 2C & 2D: Calculate R1b and R1tot


% Acute sleep deprivation
R1tot_acute = zeros(size(time));
R1b_acute = zeros(size(time));
R1tot_acute(1) = 585.5;  % Initial value guess

for i = 1:length(time)
    % Calculate R1b from Equation 11
    R1b_acute(i) = 0.5*(Atot_acute(i) + R1tot_acute(i) + Kd1/(1-beta) - ...
        sqrt((Atot_acute(i) + R1tot_acute(i) + Kd1/(1-beta))^2 - 4*Atot_acute(i)*R1tot_acute(i)));
    
    % Update R1tot using Equation 10
    dR1tot = (R1b_acute(i) - gamma*R1tot_acute(i))/lambda;
    R1tot_acute(i+1) = R1tot_acute(i) + dR1tot*dt;
   R1tot_acute = R1tot_acute(1:length(time));

end

% Chronic sleep restriction
R1tot_chronic = zeros(size(time));
R1b_chronic = zeros(size(time));
R1tot_chronic(1) = 585.5;

for i = 1:length(time)
    R1b_chronic(i) = 0.5*(Atot_chronic(i) + R1tot_chronic(i) + Kd1/(1-beta) - ...
        sqrt((Atot_chronic(i) + R1tot_chronic(i) + Kd1/(1-beta))^2 - 4*Atot_chronic(i)*R1tot_chronic(i)));
    
    dR1tot = (R1b_chronic(i) - gamma*R1tot_chronic(i))/lambda;
    R1tot_chronic(i+1) = R1tot_chronic(i) + dR1tot*dt;
R1tot_chronic = R1tot_chronic(1:length(time));
  
end

R1b_acute_end = R1b_acute(time -11 == 96) * ones(size(time));
R1tot_acute_end = R1tot_acute(time -11 == 96) * ones(size(time));

%% Figure of A_2B
R2b_acute = zeros(size(time));
R2b_chronic = zeros(size(time));
for i=1:length(time)
    R2b_acute(i)=Au_acute(i)*(beta/(1-beta));
    R2b_chronic(i)=Au_chronic(i)*(beta/(1-beta));
end 
%% Sleep drive calculation

%Parameters
a=3.25;
omega=2*pi/24;
phi=7.95;
p_max=60;
D_mid=583.2;
D_s=5.872;
D_wake= 555.4;
D_sleep=572.7;

% D= R1b + a*cos(omega*(t-phi))

 D_acute=zeros(size(time));
 D_chronic=zeros(size(time));
 p_acute=zeros(size(time));
p_chronic=zeros(size(time));
for i=1:length(time)
    D_acute(i)= R1b_acute(i) + a*cos(omega*(time(i)-phi));
    D_chronic(i)= R1b_chronic(i) + a*cos(omega*(time(i)-phi));
end 


%% Performance calculation
for i=1:length(time)
    p_acute(i)=p_max/(1+exp((D_mid-D_acute(i))/D_s));
    p_chronic(i)=p_max/(1+exp((D_mid-D_chronic(i))/D_s));
end 





%% Plotting

% Create separate figures for each plot
% Figure for Unbound Adenosine Concentration
figure('Units','normalized','Position',[0.2 0.1 0.6 0.8]);
subplot(4,1,1)
plot((time(time >= 15 & time -11 <= 96)-11)/24, Au_acute(time >= 15 & time -11 <= 96), 'r', 'LineWidth', 2);
hold on 
plot((time(time>=11)-11)/24, Au_chronic(time>=11), 'k', 'LineWidth', 2);
hold on
plot((time-11)/24, Au_acute_end, 'r --', 'LineWidth', 1);
ylim([0, 80]);
yticks(0:20:80)
xlim([0 8]);
% Hide x-axis tick labels but keep tick marks visible for vertical grid lines
xticks(0:1:8);            
ax = gca;
ax.XTickLabel = [];      % remove tick labels
ax.XGrid = 'on';         % enable vertical grid lines
ax.GridLineStyle = '-';
ax.GridColor = [0.15 0.15 0.15];
ax.GridAlpha = 0.15;
ylabel('A_{u} (nM)', 'FontSize', 40);
grid on;
set(gca, 'FontSize', 40); % Set font size for x and y ticks
text(0.1, 71, 'A', 'FontSize', 40, 'FontWeight', 'bold', 'Color', 'k');
% Figure for Total Adenosine Concentration

subplot(4,1,2)
plot((time(time >= 15 & time <= 107)-11)/24, Atot_acute(time >= 15 & time <= 107), 'r', 'LineWidth', 2);
hold on
plot((time(time>=11)-11)/24, Atot_chronic(time>=11), 'k', 'LineWidth', 2);
hold on
plot((time-11)/24, Atot_acute_end, 'r --', 'LineWidth', 1);
ylim([600, 900]);
xlim([0 8]);
ylabel('A_{tot} (nM)', 'FontSize', 40);
xticks(0:1:8);            
ax = gca;
ax.XTickLabel = [];      % remove tick labels
ax.XGrid = 'on';         % enable vertical grid lines
ax.GridLineStyle = '-';
ax.GridColor = [0.15 0.15 0.15];
ax.GridAlpha = 0.15;
grid on;
set(gca, 'FontSize', 40); % Set font size for x and y ticks
text(0.1, 866, 'B', 'FontSize', 40, 'FontWeight', 'bold', 'Color', 'k');

subplot(4,1,3)
plot((time(time >= 15 & time <= 107)-11)/24, R1b_acute(time >= 15 & time <= 107), 'r', 'LineWidth', 2);
hold on
plot((time(time>=11)-11)/24, R1b_chronic(time>=11), 'k', 'LineWidth', 2);
hold on 
plot((time-11)/24, R1b_acute_end, 'r --', 'LineWidth', 1);
ylim([560, 590]);
xlim([0 8]);
ylabel('R_{1,b} (nM)', 'FontSize', 30);
xticks(0:1:8);            
ax = gca;
ax.XTickLabel = [];      % remove tick labels
ax.XGrid = 'on';         % enable vertical grid lines
ax.GridLineStyle = '-';
ax.GridColor = [0.15 0.15 0.15];
ax.GridAlpha = 0.15;
grid on;
set(gca, 'FontSize', 40); % Set font size for x and y ticks
text(0.1, 586, 'C', 'FontSize', 40, 'FontWeight', 'bold', 'Color', 'k');

% Figure for Total Adenosine Receptor Concentration
subplot(4,1,4)
plot((time(time >= 15 & time <= 107)-11)/24, R1tot_acute(time >= 15 & time <= 107), 'r', 'LineWidth', 2);
hold on
plot((time(time>=11)-11)/24, R1tot_chronic(time>=11), 'k', 'LineWidth', 2);
hold on 
plot((time-11)/24, R1tot_acute_end, 'r --', 'LineWidth', 1);
ylim([585, 590]);
xlim([0 8]);
xlabel('Time (days)', 'FontSize', 40);
ylabel('R_{1,tot} (nM)', 'FontSize', 40);
xticks(0:1:8);            
ax = gca;
ax.XGrid = 'on';         % enable vertical grid lines
ax.GridLineStyle = '-';
ax.GridColor = [0.15 0.15 0.15];
ax.GridAlpha = 0.15;
legend('0 hours time in bed ','4 hours time in bed','Location','southeast','Fontsize',30)
grid on;
set(gca, 'FontSize', 40); % Set font size for x and y ticks
text(0.1, 589.4, 'D', 'FontSize', 40, 'FontWeight', 'bold', 'Color', 'k');

