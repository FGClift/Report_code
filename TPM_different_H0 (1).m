
%parameters
dt= 0.001; %time step
T=80; %Total time (hours)
t= 0:dt:T;


%Circadian process equation

alpha= 2.07;%alpha is a phase shift that specifies the distance from the circadian max.
%phi chosen so switch from sleep to wake occurs at appropriate time
omega= 2*pi/24;
C=0.97*sin(omega*(t-alpha))+0.22*sin(2*omega*(t-alpha))+0.07*sin(3*omega*(t-alpha))+0.03*sin(4*omega*(t-alpha))+0.001*sin(5*omega*(t-alpha));


%lower threshold parameters
H0minus = 0.09; % Varying this parameter alters lower threshold value
H0plus= 0.75;
a= 0.07;
Hminus = H0minus + a*C;
Hplus= H0plus+a*C;



%Homeostatic process (H) parameters
tau_s= 4.2; %sleep time constant-leads to faster dissipation of sleep pressure
tau_w= 18.2; %awake time constant- slower build up of sleep pressure 
U_t =1; %mu is the upper asymptote i.e the value H would reach if no switch to sleep occured
H0= 0.5;




% initialize H before the loop
H= zeros(size(t));
H(1)=H0;
wake = 1;

for i=2:length(t)

    if wake==1 
        Hnew = U_t+(H(i-1)-U_t)*exp(-dt/tau_w); %when awake update Hnew with wake equation
    elseif wake==0
        Hnew = H(i-1) * exp(-dt/tau_s); %when asleep update Hnew with sleep equation 
   
    end

    if Hplus(i)>Hnew %H is below upper limit
        if Hminus(i)<Hnew
            H(i) = Hnew; %if H is above the lower limit then H assigned to either awake or sleep equation
        else
            wake = 1; %if its below then you must need to switch from alseep to awake
            H(i) = U_t+(H(i-1)-U_t)*exp(-dt/tau_w); %awake
        end
    else
        wake = 0; 
        H(i) = H(i-1) * exp(-dt/tau_s); %otherwise only sleep when already asleep 
    end
end 


%defining shaded region
sleep=(zeros(size(t)));
sleep(1)= 0;
for i=2:length(t)
    if H(i) == H(i-1) * exp(-dt/tau_s) % if asleep
        sleep(i) = 1; % mark as asleep
    else
        sleep(i) = 0; % mark as awake
    end
end

% Define x coordinates for sleep periods
x_rect1 = [];
x_rect2 = [];
x_rect3 = [];% Initialize x coordinates for sleep periods
for i = 1:length(t)
    if sleep(i+1) == 1 && sleep(i) == 0 % First time entering sleep
        x_rect1([1,4]) = t(i);
    elseif sleep(i+1) == 0 && sleep(i) == 1 % First time waking up
        x_rect1([2,3]) = t(i);
        if length(x_rect1) == 4 %Only keep pairs of start and end times
            break; % Exit loop if we have two pairs
        end
    end 
end 
for i=14860:length (t)
     if sleep(i+1) == 1 && sleep(i) == 0 % First time entering sleep
        x_rect2([1,4]) = t(i);
    elseif sleep(i+1) == 0 && sleep(i) == 1 % First time waking up
        x_rect2([2,3]) = t(i);
        if length(x_rect2) == 4 %Only keep pairs of start and end times
            break; % Exit loop if we have two pairs
        end
    end
end

for i=51890:length (t)-1
     if sleep(i+1) == 1 && sleep(i) == 0 % First time entering sleep
        x_rect3([1,4]) = t(i);
     elseif sleep(i+1) == 0 && sleep(i) == 1 % First time waking up
        x_rect3([2,3]) = t(i);
         if length(x_rect3) == 4 %Only keep pairs of start and end times
            break; % Exit loop if we have two pairs
        end
     end 
end 
 

y_rect1= [1,1,0,0];
y_rect2= [1,1,0,0];
y_rect3= [1,1,0,0];


%Plotting results 
close all;

plot(t, Hplus, 'black');
hold on
plot(t, Hminus, 'black')
hold on
plot(t,H,linewidth=3)
xlabel('Time (hours)','FontSize',34)
ylabel  ('Process S','FontSize',34)
yticks([0 0.5 1]);
set(gca, 'FontSize', 34)
hold on
fill(x_rect1,y_rect1,'k','FaceAlpha',0.1)
fill(x_rect2,y_rect2,'k','FaceAlpha',0.1)
fill(x_rect3,y_rect3,'k','FaceAlpha',0.1)
hold off