%parameters
dt= 0.001; %time step
T=48; %Total time (hours)
t= 0:dt:T;


%Circadian process equation

alpha= -0.6;%alpha is a phase shift that specifies the distance from the circadian max.
%alpha chosen so switch from sleep to wake occurs at appropriate time
omega= 2*pi/24;
a= -0.22; %amplitude of the sine wave
C=a*(0.97*sin(omega*(t-alpha))+0.22*sin(2*omega*(t-alpha))+0.07*sin(3*omega*(t-alpha))+0.03*sin(4*omega*(t-alpha))+0.001*sin(5*omega*(t-alpha)));


%lower threshold parameters
H0minus = 0.17;
H0plus= 0.67;
Hminus = H0minus*ones(size(t));
Hplus= H0plus*ones(size(t));



%Homeostatic process (H) parameters
tau_s= 4.2; %sleep time constant-leads to faster dissipation of sleep pressure
tau_w= 18.2; %awake time constant- slower build up of sleep pressure 
U_t =1; %U_t is the upper asymptote i.e the value H would reach if no switch to sleep occured
H0= 0.6;





% initialize H before the loop
H= zeros(size(t));
H(1)=H0;
wake = 0;

I=zeros(size(t));
I(1)=H(1)-C(1);

for i=2:length(t)

    if wake==1 
        Hnew = U_t+(H(i-1)-U_t)*exp(-dt/tau_w); %when awake update Hnew with wake equation
    elseif wake==0
        Hnew = H(i-1) * exp(-dt/tau_s); %when asleep update Hnew with sleep equation 
    end

    if Hplus(i)>Hnew-C(i) %S-C is below upper limit
        if Hminus(i)<Hnew-C(i)
            H(i) = Hnew; %if S-C is above the lower limit then S-C is same as one before 
            I(i)=H(i)-C(i);
        else
            wake = 1; %if its below then you must need to switch from alseep to awake
            H(i)=U_t+(H(i-1)-U_t)*exp(-dt/tau_w); 
            I(i)=H(i)-C(i);
        end
    else
        wake = 0; 
        H(i)=H(i-1) * exp(-dt/tau_s);
        I(i) = H(i)-C(i); %otherwise only sleep when already asleep 
    end

end 

%defining shaded region
sleep=(zeros(size(t)));
sleep(1)= 1;
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

for i = 1:length(t)
x_rect1([1, 4])= 0;
    if sleep(i+1) == 0 && sleep(i) == 1 % First time waking up
        x_rect1([2,3]) = t(i);
        if length(x_rect1) == 4 %Only keep pairs of start and end times
            break; % Exit loop if we have two pairs
        end
    end 
end 


for i=9000:length(t)
if sleep(i+1) == 1 && sleep(i) == 0 % First time entering sleep
        x_rect2([1,4]) = t(i);
    elseif sleep(i+1) == 0 && sleep(i) == 1 % First time waking up
        x_rect2([2,3]) = t(i);
        if length(x_rect2) == 4 %Only keep pairs of start and end times
            break; % Exit loop if we have two pairs
        end
    end
end

y_rect1= [1,1,0,0];
y_rect2= [1,1,0,0];

%Plotting results 
close all;

plot(t, Hplus, 'black');
hold on
plot(t, Hminus, 'black')
hold on
plot(t,I,linewidth=2)
ylim([0,1])
yticks('')
xlim([0,48])
xticks(0:12:48)
xlabel('Time (hours)','FontSize',30)
ylabel('Performance')
set(gca, 'FontSize', 50); 
hold on
fill(x_rect1,y_rect1,'k','FaceAlpha',0.1)
fill(x_rect2,y_rect2,'k','FaceAlpha',0.1)
legend('','','S(t)-C(t)','sleep','','Fontsize',40,'Location','Northeast')
hold off

