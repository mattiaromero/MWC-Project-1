clear; clc; close all;

%% Chapter 1 - Time Series

%% 2.1 Pre-processing 

%% Load the data of the raw signal 
x = load('CalibP1.txt', '-ascii', 'Frequency');

% Sampling frequency
dt = 1/4; 

% Duration of the signal in minutes
D = length(x); %0.25s 
D = D/(4*60); %min

%% Creating a vector containing the time-series of total water depth (h) 

% Converting raw data into pressure 
a1 = 19.97370; %Pa*mV^{-1}
b1 = -49.95959; %Pa

p = a1*x+b1; %Pa

% The pressure measured by the sensor is related to the height of the 
% water column above the sensor (ha) 
ro = 1025; %Kg*m^{3} 
g = 9.81; %m*s^{2} 

ha = p/(ro*g); %m

% Total depth of the water column (h)
hb = 1.45; %m 

h = ha+hb; %m

% Mean water depth at the sensor (hmean) 
hmean = mean(h);

% Creating a time vector in seconds corresponding to the time-series
t = (linspace(0,D*60,D*4*60))';

% Plot of the evolution of h as a function of time. 
figure; %CAREFUL FIX X AXIS to put in in Seconds

plot(t,h);
hold on;
plot(t,linspace(hmean,hmean,length(t)),'LineWidth',2);
title('Time evolution of the water depth h [m]');
xlabel('Time [s]','FontWeight','bold');
ylabel('Water depth [m]','FontWeight','bold');
legend('Water depth','Mean water depth');
grid on;
xlim([0 3600]);
savefig('Matlab1_i');

% It can be seen in the figure that the water level is increasing with time. 
% Can you explain what is happening? 
% The tide is rising

% Remove the tidal trend from the signal 
hdetrend = detrend(h);

% Plot of the evolution of h as a function of time. 
figure;
plot(t,hdetrend);
hold on;
plot(t,h-hmean);
title('Time evolution of the water depth h [m]');
xlabel('Time [s]','FontWeight','bold');
ylabel('Water level variation [m]','FontWeight','bold');
legend('Tidally averaged water level variation','h-h_{mean}');
grid on;
xlim([0 3600]);
savefig('Matlab1_ii');

%% Remark
% In this first part we neglected the depth-attenuation of pressure in
% our calibration process. Additional corrections can be included to take these
% effects into account.
%% 2.2 Wave statistics
%Loading the files LowTide.txt and highTide.txt
lowTide = load('lowTide.txt');
highTide = load('highTide.txt');
lowP1 = lowTide(:,1);
lowP3 = lowTide(:,2);
lowP6 = lowTide(:,5);

highP1 = highTide(:,1);
highP3 = highTide(:,2);
highP6 = highTide(:,5);

% Sampling frequency
dt = 1/2; 

% Duration of the signal in minutes
D = length(lowP1); 
D = D/(2*60); %min

%Creating the new time vector
tt = (linspace(0,D*60,D*2*60))'; %new t vector for this section

%Plot high and low tide waves at P1,P3 and P6
figure;
subplot(3,2,1);
plot(tt,lowP1);
ylabel('Water level variation [m]','FontWeight','bold');
title("Low tide at P1"); 
grid on;
xlim([0 D*60]);
ylim([-1.5 1.75]);

subplot(3,2,2);
plot(tt,highP1);
title("High tide at P1"); 
grid on;
xlim([0 D*60]);
ylim([-1.5 1.75]);

subplot(3,2,3);
plot(tt,lowP3);
ylabel('Water level variation [m]','FontWeight','bold');
title("Low tide at P3"); 
grid on;
xlim([0 D*60]);
ylim([-1.5 1.75]);

subplot(3,2,4);
plot(tt,highP3);
title("High tide at P3"); 
grid on;
xlim([0 D*60]);
ylim([-1.5 1.75]);

subplot(3,2,5);
plot(tt,lowP6);
xlabel('Time [s]','FontWeight','bold');
ylabel('Water level variation [m]','FontWeight','bold');
title("Low tide at P6"); 
grid on;
xlim([0 D*60]);
ylim([-1.5 1.75]);

subplot(3,2,6);
plot(tt,highP6);
xlabel('Time [s]','FontWeight','bold');
title("High tide at P6"); 
grid on;
xlim([0 D*60]);
ylim([-1.5 1.5]);
savefig('Matlab1_iii');