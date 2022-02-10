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
t = (linspace(0,D*4*60-1,D*4*60))';

% Plot of the evolution of h as a function of time. 
figure;
plot(t,h);
hold on;
plot(t,linspace(hmean,hmean,length(t)),'LineWidth',2);
title('Time evolution of the water depth h [m]');
xlabel('Time [s]','FontWeight','bold');
ylabel('Water depth [m]','FontWeight','bold');
legend('Water depth','Mean water depth');
grid on;
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
ylabel('Water depth [m]','FontWeight','bold');
legend('Tidally averaged water depth','h-h_{mean}');
grid on;
savefig('Matlab1_ii');

%% Remark
% In this first part we neglected the depth-attenuation of pressure in
% our calibration process. Additional corrections can be included to take these
% effects into account.

