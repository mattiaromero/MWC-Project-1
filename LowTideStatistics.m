% Wave-by-wave analysis of the low tide data, Egmond Coast 3d dataset
% The objective of the script is to compute wave statistics at different
% locations along a cross-shore transect and examine their cross-shore
% evolution


%% -------------------------------------
%           Initialisation
% -------------------------------------
clear all
close all

% Load data
data = load('lowTide.txt');
prof = load("prof1018.txt");

% Constants
g = 9.81;                    % acceleration of gravity (m/s^2)
rho = 1025;                  % water density (kg/m^3)
fs = 2;                      % sampling frequency (Hz)
Npos = 5;                    % number of cross-shore positions considered

% Initialisation vectors 
% These vectors will be used to store the wave statistics at each position
Hrms_tot = zeros(Npos,1);  % root mean square height (m)
H13_tot  = zeros(Npos,1);  % significant wave height (m)
Hm_tot   = zeros(Npos,1);  % mean wave height (m)

%% --------------------------------------
%     Computation of wave statistics
% --------------------------------------

for i=1:Npos  % loop on the positions
     %Frecuency is equal to 2 because the sampling ratio has aperiod of 0.5sec
    wave = zero_crossing(data(:,i),2);
    Hrms_tot(i) = rms_height(wave(:,1));
    H13_tot(i) = significant_height(wave(:,1));
    Hm_tot(i) = mean(wave(:,1)); %I DONT KNOW IF THIS IS CORRECT (Idk if by saying the water height they mean the values of the lowtide.txt)!!!!

end

% --------------------------------------
%                  Output
% --------------------------------------

% visualisation of outputs
% ?
% ?
positions = [4478, 4765, 4790, 4814, 4835] %Positions where the sensors are located
subplot(2,1,1)
plot(positions, Hrms_tot,"or")
title("Mean of RMS, wave height and H 1/3 ")

hold on
plot(positions, H13_tot,"*b")
plot(positions, Hm_tot,"*g") 
xlim([4300,5000])
legend(["Hrms","H1/3","Height"])
xlabel("Position (m)")
ylabel("height (m)")

subplot(2,1,2)
plot(prof(:,1),prof(:,2),"black")
xlim([4300,5000])
title("Depth profile")

xlabel("Position (m)")
ylabel("Elevation (m)")