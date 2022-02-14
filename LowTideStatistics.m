% Wave-by-wave analysis of the low tide data, Egmond Coast 3d dataset
% The objective of the script is to compute wave statistics at different
% locations along a cross-shore transect and examine their cross-shore
% evolution


% -------------------------------------
%            Initialisation
% -------------------------------------
clear all
close all

% load data
data = load('lowTide.txt');

% constants
g = 9.81;                    % acceleration of gravity (m/s^2)
rho = 1025;                  % water density (kg/m^3)
fs = 2;                      % sampling frequency (Hz)
Npos = 5;                    % number of cross-shore positions considered
% More? 

% Initialisation vectors 
% These vectors will be used to store the wave statistics at each position
Hrms_tot = zeros(Npos,1);  % root mean square height (m)
H13_tot  = zeros(Npos,1);  % significant wave height (m)
Hm_tot   = zeros(Npos,1);  % mean wave height (m)

% --------------------------------------
%     Computation of wave statistics
% --------------------------------------

for i=1:Npos  % loop on the positions
     %Frecuency is equal to 2 because the sampling ratio has aperiod of 0.5sec
    wave = zero_crossing(data(:,i),2)
    Hrms(i) = rms_height(data(:,i))
    H13(i) = significant_height(data(:,i))
    H(i) = mean(data(:,i)) %I DONT KNOW IF THIS IS CORRECT!!!!

end

% --------------------------------------
%                  Output
% --------------------------------------

% visualisation of outputs
% ?
% ?