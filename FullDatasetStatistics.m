%FULL DATASET STATISTICS
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
datalow = load('lowTide.txt');
datamid = load('midTide.txt');
datahigh = load('highTide.txt');
data = [datalow, datamid, datahigh]

prof = load("prof1018.txt");

% constants
g = 9.81;                    % acceleration of gravity (m/s^2)
rho = 1025;                  % water density (kg/m^3)
fs = 2;                      % sampling frequency (Hz)
Npos = 5;                    % number of cross-shore positions considered
% More? 

% Initialisation vectors 
% These vectors will be used to store the wave statistics at each position
Hrms_tot = zeros(5,3);  % root mean square height (m)
H13_tot  = zeros(5,3);  % significant wave height (m)
Hm_tot   = zeros(5,3);  % mean wave height (m)

% --------------------------------------
%     Computation of wave statistics
% --------------------------------------

for j=1:3 %loop for low, mid and high tides:
    for i=1:Npos  % loop on the positions
         %Frecuency is equal to 2 because the sampling ratio has a period of 0.5sec
        wave = zero_crossing(data(:,i+(j-1)*5),2);
        Hrms_tot(i,j) = rms_height(wave(:,1));
        H13_tot(i,j) = significant_height(wave(:,1));
        Hm_tot(i,j) = mean(wave(:,1)); %I DONT KNOW IF THIS IS CORRECT (Idk if by saying the water height they mean the values of the lowtide.txt)!!!!
    
    end


% --------------------------------------
%                  Output
% --------------------------------------

% visualisation of outputs
% ?
% ?
color = ["*r","*b","*g"];
colorp = ["r","b","g"];   


    figure(1)
    plot(Hrms_tot(:,j),H13_tot(:,j),color(j))
    hold on

    xlabel("Root mean square wave height (m)")
    ylabel("Significant wave height (m)")
    legend("Low tide","Mid tide","High tide")
end
for j = 1:3
    p1(j,:) = polyfit(Hrms_tot(:,j),H13_tot(:,j),1); %in p we store the slope p(1) and the intercept p(2)
    x1 = linspace(0.35,1.55);
    y1 = polyval(p1(j,:),x1);
    plot(x1,y1,colorp(j))


    %Now:  Verify that we have Hm ∼ 0.89Hrms.
    %There are some problems with the Hm
    p2(j,:) = polyfit(Hrms_tot(:,j),Hm_tot(:,j),1);
    %x1 = linspace(0,1.55); % I am not sure if Hm is what I said
    %y1 = polyval(p2(j,:),x1);
    %plot(x1,y1,colorp(j))

end
fprintf("The values of the slope H13/Hrms for low, mid and high tides are:\n")%This is the fit of the significant wave height as a function of the rms wave height %Not that close to sqrt(2) 
disp(p1(:,1))

fprintf("The values of the slope Hm/Hrms for low, mid and high tides are:\n")
disp(p2(:,1))

save("StatisticsEgmond","Hm_tot","Hrms_tot","H13_tot");

