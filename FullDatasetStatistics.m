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
Hrms_tot = zeros(Npos,1);  % root mean square height (m)
H13_tot  = zeros(Npos,1);  % significant wave height (m)
Hm_tot   = zeros(Npos,1);  % mean wave height (m)

% --------------------------------------
%     Computation of wave statistics
% --------------------------------------

for j=1:3 %loop for low, mid and high tides:
    for i=1:Npos  % loop on the positions
         %Frecuency is equal to 2 because the sampling ratio has aperiod of 0.5sec
        wave = zero_crossing(data(:,i+(j-1)*5),2)
        Hrms(i,j) = rms_height(data(:,i+(j-1)*5))
        H13(i,j) = significant_height(data(:,i+(j-1)*5))
        Hm(i,j) = mean(data(:,i+(j-1)*5)) %I DONT KNOW IF THIS IS CORRECT (Idk if by saying the water height they mean the values of the lowtide.txt)!!!!
    
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
    plot(Hrms(:,j),H13(:,j),color(j))
    hold on

    xlabel("Significant wave height")
    ylabel("Root mean square wave height")
    legend("Low tide","Mid tide","High tide")
end
for j = 1:3
    p1(j,:) = polyfit(Hrms(:,j),H13(:,j),1); %in p we store the slope p(1) and the intercept p(2)
    x1 = linspace(0.1,0.55);
    y1 = polyval(p1(j,:),x1);
    plot(x1,y1,colorp(j))


    %Now:  Verify that we have Hm ∼ 0.89Hrms.
    %There are some problems with the Hm
    p2(j,:) = polyfit(Hm(:,j),Hrms(:,j),1)
    x1 = linspace(0,0.35); % I am not sure if Hm is what I said
    y1 = polyval(p2(j,:),x1);
    plot(x1,y1,colorp(j))

end
disp(p1)%This is the fit of the significant wave height as a function of the rms wave height %Not that close to sqrt(2) 
disp(p2)





%This things below i think we dont have to do them for the this script (mid and high)
% positions = [4478, 4765, 4790, 4814, 4835] %Positions where the sensors are located
% subplot(2,1,1)
% plot(positions, Hrms,"or")
% title("Mean of RMS, wave height and H 1/3 ")
% 
% hold on
% plot(positions, H13,"*b")
% plot(positions, Hm,"*g") 
% legend(["Hrms","H1/3","Height"])
% xlabel("Position (m)")
% 
% subplot(2,1,2)
% plot(prof(:,1),prof(:,2),"black")
% title("Depth profile")
% 
% xlabel("Position (m)")
% ylabel("Elevation (m)")