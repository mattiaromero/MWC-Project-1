function [wave]=zero_crossing(data,frequency)
%function [wave]=zero_crossing(data,frequency)
% Zero crossing analysis of wave data
% input   data: input array of water elevation in m. Any linear trend or
%         mean will be removed.
%         frequency: sampling frequency of data in Hz
% output  wave: array containing the individual wave characteristics.
%         It has the following columns 1: wave height (m), 2:
%         crest elevation (m), 3: trough elevation (m), 4: period (s).
%
% A threshold is defined to filter out too small waves. It is equal to 1% of the maximal wave height Hmax of the signal.
% Waves with a crest or through elavation less than this threshold are not considered.
%
% adapted from Urs Neumeier, version 1.06

if frequency <= 0
    error('Frequency must be greather than zero')
end

% The function was written for zero upward-crossing.
% To have zero downward-crossing (recommended) leave the next line uncommented
data=-data;

% ---- Calculation of the zero crossing
data=detrend(data);                 % find zero crossing avoiding zero values
d0=data(data~=0);
back0=1:length(data);
back0=back0(data~=0);
f=find(d0(1:end-1).*d0(2:end)<0);
crossing=back0(f);
if data(1)>0                        % reject first crossing if it is downward
    crossing(1)=[];
end
crossing=crossing(1:2:end);         % this are the zero up-ward crossing


% ---- Calculate crest, trough and period of each wave
wave=zeros(length(crossing)-1,4);   % initialisation of the array wave, 4 columns with 
                                    % wave height, wave crest, wave trough and wave period
for i=1:length(crossing)-1         
    wave(i,2)= max(data(crossing(i):crossing(i+1)));    % crest height
    wave(i,3)= -min(data(crossing(i):crossing(i+1)));   % trough height
end

if size(wave,1) >= 1   % if no wave was found, do nothing
    wave(:,4)=diff(crossing')/frequency;                % wave period
    
    % Definition of a threshold to filter out the waves which are too small 
    threshold=0.01*max(wave(:,2)+wave(:,3));            % threshold = 1% Hmax
    
    i=0;                                % remove waves that are too small
    while i < size(wave,1)              % by joining them to adjacent wave
        i=i+1;
        if wave(i,2)<threshold
            if i~=1
                wave(i-1,2:4)=[max(wave(i-1:i,2:3)) sum(wave(i-1:i,4))];
            end
            wave(i,:)=[];
        elseif wave(i,3)<threshold
            if i~=size(wave,1)
                wave(i,2:4)=[max(wave(i:i+1,2:3)) sum(wave(i:i+1,4))];
                wave(i+1,:)=[];
            else
                wave(i,:)=[];
            end
        end
    end
    
    % wave has 1: wave height, 2: wave crest (Hcm), 3: wave trough (Htm), 4: period.
    wave(:,1)= sum(wave(:,2:3),2);
    
end


