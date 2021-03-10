% Analysis for Ctrax'd Fly Bowl Data
% Abridged Ver1.0 for Emily Chen 
% 02/05/21 SJL 
% Updated 03/10/21 SJL 

%  imports and formats Ctrax data (after calculating per_frame_stats)
%  uses x,y position to plot velocity, walking fraction, total distance 

%  base data structure:
%  data = 1 x n cell, n = # videos
%  data{1,1}, struct for each video, each row is a fly
%  data{1,1}(1).x_mm x position for first fly
% ---------------------------------------------------------------------
%% Import and Format Ctrax Data 

% after tracking and generating per_frame_stats_xx.mat file, this script
% will import and group files (i.e. all files for same genotype) 

% start in directory where per_frame_stats_xx.mat files are 

% initialize variables 
name = 'abc'; % group name 
trial = 1.2; % trial type, will load relevant experimental settings 
save_option = 1 % 1 = save, 0 = don't save 

uipickfiles % choose set of files to group together 
test = ans; 

if trial == 1.2 % no stim (current) 
load('/Volumes/GoogleDrive/My Drive/Tuthill Lab Shared/SuYee/fly bowl/Nov 12/noatr_ctrl_males_led1_ts.mat') 
% ^ change to location of ts file from experiments 
end 

for i = 1:length(test) % iterate through each file 
    file = test{i}; 
    load(file); 
    if exist('all_trx') == 0
        all_trx = {trx}; 
    elseif exist('all_trx') == 1
        all_trx = [all_trx {trx}]; % group all_trx data from each file together 
    end 

data = all_trx; % saves each file as a cell in variable 'data' 
end 

if save_option == 1 
    save(name) 
end 

%% 

% 1. Plot stimuli/camera signals for reference 

toplot = 1; % set to 0, to turn off plotting 

figure
plot(t,ch(1,:)-2,'k','LineWidth',1) % piezo (black)
hold on
plot(t,ch(2,:)-4,'r','LineWidth',1.5) % led (red)
hold on
plot(t,ch(3,:)-6,'g','LineWidth',1.5) % camera (green)
xlabel('Time (s)','FontSize',12)
ylabel('Voltage (V)')
set(gca,'FontSize',12)
legend('Piezo', 'LED', 'Camera')

% 2. Initialize variables, load and package all data 

% Initialize some trial type specific variables, package data into workable
% format using packagedata function 

% packagedata output: 
            % data_matrix(1) = x_mm 
            % data_matrix(2) = y_mm 
            % data_matrix(3) = timestamps 

switch(trial) 
    
    case 1.2 %no stim, but use same variables as led  
    start_time = 0;
    frame_length = 2334;
    fps = 33.3;
    [data_matrix] = packagedata(data, frame_length, fps, start_time);

    otherwise 
        fprintf('no trial type found') 
    
end


% 3. Clean up data due to tracking error and calculate speed 

speed = [];
speed_mat = [];
distance = []; 
distance_mat = []; 

    % Define file name 
    for i = 1:length(data)
        matname = data{1,i}.matname; 
        underlineLocation = strfind(matname, name);
        filename = matname(underlineLocation:end-4); %defines filename
        trx = data{1,i};

    % Clean up data - for flies w/ incomplete tracking, pad spaces with nans 
        for flies = 1:length(trx) 
            nan_filler = nan(1,frame_length); % create frame length line of nans 
            if length(trx(flies).timestamps) < frame_length && length(trx(flies).timestamps)>1 ...
                    && trx(flies).timestamps(1) == 0 % if fly was tracked at time 0 for more than 1 frame, but not til end 
                nan_filler(1:length(trx(flies).timestamps)) = trx(flies).x_mm; % create new filler with x/y data and nan padding at end 
                nan_filler(1:length(trx(flies).timestamps)) = trx(flies).y_mm;
                trx(flies).x_mm = nan_filler; %replace x/y data with nan padded data 
                trx(flies).y_mm = nan_filler;
            end
        end

    % Calculate speed using Euclidean distance * fps 
    % Calculate distance using Euclidean distance 
        for ii = 1:length(data{1,i})
            speedcurr = sqrt(diff(trx(ii).x_mm).^2 + diff(trx(ii).y_mm).^2)*trx(ii).fps; % calculate speed 
            distancecurr = sqrt(diff(trx(ii).x_mm).^2 + diff(trx(ii).y_mm).^2);
            
    % filter out big jumps due to tracking error 
            for iii = 1:length(speedcurr) 
                if speedcurr(iii) > 80 % cut off of 80 mm/s, as this is impossible 
                    speedcurr(iii) = NaN; 
                end 
                if distancecurr(iii) > 2 
                    distancecurr(iii) = NaN; 
                end 
            end 
            
     % start building speed variable with data from all flies 
            if length(speedcurr) >= frame_length-1 
                speed{1,i}{ii} = {speedcurr};
                
     % for flies with incomplete tracking that couldn't be cleaned up, throw out      
            elseif length(speedcurr) < frame_length - 1 
                speed{1,i}{ii} = NaN; %NaN(1,frame_length-1); 
            end 
     % repeat for distance variable 
            if length(distancecurr) >= frame_length-1 
                distance{1,i}{ii} = {distancecurr}; 
                
            elseif length(speedcurr)< frame_length - 1 
                distance{1,i}{ii} = NaN; 
            end 
        end
    end 

% Convert speed / distance variables (in cell format) to giant num flies x frames
% matrix (speed_mat, distance_mat) 

for i = 1:length(speed) 
    for ii = 1:length(speed{1,i}) 
        if iscell(speed{1,i}{ii}) == 1 
        start = find(data{1,i}(ii).timestamps <= start_time,1, 'last'); %find last index where timestamps less than stim start time 
        speed_mat = [speed_mat; speed{1,i}{ii}{1,1}(1:frame_length-1)]; %build matrix of all speed data (easier to access than cells) 
        elseif iscell(speed{1,i}{ii}) == 0 
            speed_mat = [speed_mat; NaN(1,frame_length-1)]; 
        end 
    end 
end 

for i = 1:length(distance) 
    for ii = 1:length(distance{1,i}) 
        if iscell(distance{1,i}{ii}) == 1 
        distance_mat = [distance_mat; distance{1,i}{ii}{1,1}(1:frame_length-1)]; %build matrix of all distance data (easier to access than cells) 
        elseif iscell(distance{1,i}{ii}) == 0 
            distance_mat = [distance_mat; NaN(1,frame_length-1)]; 
        end 
    end 
end 

%% Plot Speed as 1) Time-Series and 2) PDF / Box Plot / Scatter 

close all 

switch(trial) 
    
    case{1.2} 
    
    % plot speed as time-series 
    figure
    stdshade(nanmean(speed_mat),0.2, 'b', 1:length(speed_mat),5); 
    xlim([0 1800])
    xticks([0 166.5 333 500 666 832.5 999 1165.5 1332 1498.5 1665 1831.5 1998 2164.5 2331 ])
    xticklabels({'0','5','10','15', '20', '25', '30', '35', '40', '45', '50', '55', '60', '65', '70'})
    ylim([0 45])
    xlabel('time(s)', 'fontlsize', 14)
    ylabel('avg speed (mm/s)', 'fontsize', 14)
  
    for flies = 1:size(speed_mat,1) 
        individual_avg = nanmean(speed_mat,2); 
    end 
    
    % plot speed as PDF, Summary (Box Plot), and Raw Observations (Scatter)
 
    figure 
    raincloud_plot(individual_avg, 'box_on', 1)
    xlim([-10 45])
    xlabel('speed (mm/s)') 

end 

%% Plot Distance as 1) Time-Series and 2) PDF / Box Plot / Scatter 

close all 

switch(trial) 
    
    case{1.2} 
    
    % plot speed as time-series 
    figure(1)
    stdshade(nanmean(distance_mat),0.2, 'b', 1:length(distance_mat),5); 
    xlim([0 1800])
    xticks([0 166.5 333 500 666 832.5 999 1165.5 1332 1498.5 1665 1831.5 1998 2164.5 2331 ])
    xticklabels({'0','5','10','15', '20', '25', '30', '35', '40', '45', '50', '55', '60', '65', '70'})
    ylim([0 1])
    xlabel('time(s)', 'fontsize', 14)
    ylabel('avg distance (mm)', 'fontsize', 14)
    title('Avg Distance Travelled of Group') 
  
    for flies = 1:size(distance_mat,1) 
        total_dist = sum(distance_mat,2); 
    end 
   
    % plot speed as PDF, Summary (Box Plot), and Raw Observations (Scatter
    
    figure(2)
    raincloud_plot(total_dist, 'box_on', 1)
    %xlim([-10 45])
    xlabel('distance (mm)') 
    title('Total Distance Travelled') 

end 


%% Maximum Speed - finds maximum speed per fly 

for flies = 1:size(speed_mat,1)
    individual_max = max(speed_mat,[ ], 2);
end

% plot speed as PDF, Summary (Box Plot), and Raw Observations (Scatter)
    
figure
raincloud_plot(individual_max, 'box_on', 1)
xlim([-10 80])
xlabel('max speed (mm/s)')

%% Walking Fraction  

for i = 1:size(speed_mat,1)
    walk_frac(i) = sum(speed_mat(i,:)>3)/length(speed_mat(i,:))*100;
end

figure
raincloud_plot(walk_frac, 'box_on', 1)
xlim([-30 110])
xlabel('walking frac')
    
