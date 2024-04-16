% Brendan M Unikewicz, PhD Student
% Andre M Pincot, MSc Student
% Tal Cohen, Asc. Professor
% MIT, Dept. Mechanical Engineering
% MIT, Dept. Civil & Environmental Engineering
% Date of Creation: 03/22/2024
% Code Purpose: creating notable figures and outcomes from validation

%% Code Start
clear; clc; format longg; close all; addpath('subFunctions');

%% Define Colors
redColor = [0.98, 0.40, 0.35]; 
greenColor = [0.45, 0.80, 0.69]; 
blueColor = [0.55, 0.60, 0.79]; 
c = [redColor; greenColor; blueColor];

fontSize = 14; 

%% Loading from preProcessedData 
preProcessed_data = [pwd,'\data_preProcessing'];
matFiles = dir(fullfile(preProcessed_data, '*.mat'));
for dataIndex = 1:length(matFiles)
    filePath = fullfile(preProcessed_data, matFiles(dataIndex).name);
    load(filePath);
end

%% Showing all Results -- Pressure versus Stretch
fs = 500;  
stretch_frac = 4.5; % Evaluate at fracture stretch of value, #

% Evaluate 43-1 results
for i=1:length(p_43_1)
    t_43_1{i} = (1/fs):(1/fs):(length(p_43_1{i})/fs); % Time scale
    diff_stretch_frac_43{i} = abs(stretch_43_1{i}-stretch_frac); % Stretch change relative to stretch_frac
    [~, index_stretch_frac_43{i}] = min(diff_stretch_frac_43{i}); % Minimum change in stretch
    closest_stretch_frac_43{i} = stretch_43_1{i}(index_stretch_frac_43{i}); % Closest stretch to stretch_frac
    p_stretch_frac_43{i} = p_43_1{i}(index_stretch_frac_43{i}); % Pressure at minimum change in stretch
    [pks_43{i},locs_43{i}] = findpeaks(p_43_1{i}); % Beginning findpks analysis for p_c
    [max_peak_43{i}, idx_43{i}] = max(pks_43{i}); % Max peak value
    V_T_43_1{i} = (4/3)*pi*(stretch_43_1{i}*initialDefect_43_1(i)).^3; % volume in uL
    V_c_43_1{i} = (4/3)*pi*(stretch_43_1{i}(locs_43{i}(idx_43{i}))*initialDefect_43_1(i)).^3; % volume at p_c
    p_c_stretch_43{i} = stretch_43_1{i}(locs_43{i}(idx_43{i})); % stretch at p_c
    U_43_V_c{i} = trapz((V_T_43_1{i}(1:(locs_43{i}(idx_43{i})))*1e-9),1e3*p_43_1{i}(1:(locs_43{i}(idx_43{i})))); % energy at p_c
    U_43_V_T_cumul{i} = cumtrapz((V_T_43_1{i}(1:end)*1e-9),1e3*p_43_1{i}(1:end)); % cumulative energy at volume total
    U_43_V_T{i} = trapz((V_T_43_1{i}(1:end)*1e-9),1e3*p_43_1{i}(1:end)); % energy at volume total
end

% cell2mat for min/max/mean for findpks analysis (p_c / V_c)
p_c_stretch_43 = cell2mat(p_c_stretch_43); 
p_c_43 = cell2mat(max_peak_43); 

% cell2mat for Energy analysis (p_c / V_c, & V_T)
U_43_V_c = cell2mat(U_43_V_c); 
U_43_V_T = cell2mat(U_43_V_T);
V_c_43_1 = cell2mat(V_c_43_1); 

% Evaluate 45-1 results
for i=1:length(p_45_1)
    t_45_1{i} = (1/fs):(1/fs):(length(p_45_1{i})/fs); % Time scale
    diff_stretch_frac_45{i} = abs(stretch_45_1{i}-stretch_frac); % Stretch change relative to stretch_frac
    [~, index_stretch_frac_45{i}] = min(diff_stretch_frac_45{i}); % Minimum change in stretch
    closest_stretch_frac_45{i} = stretch_45_1{i}(index_stretch_frac_45{i}); % Closest stretch to stretch_frac
    p_stretch_frac_45{i} = p_45_1{i}(index_stretch_frac_45{i}); % Pressure at minimum change in stretch
    [pks_45{i},locs_45{i}] = findpeaks(p_45_1{i}); % Beginning findpks analysis for p_c
    [max_peak_45{i}, idx_45{i}] = max(pks_45{i}); % Max peak value
    V_T_45_1{i} = (4/3)*pi*(stretch_45_1{i}*initialDefect_45_1(i)).^3; % volume in uL
    V_c_45_1{i} = (4/3)*pi*(stretch_45_1{i}(locs_45{i}(idx_45{i}))*initialDefect_45_1(i)).^3; % volume at p_c
    p_c_stretch_45{i} = stretch_45_1{i}(locs_45{i}(idx_45{i})); % stretch at p_c
    U_45_V_c{i} = trapz((V_T_45_1{i}(1:(locs_45{i}(idx_45{i})))*1e-9),1e3*p_45_1{i}(1:(locs_45{i}(idx_45{i})))); % energy at p_c
    U_45_V_T_cumul{i} = cumtrapz((V_T_45_1{i}(1:end)*1e-9),1e3*p_45_1{i}(1:end)); % cumulative energy at volume total
    U_45_V_T{i} = trapz((V_T_45_1{i}(1:end)*1e-9),1e3*p_45_1{i}(1:end)); % energy at volume total
end

% cell2mat for min/max/mean for findpks analysis (p_c / V_c)
p_c_stretch_45 = cell2mat(p_c_stretch_45); 
p_c_45 = cell2mat(max_peak_45); 

% cell2mat for Energy analysis (p_c / V_c, & V_T)
U_45_V_c = cell2mat(U_45_V_c); 
U_45_V_T = cell2mat(U_45_V_T);
V_c_45_1 = cell2mat(V_c_45_1); 

for i=1:length(p_47_1)
    t_47_1{i} = (1/fs):(1/fs):(length(p_47_1{i})/fs); % Time scale
    diff_stretch_frac_47{i} = abs(stretch_47_1{i}-stretch_frac); % Stretch change relative to stretch_frac
    [~, index_stretch_frac_47{i}] = min(diff_stretch_frac_47{i}); % Minimum change in stretch
    closest_stretch_frac_47{i} = stretch_47_1{i}(index_stretch_frac_47{i}); % Closest stretch to stretch_frac
    p_stretch_frac_47{i} = p_47_1{i}(index_stretch_frac_47{i}); % Pressure at minimum change in stretch
    [pks_47{i},locs_47{i}] = findpeaks(p_47_1{i}); % Beginning findpks analysis for p_c
    [max_peak_47{i}, idx_47{i}] = max(pks_47{i}); % Max peak value
    V_T_47_1{i} = (4/3)*pi*(stretch_47_1{i}*initialDefect_47_1(i)).^3; % volume in uL
    V_c_47_1{i} = (4/3)*pi*(stretch_47_1{i}(locs_47{i}(idx_47{i}))*initialDefect_47_1(i)).^3; % volume at p_c
    p_c_stretch_47{i} = stretch_47_1{i}(locs_47{i}(idx_47{i})); % stretch at p_c
    U_47_V_c{i} = trapz((V_T_47_1{i}(1:(locs_47{i}(idx_47{i})))*1e-9),1e3*p_47_1{i}(1:(locs_47{i}(idx_47{i})))); % energy at p_c
    U_47_V_T_cumul{i} = cumtrapz((V_T_47_1{i}(1:end)*1e-9),1e3*p_47_1{i}(1:end)); % cumulative energy at volume total
    U_47_V_T{i} = trapz((V_T_47_1{i}(1:end)*1e-9),1e3*p_47_1{i}(1:end)); % energy at volume total
end

% cell2mat for min/max/mean for findpks analysis (p_c / V_c)
p_c_stretch_47 = cell2mat(p_c_stretch_47); 
p_c_47 = cell2mat(max_peak_47); 

% cell2mat for Energy analysis (p_c / V_c, & V_T)
U_47_V_c = cell2mat(U_47_V_c); 
U_47_V_T = cell2mat(U_47_V_T);
V_c_47_1 = cell2mat(V_c_47_1); 

%% Plotting all Data for 43-1, 45-1, 47-1
figure; 
subplot(2,1,2);
for i=1:length(p_43_1)
    plot(stretch_43_1{i},p_43_1{i},'Color',redColor,'LineWidth',1.5); hold on; 
    plot(stretch_45_1{i},p_45_1{i},'Color',greenColor,'LineWidth',1.5); hold on; 
    plot(stretch_47_1{i},p_47_1{i},'Color',blueColor,'LineWidth',1.5); 
end
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
xlabel('Stretch ($$\lambda$$)', 'Interpreter', 'latex','FontSize',fontSize); ylabel('Pressure (kPa)', 'Interpreter', 'latex','FontSize',fontSize);
xlim([0, 5]); ylim([0,50]); legend('43:1', '45:1', '47:1','Interpreter','latex','FontSize',fontSize);
hold off; 

subplot(2,1,1);
for i=1:length(p_43_1)
    plot(((4/3)*1e3*pi*(stretch_43_1{i}*initialDefect_43_1(i)).^3) ./ 1000,p_43_1{i},'Color',redColor,'LineWidth',1.5); hold on; 
    plot(((4/3)*1e3*pi*(stretch_45_1{i}*initialDefect_45_1(i)).^3) ./ 1000,p_45_1{i},'Color',greenColor,'LineWidth',1.5); hold on; 
    plot(((4/3)*1e3*pi*(stretch_47_1{i}*initialDefect_47_1(i)).^3) ./ 1000,p_47_1{i},'Color',blueColor,'LineWidth',1.5); hold on; 
end
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
xlabel('Volume ($\mu$L)', 'Interpreter', 'latex','FontSize',fontSize); ylabel('Pressure (kPa)','Interpreter','latex','FontSize',fontSize);
xlim([0, 9]); ylim([0,50]); legend('43:1', '45:1', '47:1','Interpreter','latex','FontSize',fontSize-1);
hold off; 

%% Signal Smoothing to determine peak-to-peak noise values
binStart_stretch = 0:0.25:6.25; % Starts at 0, ends at 6.25, increment by 0.25
binEnd_stretch = 0.25:0.25:6.50; % Starts at 0.25, ends at 6.50, increment by 0.25

binStart_time = 0:0.5:30.5; % Time-bins (s)
binEnd_time = 0.5:0.5:31.0;

% Calculating noise
for i=1:length(p_43_1)
    smooth_43_1{i} = movmean(p_43_1{i},99); % moving mean smoothed signal
    noise_43_1{i} = p_43_1{i} - smooth_43_1{i}; % subtract out smoothed signal
    snr_43_1{i} = p_43_1{i} / abs(noise_43_1{i}); % signal-to-noise ratio
    meanNoise_43_1{i} = mean(noise_43_1{i}); % average noise
    meanSNR_43_1{i} = mean(snr_43_1{i}); % average signal-to-noise ratio
end

% Binning procedure for plotting (stretch)
for i = 1:length(stretch_43_1)
    for j = 1:length(binStart_stretch)
        binIndices = find(stretch_43_1{i} >= binStart_stretch(j) & stretch_43_1{i} < binEnd_stretch(j));
        stretch_43_bins{i, j} = stretch_43_1{i}(binIndices);
        noise_43_bins{i, j} = noise_43_1{i}(binIndices);
    end
end

% Binning procedure for plotting (stretch)
for i = 1:length(t_43_1)
    for j = 1:length(binStart_time)
        binIndices = find(t_43_1{i} >= binStart_time(j) & t_43_1{i} < binEnd_time(j));
        t_43_bins{i, j} = t_43_1{i}(binIndices);
        noiseTime_43_bins{i, j} = noise_43_1{i}(binIndices);
    end
end

meanNoise_43_bins = NaN(length(stretch_43_1), length(binStart_stretch)-1);
meanNoise_time_43_bins = NaN(length(t_43_1), length(binStart_time)-1);

for i = 1:length(stretch_43_1)
    for j = 1:length(binStart_stretch)-1  % Adjust to exclude the last binStart, since it has no end
        if ~isempty(noise_43_bins{i,j})
            meanNoise_43_bins(i,j) = mean(abs(noise_43_bins{i,j}) * 2);
        end
    end
end

for i = 1:length(t_43_1)
    for j = 1:length(binStart_time)-1  % Adjust to exclude the last binStart, since it has no end
        if ~isempty(noiseTime_43_bins{i,j})
            meanNoise_time_43_bins(i,j) = mean(abs(noiseTime_43_bins{i,j}) * 2);
        end
    end
end

avg_meanNoise_43_bins = mean(meanNoise_43_bins, 1, 'omitnan');  
avg_meanNoise_time_43_bins = mean(meanNoise_time_43_bins, 1, 'omitnan');  

barCenters = (binStart_stretch(2:end) + binEnd_stretch(1:end-1)) / 2; 
barCenters_time = (binStart_time(2:end) + binEnd_time(1:end-1)) / 2; 

% 45:1
for i=1:length(p_45_1)
    smooth_45_1{i} = movmean(p_45_1{i},99); 
    noise_45_1{i} = p_45_1{i} - smooth_45_1{i}; % A vector
    snr_45_1{i} = p_45_1{i} / abs(noise_45_1{i});
    meanNoise_45_1{i} = mean(noise_45_1{i}); % Or use rms function.
    meanSNR_45_1{i} = mean(snr_45_1{i});
end

for i = 1:length(stretch_45_1)
    for j = 1:length(binStart_stretch)
        binIndices = find(stretch_45_1{i} >= binStart_stretch(j) & stretch_45_1{i} < binEnd_stretch(j));
        stretch_45_bins{i, j} = stretch_45_1{i}(binIndices);
        noise_45_bins{i, j} = noise_45_1{i}(binIndices);
    end
end

for i = 1:length(t_45_1)
    for j = 1:length(binStart_time)
        binIndices = find(t_45_1{i} >= binStart_time(j) & t_45_1{i} < binEnd_time(j));
        t_45_bins{i, j} = t_45_1{i}(binIndices);
        noiseTime_45_bins{i, j} = noise_45_1{i}(binIndices);
    end
end

meanNoise_45_bins = NaN(length(stretch_45_1), length(binStart_stretch)-1);
meanNoise_time_45_bins = NaN(length(t_45_1), length(binStart_time)-1);

for i = 1:length(stretch_45_1)
    for j = 1:length(binStart_stretch)-1  % Adjust to exclude the last binStart, since it has no end
        if ~isempty(noise_45_bins{i,j})
            meanNoise_45_bins(i,j) = mean(abs(noise_45_bins{i,j}) * 2);
        end
    end
end

for i = 1:length(t_45_1)
    for j = 1:length(binStart_time)-1  % Adjust to exclude the last binStart, since it has no end
        if ~isempty(noiseTime_45_bins{i,j})
            meanNoise_time_45_bins(i,j) = mean(abs(noiseTime_45_bins{i,j}) * 2);
        end
    end
end

avg_meanNoise_45_bins = mean(meanNoise_45_bins, 1, 'omitnan');  
avg_meanNoise_time_45_bins = mean(meanNoise_time_45_bins, 1, 'omitnan');  

barCenters = (binStart_stretch(2:end) + binEnd_stretch(1:end-1)) / 2; 
barCenters_time = (binStart_time(2:end) + binEnd_time(1:end-1)) / 2; 

for i=1:length(p_47_1)
    smooth_47_1{i} = movmean(p_47_1{i},99); 
    noise_47_1{i} = p_47_1{i} - smooth_47_1{i}; % A vector
    snr_47_1{i} = p_47_1{i} / abs(noise_47_1{i});
    meanNoise_47_1{i} = mean(noise_47_1{i}); % Or use rms function.
    meanSNR_47_1{i} = mean(snr_47_1{i});
end

for i = 1:length(stretch_47_1)
    for j = 1:length(binStart_stretch)
        binIndices = find(stretch_47_1{i} >= binStart_stretch(j) & stretch_47_1{i} < binEnd_stretch(j));
        stretch_47_bins{i, j} = stretch_47_1{i}(binIndices);
        noise_47_bins{i, j} = noise_47_1{i}(binIndices);
    end
end

for i = 1:length(t_47_1)
    for j = 1:length(binStart_time)
        binIndices = find(t_47_1{i} >= binStart_time(j) & t_47_1{i} < binEnd_time(j));
        t_47_bins{i, j} = t_47_1{i}(binIndices);
        noiseTime_47_bins{i, j} = noise_47_1{i}(binIndices);
    end
end

meanNoise_47_bins = NaN(length(stretch_47_1), length(binStart_stretch)-1);
meanNoise_time_47_bins = NaN(length(t_47_1), length(binStart_time)-1);

for i = 1:length(stretch_47_1)
    for j = 1:length(binStart_stretch)-1  % Adjust to exclude the last binStart, since it has no end
        if ~isempty(noise_47_bins{i,j})
            meanNoise_47_bins(i,j) = mean(abs(noise_47_bins{i,j}) * 2);
        end
    end
end

for i = 1:length(t_47_1)
    for j = 1:length(binStart_time)-1  % Adjust to exclude the last binStart, since it has no end
        if ~isempty(noiseTime_47_bins{i,j})
            meanNoise_time_47_bins(i,j) = mean(abs(noiseTime_47_bins{i,j}) * 2);
        end
    end
end

avg_meanNoise_47_bins = mean(meanNoise_47_bins, 1, 'omitnan');  
avg_meanNoise_time_47_bins = mean(meanNoise_time_47_bins, 1, 'omitnan');  

barCenters = (binStart_stretch(2:end) + binEnd_stretch(1:end-1)) / 2; 
barCenters_time = (binStart_time(2:end) + binEnd_time(1:end-1)) / 2; 

%% Mean Noise Distribution for 43, 45, 47:1 Samples -- Stretch
figure;
noise_hGroup = [avg_meanNoise_43_bins; avg_meanNoise_45_bins; avg_meanNoise_47_bins]; 
noise_hGroup_time = [avg_meanNoise_time_43_bins; avg_meanNoise_time_45_bins; avg_meanNoise_time_47_bins];
b = bar(barCenters, noise_hGroup','grouped'); hold on; 
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
legend('43:1', '45:1', '47:1', 'Interpreter', 'latex','FontSize',fontSize-1);
set(b(1),'FaceColor',redColor,'EdgeColor',redColor); set(b(2),'FaceColor',greenColor,'EdgeColor',greenColor); set(b(3),'FaceColor',blueColor,'EdgeColor',blueColor);
xlabel('Stretch ($\lambda$)', 'Interpreter', 'latex','FontSize',fontSize);
ylabel('Mean Noise (kPa)', 'Interpreter', 'latex','FontSize',fontSize);
hold off; 

%% Mean Noise Distribution for 43, 45, 47:1 Samples -- Time
figure;
noise_hGroup_time = [avg_meanNoise_time_43_bins; avg_meanNoise_time_45_bins; avg_meanNoise_time_47_bins];
b_t = bar(barCenters_time, noise_hGroup_time', 'grouped'); hold on; 
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
legend('43:1', '45:1', '47:1','Interpreter','latex','FontSize',fontSize-1);
set(b_t(1),'FaceColor',redColor); set(b_t(2),'FaceColor',greenColor); set(b_t(3),'FaceColor',blueColor)
xlim([0 15]);
xlabel('Time (s)', 'Interpreter', 'latex','FontSize',fontSize);
ylabel('Mean Noise (kPa)', 'Interpreter', 'latex','FontSize',fontSize);
hold off; 

%% Elastic Moduli
data_iE = [E_i_43_1(:), E_i_45_1(:), E_i_47_1(:)]; 

%% Initial Defect Volumes
V_i_43_1 = (4./3).*pi.*(initialDefect_43_1.^3);
V_i_45_1 = (4./3).*pi.*(initialDefect_45_1.^3);
V_i_47_1 = (4./3).*pi.*(initialDefect_47_1.^3);

%% Defining Energies for Box-plotting
data_U_V_c = [U_43_V_c(:), U_45_V_c(:), U_47_V_c(:)]; 
data_U_V_T = [U_43_V_T(:), U_45_V_T(:), U_47_V_T(:)]; 

%% Showing neo-Hookean Fitting Process -- 43:1
L=length(stretch_43_1); 

figure;
for i = 1:L
    plot(stretch_43_1{i}, p_43_1{i}, 'Color',redColor, 'LineWidth', 1.5); 
    hold on; 
    plot(fittingStretch_43_1{i}, pressure_regionFitting_43_1{i}, 'Color',blueColor, 'LineWidth', 1.5);
    hold on;
    plot(stretch_43_1{i}, fullFit_dataY_43_1{i}, 'k--','LineWidth', 1.5);
    hold on;
end
xlim([1 3]); ylim([0 60]);
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
xlabel('Stretch ($\lambda$)','Interpreter','latex','FontSize',fontSize); ylabel('Pressure (kPa)','Interpreter','latex','FontSize',fontSize);
legend('43:1','Fitting Region','neo-Hookean Fit','Location','northwest','Interpreter','latex','FontSize',fontSize-1);
hold off; 

%% Plotting FFT of 250uL vs 010uL syringe for 43-1 
O = length(stretch_43_Wo); 
fs = 500; 
T = 1/fs; 

figure; 
subplot(2,1,1); 
for i=O
    plot(stretch_43_Wo{i},p_43_Wo{i},'k','LineWidth',1.5); hold on; 
end

for i=1:L
    plot(stretch_43_1{i},p_43_1{i},'Color',redColor,'LineWidth',1.5); hold on; 
end
% Duplicate Wo to overlay first
for i=O
    plot(stretch_43_Wo{i},p_43_Wo{i},'k','LineWidth',1.5); hold on; 
end
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
xlabel('Stretch ($\lambda$)', 'Interpreter', 'latex','FontSize',fontSize); ylabel('Pressure (kPa)', 'Interpreter', 'latex','FontSize',fontSize);
legend('43:1 (250$\mu$L)','43:1 (10$\mu$L)','Interpreter','latex','FontSize',fontSize-1);
xlim([0.5, 4]); ylim([0, 50]);
hold off; 

subplot(2,1,2); 
for i = 1:L
    % Compute FFT
    L1{i} = length(p_43_1{i});
    t1{i} = (0:L1{i}-1)*T;        % Time vector
    y1{i} = fft(p_43_1{i});
    % Compute two-sided spectrum and single-sided spectrum
    P21{i} = abs(y1{i}/L1{i});
    P11{i} = P21{i}(1:floor(L1{i}/2)+1);
    P11{i}(2:end-1) = 2*P11{i}(2:end-1);
    % Define frequency domain
    f1{i} = fs*(0:(L1{i}/2))/L1{i};
    % Plot single-sided amplitude spectrum.
    plot(f1{i}, P11{i}, 'Color', redColor, 'LineWidth',1.5); hold on; 
end

for i = 1:O
    % Compute FFT
    L2{i} = length(p_43_Wo{i});
    t2{i} = (0:L2{i}-1)*T;        % Time vector
    y2{i} = fft(p_43_Wo{i});
    % Compute two-sided spectrum and single-sided spectrum
    P22{i} = abs(y2{i}/L2{i});
    P12{i} = P22{i}(1:L2{i}/2+1);
    P12{i}(2:end-1) = 2*P12{i}(2:end-1);
    % Define frequency domain
    f2{i} = fs*(0:(L2{i}/2))/L2{i};
    % Plot single-sided amplitude spectrum.
    plot(f2{i}, P12{i}, 'k', 'LineWidth',1.5); hold on; 
end
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
xlabel('Frequency (Hz)', 'Interpreter', 'latex','FontSize',fontSize);
ylabel('Pressure (kPa)', 'Interpreter', 'latex','FontSize',fontSize);
xlim([0 15]); ylim([0 2])
hold on; 

% Define the positions for the insets in normalized units
insetPosition1 = [0.25 0.2 0.2 0.2]; 
insetPosition2 = [0.65 0.2 0.2 0.2]; 

% Create the first inset plot for 0.225Hz
axes('Position', insetPosition1);
box on;
for i = 1:L
    plot(f1{i}, P11{i}, 'Color',redColor, 'LineWidth',1.5);
    hold on; ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
end
for i = 1:O
    plot(f2{i}, P12{i}, 'k', 'LineWidth', 1);
    ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
end
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
xlim([0.0 0.46]); ylim([0.5 2.25]);

% Create the second inset plot for 11.35Hz
axes('Position', insetPosition2);
box on; 
for i = 1:L
    plot(f1{i}, P11{i}, 'Color',redColor, 'LineWidth',1.5);
    hold on; ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
end
for i = 1:O
    plot(f2{i}, P12{i}, 'k', 'LineWidth',1.5);
    ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
end
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
xlim([10.6 12.2]); ylim([0 0.2]); 

%% Box-Plot for Instantaneous Elastic Modulus, E_i 
numRows_E_i = size(data_iE, 1); % Number of rows
numGroups_E_i = size(data_iE, 2); % Number of groups, which is 3 in your case

% Create group_inx for a long format
group_inx_E_i = repmat(1:numGroups_E_i, numRows_E_i, 1);
group_inx_E_i = group_inx_E_i(:); % Make it a single column

% Reshape your data to long format
data_long_E_i = reshape(data_iE, [], 1); % Transpose data_iE to stack columns vertically

condition_names = {'43:1','45:1','47:1'}; 
group_names = {'43:1','45:1','47:1'}; 

figure; 
h = daboxplot(data_long_E_i,'groups',group_inx_E_i,'colors', c,'xtlabels', condition_names,...
    'colors',c,'fill',0,'whiskers',0,'scatter',2,'jitter',1,'outsymbol','k*',...
    'outliers',1,'scattersize',48,'flipcolors',1,'boxspacing',1.2,...
    'legend',group_names); 
ylabel('Instantaneous Elastic Modulus (kPa)','Interpreter','latex','FontSize',fontSize);
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
%xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
set(gca, 'FontSize', fontSize); box on; 
xlabel('PDMS Ratio','Interpreter','latex','FontSize',fontSize);
legend(fliplr(group_names),'Interpreter','latex','FontSize',fontSize-1);

%% Energy vs. Sample Distribution (V_c)
numRows_U_V_c = size(data_U_V_c, 1); % Number of rows
numGroups_U_V_c = size(data_U_V_c, 2); % Number of groups, which is 3 in your case

% Create group_inx for a long format
group_inx_U_V_c = repmat(1:numGroups_U_V_c, numRows_U_V_c, 1);
group_inx_U_V_c = group_inx_U_V_c(:); % Make it a single column

% Reshape your data to long format
data_long_U_V_c = reshape(data_U_V_c.*1e6, [], 1); % Transpose data_iE to stack columns vertically

condition_names = {'43:1','45:1','47:1'}; 
group_names = {'43:1','45:1','47:1'}; 

figure; 
subplot(2,1,2)
h = daboxplot(data_long_U_V_c,'groups',group_inx_U_V_c,'colors', c,'xtlabels', condition_names,...
    'colors',c,'fill',0,'whiskers',0,'scatter',2,'jitter',1,'outsymbol','*',...
    'outliers',1,'scattersize',48,'flipcolors',1,'boxspacing',1.2,...
    'legend',group_names); 
ylabel('Energy ($\mu$J)','Interpreter','latex','FontSize',fontSize);
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
xl_U_V_c = xlim; xlim([xl_U_V_c(1)-0.1, xl_U_V_c(2)+0.2]); % make more space for the legend
set(gca, 'FontSize', fontSize); box on; xlabel('PDMS Ratio','Interpreter','latex','FontSize',fontSize);
legend(fliplr(group_names),'Interpreter','latex','FontSize',fontSize-1);

%% Energy vs. Sample Distribution (V_T)
numRows_U_V_T = size(data_U_V_T, 1); % Number of rows
numGroups_U_V_T = size(data_U_V_T, 2); % Number of groups, which is 3 in your case

% Create group_inx for a long format
group_inx_U_V_T = repmat(1:numGroups_U_V_c, numRows_U_V_T, 1);
group_inx_U_V_T = group_inx_U_V_T(:); % Make it a single column

% Reshape your data to long format
data_long_Q_V_T = reshape(data_U_V_T.*1000, [], 1); % Transpose data_iE to stack columns vertically

condition_names = {'43:1','45:1','47:1'}; 
group_names = {'43:1','45:1','47:1'}; 

subplot(2,1,1); 
h = daboxplot(data_long_Q_V_T,'groups',group_inx_U_V_T,'colors', c,'xtlabels', condition_names,...
    'colors',c,'fill',0,'whiskers',0,'scatter',2,'jitter',1,'outsymbol','k*',...
    'outliers',1,'scattersize',48,'flipcolors',1,'boxspacing',1.2,...
    'legend',group_names); 
ylabel('Energy (mJ)','Interpreter','latex','FontSize',fontSize);
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
xl_U_V_T = xlim; xlim([xl_U_V_T(1)-0.1, xl_U_V_T(2)+0.2]); % make more space for the legend
set(gca, 'FontSize', fontSize); box on; xlabel('PDMS Ratio','Interpreter','latex','FontSize',fontSize);
legend(fliplr(group_names),'Interpreter','latex','FontSize',fontSize-1);

%% Energy vs. Volume (V_c & V_T)
meanMean_volFracture = mean([mean(V_c_43_1),mean(V_c_45_1),mean(V_c_47_1)]); 
mean_V_c_43 = mean(V_c_43_1); 
mean_V_c_45 = mean(V_c_45_1); 
mean_V_c_47 = mean(V_c_47_1);

figure; 
subplot(2,1,1)
for i=1:length(p_43_1)
    plot(V_T_43_1{i},U_43_V_T_cumul{i}.*1000,'Color',redColor,'LineWidth',1.5); hold on; 
    plot(V_T_45_1{i},U_45_V_T_cumul{i}.*1000,'Color',greenColor,'LineWidth',1.5); hold on; 
    plot(V_T_47_1{i},U_47_V_T_cumul{i}.*1000,'Color',blueColor,'LineWidth',1.5); hold on; 
    xline(meanMean_volFracture,'--k','LineWidth',1.5); hold on; 
end 
legend('43:1','45:1','47:1','$<V_c>$','Interpreter','latex','FontSize',fontSize-1,'Location','northwest');
xlim([0 9]); 
xlabel('Volume ($\mu$L)','Interpreter','latex','FontSize',fontSize); 
ylabel('Energy (mJ)','Interpreter','latex','FontSize',fontSize); 
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 

subplot(2,1,2)
for i=1:length(p_43_1)
    plot(V_T_43_1{i}*1000,U_43_V_T_cumul{i}.*1e6,'Color',redColor,'LineWidth',1.5); hold on; 
    plot(V_T_45_1{i}*1000,U_45_V_T_cumul{i}.*1e6,'Color',greenColor,'LineWidth',1.5); hold on; 
    plot(V_T_47_1{i}*1000,U_47_V_T_cumul{i}.*1e6,'Color',blueColor,'LineWidth',1.5); hold on; 
    xline(mean_V_c_43*1000,'--','Color',redColor,'LineWidth',1.5); hold on; 
    xline(mean_V_c_45*1000,'--','Color',greenColor,'LineWidth',1.5); hold on; 
    xline(mean_V_c_47*1000,'--','Color',blueColor,'LineWidth',1.5); hold on; 
end
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
xlim([0 0.4*1000]); xlabel('Volume (nL)','Interpreter','latex','FontSize',fontSize);
ylabel('Energy ($\mu$J)','Interpreter','latex','FontSize',fontSize);
legend('43:1','45:1','47:1','43:1 $V_c$','45:1 $V_c$','47:1 $V_c$','Interpreter','latex','FontSize',fontSize-1,'Location','northwest');

%% Defining Relevant Data Variables
varList = {'p_c_43','p_c_45','p_c_47','p_c_stretch_43','p_c_stretch_45','p_c_stretch_47', ...
    'V_c_43_1','V_c_45_1','V_c_47_1','V_i_43_1','V_i_45_1','V_i_47_1', ...
    'V_T_43_1','V_T_45_1','V_T_47_1','t_43_1','t_45_1','t_47_1', ...
    'U_43_V_c','U_45_V_c','U_47_V_c','U_43_V_T','U_45_V_T','U_47_V_T', ...
    'U_43_V_T_cumul', 'U_45_V_T_cumul', 'U_47_V_T_cumul', ...
    'p_stretch_frac_43','p_stretch_frac_45','p_stretch_frac_47', ...
    'p_c_stretch_43','p_c_stretch_45','p_c_stretch_47'};

folderName = 'data_postProcessing';

% Check if the folder exists, if not, create it
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

%% Saving Relevant Data as .mat files for further plotting / post-processing
for i = 1:length(varList)
    varName = varList{i}; 
    fileName = fullfile(folderName, [varName '.mat']);
    eval(['save(fileName, ''' varName ''');']);
end

