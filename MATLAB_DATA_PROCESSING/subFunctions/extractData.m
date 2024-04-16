% This function extracts and processes data from the given directories for PDMS and PBS, 
% and then plots the information in a set of subplots.
%
% Inputs:
% fs - desired sampling frequency for resampling
% directoryPDMS - path to the directory containing PDMS CSV files
% directoryPBS - path to a single CSV file for PBS
%
% Outputs:
% Sample_data - data extracted and processed from test Sample CSV files
% PBS_data - data extracted and processed from the PBS CSV file
% structure will be a cell array where each row is a different test, the
% sub table 1: all data from CSV raw
% sub table 2: raw data from CSV only while injecting
% sub table 3: all data interpolated to specified frequency
% sub table 4: interpolated to specified frequency and only while injecting

% Revision log
% 10/18/2023 - HMV - update language from PDMS to "test sample", 

function [Sample_data, PBS_data] = extractData(fs, directorySample, directoryPBS)

    % Obtain a list of all the CSV files present in the directorySample
    files = dir(fullfile(directorySample, '*.csv'));

    % Initialize PDMS_data cell array to store raw and processed data
    Sample_data = cell(length(files), 2);

   %colormapping
    %color = ["r","g","b","r","g","b","r","g","b"]; 
    colors = customcolor(size(Sample_data,1)); 

    % Load data from the PBS CSV file 
    tablePBS = readtable(directoryPBS);
    
    % Filter rows of tablePBS where VOLUME_CONTROL is equal to 1
    % ==1 is when the system is actively injection
    filt_tablePBS = tablePBS(tablePBS.VOLUME_CONTROL == 1, :);
    
    % Identify and remove rows with duplicate 'TOTAL_TEST_TIME' values in tablePBS
    [~, idx_PBS] = unique(tablePBS.TOTAL_TEST_TIME);
    tablePBS = tablePBS(idx_PBS, :);

    % Compute the desired time step and create a new timeline based on the desired sampling frequency
    dt = 1/fs; 
    fs_timePBS = (min(tablePBS.TOTAL_TEST_TIME):dt:max(tablePBS.TOTAL_TEST_TIME))';

    % Initialize a new table to store resampled data for PBS
    timePBS_const_fs = table();
    timePBS_const_fs.TOTAL_TEST_TIME = fs_timePBS;

    % List of columns in the data table that need to be interpolated
    columns_to_interpolate = {'ADC_RAW_GAINED_VALUE', 'KPA', 'VOLUME_CONTROL', 'RATE', 'VOLUME', 'VOLUME_INSTANT_TRANSFER', 'VOLUME_TOTAL'}; 
    
    % Interpolate data for each of the desired columns
    for c = 1:numel(columns_to_interpolate)
        column = columns_to_interpolate{c};
        timePBS_const_fs.(column) = interp1(tablePBS.TOTAL_TEST_TIME, tablePBS.(column), fs_timePBS);
    end

    % Filter the interpolated data where VOLUME_CONTROL is equal to 1
    filt_timePBS_const_fs = timePBS_const_fs(timePBS_const_fs.VOLUME_CONTROL == 1, :);

    % Group all PBS-related data tables into a single cell array for output
    PBS_data = {tablePBS, filt_tablePBS, timePBS_const_fs, filt_timePBS_const_fs};

% %     % Set up a 2x2 subplot framework for plotting
% %     figure(1);
    titles = {'Pressure vs Time', 'Rate vs Time', 'Volume vs Time', 'Pressure vs Radius'};
    axs = [];
    for i=1:4
        axs = [axs subplot(2, 2, i)];
    end

    % Loop through each PDMS CSV file, extract, process, and plot the data
    for f = 1:length(files)

        % Read the current PDMS CSV file
        filename = fullfile(directorySample, files(f).name);
        tableSample = readtable(filename);

        % Conversion constant from nL to mm3
        nL_to_mm3 = 0.001; 

% %         % Plot 1: Pressure vs Time
% %         plot(axs(1), tableSample.TOTAL_TEST_TIME, tableSample.KPA,'Color', colors(f,:),'LineWidth',1.5);
% %         hold(axs(1), 'on');
% %         xlabel(axs(1), 'Time (s)');
% %         ylabel(axs(1), 'Pressure (kPa)');

        % Filter rows of tablePDMS where VOLUME_CONTROL is equal to 1
        filt_tableSample = tableSample(tableSample.VOLUME_CONTROL == 1, :);

% %         % Plot 4: Pressure vs Volume
% %         plot(axs(4), (((((3/4) * (filt_tableSample.VOLUME_TOTAL * nL_to_mm3)) / pi).^(1/3))), filt_tableSample.KPA,'Color', colors(f,:),'LineWidth',1.5);
% % %         plot(axs(4), (filt_tableSample.VOLUME_TOTAL), filt_tableSample.KPA,'LineWidth',1.5);
% %         hold(axs(4), 'on');
% %         xlabel(axs(4), 'Radius (mm)');
% %         ylabel(axs(4), 'Pressure (kPa)');

% %         % Plot 3: Volume vs Time
% %         plot(axs(3), tableSample.TOTAL_TEST_TIME, tableSample.VOLUME_TOTAL,'Color', colors(f,:),'LineWidth',1.5);
% %         hold(axs(3), 'on');
% %         xlabel(axs(3), 'Time (s)');
% %         ylabel(axs(3), 'Volume');
% % 
% %         % Plot 2: Rate vs Time
% %         plot(axs(2), tableSample.TOTAL_TEST_TIME, tableSample.RATE,'Color', colors(f,:),'LineWidth',1.5);
% %         hold(axs(2), 'on');
% %         xlabel(axs(2), 'Time (s)');
% %         ylabel(axs(2), 'Rate');

        % Store the raw and filtered data from the current PDMS file
        Sample_data{f, 1} = tableSample; % raw data
        Sample_data{f, 2} = filt_tableSample; % raw but only injecting

        % Identify and remove rows with duplicate 'TOTAL_TEST_TIME' values in tablePDMS
        [~, idx] = unique(tableSample.TOTAL_TEST_TIME);
        tableSample = tableSample(idx, :);

        % Compute the desired time step and create a new timeline based on the desired sampling frequency for PDMS
        dt = 1/fs;
        fs_timePDMS = (min(tableSample.TOTAL_TEST_TIME):dt:max(tableSample.TOTAL_TEST_TIME))';

        % Initialize a new table to store resampled data for PDMS
        timeSample_const_fs = table();
        timeSample_const_fs.TOTAL_TEST_TIME = fs_timePDMS;

        % Interpolate data for each of the desired columns for PDMS
        columns_to_interpolate = {'ADC_RAW_GAINED_VALUE', 'KPA', 'VOLUME_CONTROL', 'RATE', 'VOLUME', 'VOLUME_INSTANT_TRANSFER', 'VOLUME_TOTAL'};
        for c = 1:numel(columns_to_interpolate)
            column = columns_to_interpolate{c};
            timeSample_const_fs.(column) = interp1(tableSample.TOTAL_TEST_TIME, tableSample.(column), fs_timePDMS);
        end
        
        % Filter the interpolated data where VOLUME_CONTROL is equal to 1 (expansion)
        filt_timeSample_const_fs = timeSample_const_fs(timeSample_const_fs.VOLUME_CONTROL == 1, :);

        % Store the interpolated and filtered interpolated data from the current PDMS file
        Sample_data{f, 3} = timeSample_const_fs; % all data interpolated to specified sample frequency
        Sample_data{f, 4} = filt_timeSample_const_fs; % interpolated data a specified frequency only injecting
    end
    
    % Assign titles to each subplot after processing all files
    for i = 1:length(axs)
        title(axs(i), titles{i});
        hold(axs(i), 'off'); % Release the hold to avoid adding more plots to the current axis
    end
end