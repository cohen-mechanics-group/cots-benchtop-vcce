% This function processes and truncates & calibrates data for a test sample, by subtracting out PBS data.
% It returns calibrated/truncated to only injection (TC) test sample pressure data (kPa),
% truncated Sample volume data (nL),
% and truncated (T) (not calibrated) test sample pressure data (kPa).

% Revision history
%10/16/2023-HMV- renamed references from PDMS to "test sample" or xxxSample
%instead of "PDMS" to make code more versatile. Updated pressure balance
%step of calibration to accomidate non zero initial pressure and added a
%data validation prompt if the window for data validation has a large range
%of pressures.

function [pressureSample_TC, nLSample_T, pressureSample_T, timeSample_T, pressureSample_FC, nLSample_F, timeSample_F, uniqueVol_PBS, uniquePressure_PBS] = truncate_calibrateData(dataSample, dataPBS)

% Define the volume increment for calculations
dV = 0.01;

% Extract unique VOLUME_TOTAL values from dataPBS and return the index mapping
[unique_dataPBS_V, ~, index_mapping] = unique(dataPBS{1,4}.VOLUME_TOTAL);
% Using the index mapping, create a new array of pressures from dataPBS
unique_dataPBS_P = accumarray(index_mapping, dataPBS{1,4}.KPA', [], @(x) x(1))';

% Generate a vector of volume values within the range of unique_dataPBS_V, with increments of dV
% (Garuntees only injecting data used, no time series relaxation)
xv = round((min(unique_dataPBS_V):dV:max(unique_dataPBS_V)),2)';

% Initialize cell arrays to store calibrated pressure, truncated volume, and truncated pressure for PDMS
pressureSample_TC = cell(size(dataSample, 1), 1);
nLSample_T = cell(size(dataSample, 1), 1);
pressureSample_T = cell(size(dataSample, 1), 1);
timeSample_T = cell(size(dataSample, 1), 1);
pressureSample_FC = cell(size(dataSample, 1), 1);
nLSample_F = cell(size(dataSample, 1), 1);
timeSample_F = cell(size(dataSample, 1), 1);

% Generate a colormap for plotting with different colors
%colors = jet(size(dataPDMS, 1));
colors = customcolor(size(dataSample,1));

% % % Create a figure and hold it for adding plots
% % figure(2); hold on;

% Iterate through each row of dataSample
for idx = 1:size(dataSample, 1)
    %% Calibrate injecting data
    % Extract and round off VOLUME_TOTAL values for the current test sample data row
    Sample_Vol = round(dataSample{idx,4}.VOLUME_TOTAL,2);
    % Extract pressure values for the current test sample data row
    Sample_pressure = dataSample{idx,4}.KPA;
    Sample_time=dataSample{idx,4}.TOTAL_TEST_TIME;
    
    % Define the range for volume values based on xv (data truncated to
    % only injecting in PBS unique values)
    minX = min(xv(:));
    maxX = max(xv(:));
    
    % Identify indices where Test Sample volume falls within the defined range
    ind = find(Sample_Vol >= minX & Sample_Vol <= maxX);
    
    % Truncate the test sample volume and pressure values based on identified indices
    Sample_volTruncated = Sample_Vol(ind);
    Sample_pressureTruncated = Sample_pressure(ind);
    Sample_timeTruncated=Sample_time(ind);
    
    % Find the positions of PDMS_volTruncated in xv
    [~, loc_x] = ismember(Sample_volTruncated, xv);
    
    % Extract unique VOLUME_TOTAL values from dataPBS for the second column
    [uniqueVol_PBS, ~, index_mapping] = unique(dataPBS{1,4}.VOLUME_TOTAL);
    % Using the index mapping, create a new array of pressures for this column of dataPBS
    uniquePressure_PBS = accumarray(index_mapping, dataPBS{1,4}.KPA', [], @(x) x(1))';
    
    % Interpolate the uniquePressure_PBS values over the range of xv
    xp = interp1(uniqueVol_PBS,uniquePressure_PBS,xv);
    p_offset_points=10;
    xp_offset=xp-mean(xp(1:p_offset_points),'omitnan');
    % data validation
    deltaPBal=range(xp(1:p_offset_points));
    deltaPmax=.02;
    if deltaPBal>deltaPmax
        fprintf(['Pressure during balance step is not stable (>',...
            num2str(deltaPmax),'kPa). Verify your initial region for pressure balancing the calibration step\n']);
    end
    
    % Calibrate the PDMS pressure by subtracting the interpolated water pressures
    Sample_pressure_TC = Sample_pressureTruncated - mean(xp_offset(loc_x));
    
    %% Append time series data from before and after truncation
    % Extract and round off VOLUME_TOTAL values for the current test sample data row
    Sample_VolAll = round(dataSample{idx,3}.VOLUME_TOTAL,2);
    % Extract pressure values for the current test sample data row
    Sample_pressureAll = dataSample{idx,3}.KPA;
    Sample_timeAll=dataSample{idx,3}.TOTAL_TEST_TIME;
    
    % Find maximum volume of truncated data in the untruncated but filtered dataset
    trunc_idxhigh=find((Sample_timeAll>max(Sample_timeTruncated)),1);
        trunc_idxlow=find((Sample_timeAll<Sample_timeTruncated(1)),1,'last');

    % Append the truncated data to the later data for full pressure, volume, and time series
    Sample_pressure_FC=[Sample_pressureAll(1:trunc_idxlow); Sample_pressure_TC; Sample_pressureAll(trunc_idxhigh:end)];
    Sample_vol_F=[Sample_VolAll(1:trunc_idxlow);Sample_volTruncated; Sample_VolAll(trunc_idxhigh:end)];
    Sample_time_F=[Sample_timeAll(1:trunc_idxlow);Sample_timeTruncated;Sample_timeAll(trunc_idxhigh:end)];
   
    
% %     %% Plot
% %     % Plot the calibrated and original pressure against the truncated volume with different line widths
% %     figure(2)
% %     hold on
% %     if idx==1
% %         plot(radius(Sample_volTruncated),xp_offset(loc_x),'Color', [0 0 0], 'LineWidth',1,'DisplayName','calibration curve')
% %     end
% %     plot(radius(Sample_volTruncated), Sample_pressure_TC,'Color', colors(idx,:), 'LineWidth',2,'DisplayName','corrected data');
% %     plot(radius(Sample_volTruncated), Sample_pressureTruncated,'--','Color', colors(idx,:), 'LineWidth',1.0,'DisplayName','raw data');
% % 
% %     legend('Location','bestoutside')
% %     grid on
% %     title('Calibration corrected vs radius')
% %     ylabel('Pressure [kPa]')
% %     xlabel('Radius [mm]')
% %     
% %     figure(3)
% %     hold on
% %     plot(Sample_timeTruncated, Sample_pressureTruncated,'--','Color', colors(idx,:), 'LineWidth',1,'DisplayName','raw data');
% %     plot(Sample_time_F, Sample_pressure_FC,'Color', colors(idx,:), 'LineWidth',2.0,'DisplayName','corrected data');
% %     grid on
% %     title('Calibration corrected vs time')
% %     ylabel('Pressure [kPa]')
% %     xlabel('Time [s]')
    
    %% Store for export from function
    % Store the calibrated pressure, truncated volume, and truncated pressure in their respective cell arrays
    pressureSample_TC{idx} = Sample_pressure_TC;
    nLSample_T{idx} = Sample_volTruncated;
    pressureSample_T{idx} = Sample_pressureTruncated;
    timeSample_T{idx}=Sample_timeTruncated;
    
    % Store full datasets also
    nLSample_F{idx} = Sample_vol_F;
    pressureSample_FC{idx} = Sample_pressure_FC;
    timeSample_F{idx }= Sample_time_F;
    
    % Compute the accuracy/error in volume data
    % valueAccuracy = xv(loc_x) - PDMS_volTruncated;
end

% % % Release the figure hold
% % hold off;
end
function rad=radius(nL)
rad=((nL./1000).*(3./(4.*pi()))).^(1/3);
end
