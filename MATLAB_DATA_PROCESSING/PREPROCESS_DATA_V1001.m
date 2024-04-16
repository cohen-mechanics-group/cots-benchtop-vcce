% Brendan M Unikewicz, PhD Student
% Andre M Pincot, MSc Student
% Tal Cohen, Asc. Professor
% MIT, Dept. Mechanical Engineering
% MIT, Dept. Civil & Environmental Engineering
% Date of Creation: 03/21/2024
% Code Purpose: Pre-processing 43-, 45-, & 47-1 PDMS samples for properties

%% Code Start
clear; clc; format longg; close all; addpath('subFunctions');
cd ..\TEST_DATA\CSV_FILES\2024-03-21\  

%% Define Color Scheme and fontSize
redColor = [0.98, 0.40, 0.35];
greenColor = [0.45, 0.80, 0.69]; 
blueColor = [0.55, 0.60, 0.79];

fontSize = 14; 

%% Loading Files via Directories
% PDMS Raw CSV Data from Visual Studio Code (Foldered)
PDMS_43_1 = [pwd,'\43_1_010uL'];
PDMS_45_1 = [pwd,'\45_1_010uL'];
PDMS_47_1 = [pwd,'\47_1_010uL'];
PDMS_43_Wo = [pwd,'\43_1_250uL'];

cd ..\2024-04-12\

% PBS Raw CSV Data from Visual Studio Code (Filed)
PBS_010_CV = [pwd,'\PBS_CAL_300NLPS_010uL\PBS_calibration_300nlps_010uL_V1002.csv']; 
PBS_250_CV = [pwd,'\PBS_CAL_300NLPS_250uL\PBS_calibration_300nlps_250uL_V1002.csv']; 

cd ..\..\..\MATLAB_DATA_PROCESSING

%% Extracting Pressure-Volume Data from boolean'd files
fs = 500; % digitally upcycle sampling rate to uniform 500hz from ~80hz

[PDMS_43_1, dataPBS_43_1] = extractData_wPlot(fs, PDMS_43_1, PBS_010_CV); 
[PDMS_45_1, dataPBS_45_1] = extractData(fs, PDMS_45_1, PBS_010_CV); 
[PDMS_47_1, dataPBS_47_1] = extractData(fs, PDMS_47_1, PBS_010_CV); 
[PDMS_43_Wo, dataPBS_43_Wo] = extractData(fs, PDMS_43_Wo, PBS_250_CV); 

% Truncating and calibrating PDMS data from PBS result(s) 
[p_43_1, vPDMS_T_43_1, pPDMS_T_43_1, uV_PBS_43_1, uP_PBS_43_1,xp_off,xp_mean] = truncate_calibrateData_meanDyn(PDMS_43_1, dataPBS_43_1);  
[p_45_1, vPDMS_T_45_1, pPDMS_T_45_1, uV_PBS_45_1, uP_PBS_45_1] = truncate_calibrateData_meanDyn(PDMS_45_1, dataPBS_45_1);  
[p_47_1, vPDMS_T_47_1, pPDMS_T_47_1, uV_PBS_47_1, uP_PBS_47_1] = truncate_calibrateData_meanDyn(PDMS_47_1, dataPBS_47_1);  
[p_43_Wo, vPDMS_T_43_Wo, pPDMS_T_43_Wo, uV_PBS_43_Wo, uP_PBS_43_Wo] = truncate_calibrateData_meanDyn(PDMS_43_Wo, dataPBS_43_Wo);   

%% Fluid Movement Calibration (300nlps - 10uL syringe) 
mean_1 = repmat(mean(dataPBS_43_1{1}.KPA(1:886)),length(dataPBS_43_1{1}.KPA(1:886)),1);
mean_2 = repmat(mean(dataPBS_43_1{1}.KPA(887:3022)),length(dataPBS_43_1{1}.KPA(887:3022)),1);
mean_3 = repmat(mean(dataPBS_43_1{1}.KPA(3023:6839)),length(dataPBS_43_1{1}.KPA(3023:6839)),1);

figure; 
plot(dataPBS_43_1{1}.TOTAL_TEST_TIME(1:887),dataPBS_43_1{1}.KPA(1:887),'Color',greenColor,'LineWidth',1.5); hold on; 
plot(dataPBS_43_1{1}.TOTAL_TEST_TIME(887:3022),dataPBS_43_1{1}.KPA(887:3022),'Color',blueColor,'LineWidth',1.5); hold on; 
plot(dataPBS_43_1{1}.TOTAL_TEST_TIME(1:886),mean_1,'Color',redColor,'LineWidth',1.5); hold on; 
plot(dataPBS_43_1{1}.TOTAL_TEST_TIME(3022:6839),dataPBS_43_1{1}.KPA(3022:6839),'Color',greenColor,'LineWidth',1.5); hold on; 
plot(dataPBS_43_1{1}.TOTAL_TEST_TIME(887:3022),mean_2,'Color',redColor,'LineWidth',1.5); hold on; 
plot(dataPBS_43_1{1}.TOTAL_TEST_TIME(3023:6839),mean_3,'Color',redColor,'LineWidth',1.5); hold on;
xlabel('Time (s)','Interpreter','latex','FontSize',fontSize); ylabel('Pressure (kPa)','Interpreter','latex','FontSize',fontSize); 
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
legend('Static Pressure', 'Static \& Dynamic Pressure','Mean Pressure','Interpreter','latex','FontSize',fontSize-1,'Location','southeast');
xlim([0 57]);

%% Identifying pointPairs, their location, and identifying the initial defect for later fitting
desiredPointPair_locations = 9; 

% 43:1 PDMS Sample
[maxSlope_PDMS_43_1, pointPairs_PDMS_43_1, aBar_43_1, peakPressure_val_43_1, ... 
    peakPressure_idx_43_1, aBar_split_43_1, aBar_idx_43_1, aBar_regionFitting_43_1, ... 
    aBar_totalFitting_43_1,pressure_regionFitting_43_1,pressure_totalFitting_43_1] = ... 
    pointPairs_initialDefect_43_1(desiredPointPair_locations, vPDMS_T_43_1, p_43_1); 

% 45:1 PDMS Sample
[maxSlope_PDMS_45_1, pointPairs_PDMS_45_1, aBar_45_1, peakPressure_val_45_1, ... 
    peakPressure_idx_45_1, aBar_split_45_1, aBar_idx_45_1, aBar_regionFitting_45_1, ... 
    aBar_totalFitting_45_1,pressure_regionFitting_45_1,pressure_totalFitting_45_1] = ... 
    pointPairs_initialDefect_45_1(desiredPointPair_locations, vPDMS_T_45_1, p_45_1); 

% 47:1 PDMS Sample
[maxSlope_PDMS_47_1, pointPairs_PDMS_47_1, aBar_47_1, peakPressure_val_47_1, ... 
    peakPressure_idx_47_1, aBar_split_47_1, aBar_idx_47_1, aBar_regionFitting_47_1, ... 
    aBar_totalFitting_47_1,pressure_regionFitting_47_1,pressure_totalFitting_47_1] = ... 
    pointPairs_initialDefect_47_1(desiredPointPair_locations, vPDMS_T_47_1, p_47_1); 

% 43:1 PDMS Sample -- 250uL
[maxSlope_PDMS_43_Wo, pointPairs_PDMS_43_Wo, aBar_43_Wo, peakPressure_val_43_Wo, ... 
    peakPressure_idx_43_Wo, aBar_split_43_Wo, aBar_idx_43_Wo, aBar_regionFitting_43_Wo, ... 
    aBar_totalFitting_43_Wo,pressure_regionFitting_43_Wo,pressure_totalFitting_43_Wo] = ... 
    pointPairs_initialDefect_43_Wo(desiredPointPair_locations, vPDMS_T_43_Wo, p_43_Wo); 

%% Fitting the NeoHookean Model
elasticGuess = 50; % Defining an Initial Guess for Fitting
percentile = 99; % Defining % for slope-calc

%43:1 PDMS Sample
[fitParams_43_1, fullFit_data_43_1, fullFit_dataY_43_1, fitParams_optimized_43_1] = neoHookean_fitting(percentile,pointPairs_PDMS_43_1,elasticGuess, ...
    aBar_regionFitting_43_1,pressure_regionFitting_43_1,aBar_totalFitting_43_1,aBar_43_1);

%45:1 PDMS Sample
[fitParams_45_1, fullFit_data_45_1, fullFit_dataY_45_1, fitParams_optimized_45_1] = neoHookean_fitting(percentile,pointPairs_PDMS_45_1,elasticGuess, ...
    aBar_regionFitting_45_1,pressure_regionFitting_45_1,aBar_totalFitting_45_1,aBar_45_1);

%47:1 PDMS Sample
[fitParams_47_1, fullFit_data_47_1, fullFit_dataY_47_1, fitParams_optimized_47_1] = neoHookean_fitting(percentile,pointPairs_PDMS_47_1,elasticGuess, ...
    aBar_regionFitting_47_1,pressure_regionFitting_47_1,aBar_totalFitting_47_1,aBar_47_1);

%43:1 PDMS Sample - 250uL
[fitParams_43_Wo, fullFit_data_43_Wo, fullFit_dataY_43_Wo, fitParams_optimized_43_Wo] = neoHookean_fitting(percentile,pointPairs_PDMS_43_Wo,elasticGuess, ...
    aBar_regionFitting_43_Wo,pressure_regionFitting_43_Wo,aBar_totalFitting_43_Wo,aBar_43_Wo);

%% Pulling Revelant Data (Initial Defect, Elastic Modulus, aBar, Stretch, pPDMS_TC_43_1, ...)
%43:1 PDMS Sample
for i=1:length(maxSlope_PDMS_43_1)
    initialDefect_43_1(1,i) = maxSlope_PDMS_43_1{i}.x_intercept; 
    fields_43_1 = fieldnames(fitParams_43_1{i}.pointPairData); % Get all field names
    pointPairField_43_1 = fields_43_1{contains(fields_43_1, 'pointPair')}; % Find fields that contain 'pointPair'
    E_i_43_1(1,i) = fitParams_43_1{i}.pointPairData.(pointPairField_43_1)(2); 
    stretch_43_1{i} = aBar_43_1{i} / initialDefect_43_1(1,i); 
    fittingStretch_43_1{i} = aBar_regionFitting_43_1{i} / initialDefect_43_1(1,i);
end

%45:1 PDMS Sample
for i=1:length(maxSlope_PDMS_45_1)
    initialDefect_45_1(1,i) = maxSlope_PDMS_45_1{i}.x_intercept; 
    fields_45_1 = fieldnames(fitParams_45_1{i}.pointPairData); % Get all field names
    pointPairField_45_1 = fields_45_1{contains(fields_45_1, 'pointPair')}; % Find fields that contain 'pointPair'
    E_i_45_1(1,i) = fitParams_45_1{i}.pointPairData.(pointPairField_45_1)(2); 
    stretch_45_1{i} = aBar_45_1{i} / initialDefect_45_1(1,i); 
    fittingStretch_45_1{i} = aBar_regionFitting_45_1{i} / initialDefect_45_1(1,i);
end

%47:1 PDMS Sample
for i=1:length(maxSlope_PDMS_47_1)
    initialDefect_47_1(1,i) = maxSlope_PDMS_47_1{i}.x_intercept; 
    fields_47_1 = fieldnames(fitParams_47_1{i}.pointPairData); % Get all field names
    pointPairField_47_1 = fields_47_1{contains(fields_47_1, 'pointPair')}; % Find fields that contain 'pointPair'
    E_i_47_1(1,i) = fitParams_47_1{i}.pointPairData.(pointPairField_47_1)(2);
    stretch_47_1{i} = aBar_47_1{i} / initialDefect_47_1(1,i); 
    fittingStretch_47_1{i} = aBar_regionFitting_47_1{i} / initialDefect_47_1(1,i);
end

%43:1 PDMS Sample -- Womersley Validation
for i=1:length(maxSlope_PDMS_43_Wo)
    initialDefect_43_Wo(1,i) = maxSlope_PDMS_43_Wo{i}.x_intercept; 
    fields_43_Wo = fieldnames(fitParams_43_Wo{i}.pointPairData); % Get all field names
    pointPairField_43_Wo = fields_43_Wo{contains(fields_43_Wo, 'pointPair')}; % Find fields that contain 'pointPair'
    E_i_43_Wo(1,i) = fitParams_43_Wo{i}.pointPairData.(pointPairField_43_Wo)(2); 
    stretch_43_Wo{i} = aBar_43_Wo{i} / initialDefect_43_Wo(1,i); 
    fittingStretch_43_Wo{i} = aBar_regionFitting_43_Wo{i} / initialDefect_43_Wo(1,i);
end

%% Defining Relevant Data Variables
varList = {'initialDefect_43_1','initialDefect_45_1','initialDefect_47_1','initialDefect_43_Wo', ...
    'aBar_43_1','aBar_45_1','aBar_47_1','aBar_43_Wo', ...
    'aBar_regionFitting_43_1','aBar_regionFitting_45_1','aBar_regionFitting_47_1','aBar_regionFitting_43_Wo', ...
    'p_43_1', 'p_45_1', 'p_47_1', 'p_43_Wo', ... 
    'pressure_regionFitting_43_1', 'pressure_regionFitting_45_1', 'pressure_regionFitting_47_1', 'pressure_regionFitting_43_Wo', ...
    'stretch_43_1', 'stretch_45_1', 'stretch_47_1', 'stretch_43_Wo', ...
    'fittingStretch_43_1', 'fittingStretch_45_1', 'fittingStretch_47_1', 'fittingStretch_43_Wo', ...
    'fullFit_dataY_43_1', 'fullFit_dataY_45_1', 'fullFit_dataY_47_1', 'fullFit_dataY_43_Wo', ... 
    'initialDefect_43_1', 'initialDefect_45_1', 'initialDefect_47_1', 'initialDefect_43_Wo', ...
    'E_i_43_1', 'E_i_45_1', 'E_i_47_1', 'E_i_43_Wo'};

folderName = 'data_preProcessing';

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
