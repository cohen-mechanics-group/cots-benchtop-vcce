% Brendan M Unikewicz, PhD Student
% Andre M Pincot, MSc Student
% Tal Cohen, Asc. Professor
% MIT, Dept. Mechanical Engineering
% MIT, Dept. Civil & Environmental Engineering
% Component Analysis: worse-case scenarios
% Date of Creation: 03/25/2024
% Code Purpose: checking structural components for volumetric losses

%% Code Start
clear; close all; clc; 

%% Define Colors
redColor = [0.98, 0.40, 0.35];
greenColor = [0.45, 0.80, 0.69]; 
blueColor = [0.55, 0.60, 0.79];

c = [redColor; greenColor; blueColor];

fontSize = 14; 
%% Constants
g = 9.81; 
E_aluminum = 70e9; % aluminum modulus of elasticity
peakPressure = [5e3,50e3,100e3]; % peak pressures 

m3_2_nL = 1.0E12; % convert m^3 to nL
volPeak_avg_all = 377.26; 

%% Needle Characteristics
needle_D = 0.000508; % outer diameter, m
needleArea = pi*(needle_D/2)^2; % m^2

%% Critical Load for Buckling: Syringe Plunger
radiusPlunger = (0.485e-3 / 2); % radius of the plunger
I_plunger = (pi/4)*(radiusPlunger)^4; % second moment of area of the plunger
areaPlunger = pi*radiusPlunger^2; % cross-sectional area of the plunger

forcePressure = peakPressure * areaPlunger; 
forceFriction = 1; % Newtons -- experimentally from Roth thesis
 
forceTotal = forceFriction + forcePressure; %+ forceViscousFluids; 
K = 1.2; % design factor -- 1 vs 1.2

bucklingLength_theory = sqrt((((pi^2)*(E_aluminum*I_plunger))./(K^2.*forceTotal./5))); 
exposedLength =  65e-3; 
safetyFactor_buckling = bucklingLength_theory / exposedLength; 

%% Axially Loading Plunger Displacement
totalPlunger_length = 90e-3; % 65mm?

deltaPlunger = (forcePressure .* totalPlunger_length) ./ (areaPlunger*E_aluminum); 
dV_plunger_nL = deltaPlunger.*areaPlunger.*m3_2_nL; % change i  
dV_plunger_percentage = dV_plunger_nL ./ volPeak_avg_all * 100; 

%% Compressibility of Water in Syringe
densityDistilled = 1000; % density of distilled water, kg/m^3
densityPBS = 1005; 
densitySolution = (densityDistilled*9 + densityPBS*1) / 10; 

Vo_WATER = 4e-4 / densitySolution; % 0.4 grams to fill the syringe -- 1kg = 1L
K_WATER = 2.22e9; % bulk modulus, water
dV_WATER_m3 = (peakPressure .* (Vo_WATER)) ./ (-K_WATER);
dV_WATER_nL = dV_WATER_m3 .* m3_2_nL;
dV_WATER_percentage = dV_WATER_nL ./ volPeak_avg_all * 100; 

%% Expansion of Glass Syringe Body (10uL)
nu = 0.2; % Poisson's ratio
E = 60e9; % Young's modulus
Pi = [5e3, 50e3, 100e3]; % Internal pressure
Po = 0; 
ri = (0.485e-3)/2; % Inner radius
ro = (6.60e-3)/2; % Outer radius
%r = ri; % Radial distance at which deformation is calculated
h_syringe = 0.1; 

% Calculate the radial deformation, Delta r
delta_ri = ((1 - nu) / E) .* ((Pi .* ri^2 - Po .* ro^2) ./ (ro^2 - ri^2)) .* ri + ...
          ((1 + nu) / E) * ((Pi - Po) .* ri^2 .* ro^2 ./ ((ro^2 - ri^2) * ri));

% Calculate the radial deformation, Delta r
delta_ro = ((1 - nu) / E) .* ((Pi .* ri^2 - Po .* ro^2) ./ (ro^2 - ri^2)) .* ro + ...
          ((1 + nu) / E) .* ((Pi - Po) .* ri^2 .* ro^2 ./ ((ro^2 - ri^2) * ro));

dV_syringe = pi.*h_syringe.*(delta_ro.^2 + 2.*delta_ro.*ro - 2.*delta_ri.*ri - delta_ri.^2); 
dV_syringe_nL = m3_2_nL .* dV_syringe; 
dV_syringe_percentage = dV_syringe_nL ./ volPeak_avg_all * 100; 

%% Expansion of Luer Needle (25G Needle Barrel)
nu_25G = 0.29; % Poisson's ratio
E_25G = 193e9; % Young's modulus
ri_25G = (0.3048e-3)/2; % Inner radius
ro_25G = (0.508e-3)/2; % Outer radius
h_25G = 0.0254; % 1in

% Calculate the radial deformation, Delta r
delta_ri_25G = ((1 - nu_25G) / E_25G) .* ((Pi .* ri_25G^2 - Po .* ro_25G^2) ./ (ro_25G^2 - ri_25G^2)) .* ri_25G + ...
          ((1 + nu_25G) / E_25G) .* ((Pi - Po) .* ri_25G^2 * ro_25G^2 ./ ((ro_25G^2 - ri_25G^2) * ri_25G));

% Calculate the radial deformation, Delta r
delta_ro_25G = ((1 - nu_25G) / E_25G) .* ((Pi .* ri_25G^2 - Po .* ro_25G^2) ./ (ro_25G^2 - ri_25G^2)) .* ro_25G + ...
          ((1 + nu_25G) / E_25G) .* ((Pi - Po) .* ri_25G^2 .* ro_25G^2 ./ ((ro_25G^2 - ri_25G^2) * ro_25G));

dV_25G = pi.*h_25G.*(delta_ro_25G.^2 + 2.*delta_ro_25G.*ro_25G - 2.*delta_ri_25G.*ri_25G - delta_ri_25G.^2); 
dV_25G_nL = m3_2_nL .* dV_25G; 
dV_25G_percentage = dV_25G_nL ./ volPeak_avg_all * 100; 

%% Expansion of Pressure Sensor Channel -- Polycarbonate
nu = 0.35; % Poisson's ratio
E = 2.4e9; % Young's modulus - polycarbonate
ri = (2.9e-3)/2; % Inner radius; dia=2.9mm
ro = (4.315e-3)/2; % Outer radius; dia=4.315mm
h_pressureChannel = 0.00635 * 4; % 1in length

% Calculate the radial deformation, Delta r
delta_ri = ((1 - nu) / E) .* ((Pi .* ri^2 - Po .* ro^2) ./ (ro^2 - ri^2)) .* ri + ...
          ((1 + nu) / E) * (((Pi - Po) .* ri^2 .* ro^2) ./ ((ro^2 - ri^2) * ri));

% Calculate the radial deformation, Delta r
delta_ro = ((1 - nu) / E) .* ((Pi .* ri^2 - Po .* ro^2) ./ (ro^2 - ri^2)) .* ro + ...
          ((1 + nu) / E) .* (((Pi - Po) .* ri^2 .* ro^2) ./ ((ro^2 - ri^2) * ro));

dV_pressureChannel = pi.*h_pressureChannel.*(delta_ro.^2 + 2.*delta_ro.*ro - 2.*delta_ri.*ri - delta_ri.^2); % Volume loss in m^3

dV_pressureChannel_nL = m3_2_nL .* dV_pressureChannel; % volume loss in nL
dV_pressureChannel_percentage = dV_pressureChannel_nL ./ volPeak_avg_all * 100; 

%% Expansion of Luer Needle Polypropolene (Needle Cap)
nu_PP = 0.42; % Poisson's ratio
E_PP = 1.5e9; % Young's modulus
ri_PP = (0.0061976)/2; % Inner radius
ro_PP = (0.0046228)/2; % Outer radius
h_PP = 0.003175; % 1/8in

% Calculate the radial deformation, Delta r
delta_ri_PP = ((1 - nu_PP) / E_PP) .* ((Pi .* ri_PP^2 - Po .* ro_PP^2) ./ (ro_PP^2 - ri_PP^2)) .* ri_PP + ...
          ((1 + nu_PP) / E_PP) .* ((Pi - Po) .* ri_PP^2 .* ro_PP^2 ./ ((ro_PP^2 - ri_PP^2) * ri_PP));

% Calculate the radial deformation, Delta r
delta_ro_PP = ((1 - nu_PP) / E_PP) .* ((Pi .* ri_PP^2 - Po .* ro_PP^2) ./ (ro_PP^2 - ri_PP^2)) .* ro_PP + ...
          ((1 + nu_PP) / E_PP) .* ((Pi - Po) .* ri_PP^2 .* ro_PP^2 ./ ((ro_PP^2 - ri_PP^2) * ro_PP));

dV_PP = pi.*h_PP.*(delta_ro_PP.^2 + 2.*delta_ro_PP.*ro_PP - 2.*delta_ri_PP.*ri_PP - delta_ri_PP.^2); 
dV_PP_nL = m3_2_nL .* dV_PP; 
dV_PP_percentage = dV_PP_nL ./ volPeak_avg_all * 100; 

%% Compressibility of MEMS Diaphragm in Pressure Sensor
r_MEMS = 1e-3; 
h_MEMS = 1e-4; 
A_MEMS = pi*((r_MEMS)^2); 

E_MEMS = 5e6;
v_MEMS = 0.48; % poisson's ratio MEMS  

q = (peakPressure .* A_MEMS) ./ r_MEMS; 

aRadial_MEMS = 0:(r_MEMS/1000):r_MEMS; 
flexuralRigidity = (E_MEMS*h_MEMS^3) / (12*(1-v_MEMS^2)); 

w = -q' ./ (64.*flexuralRigidity).*(((r_MEMS^2)-(aRadial_MEMS).^2).^2); 

figure; 
plot((aRadial_MEMS*1e3),(w(1,:).*1e6),'Color',redColor,'LineWidth',1.5); hold on; 
plot((aRadial_MEMS*1e3),(w(2,:).*1e6),'Color',blueColor,'LineWidth',1.5); hold on; 
plot((aRadial_MEMS*1e3),(w(3,:).*1e6),'Color',greenColor,'LineWidth',1.5); hold on;  
xlabel('Distance from Center (mm)','Interpreter','latex','FontSize',fontSize); 
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 

ylabel('Displacement ($\mu$m)','Interpreter','latex','FontSize',fontSize);
legend('5kPa','50kPa','100kPa','Interpreter','latex','Location','southeast','FontSize',fontSize-1); hold off; 

dV_MEMS_m3 = (pi.*(r_MEMS.^6).*q) ./ (192.*flexuralRigidity); 
dV_MEMS_nL = dV_MEMS_m3 .* m3_2_nL; 
dV_MEMS_percentage = dV_MEMS_nL ./ volPeak_avg_all * 100; 

%% PendoTech Calibration Visualized
direc_PendoTech_cal_upper = [pwd,'\calibration_pendoTech\upperLimit.csv']; 
direc_PendoTech_cal_lower = [pwd,'\calibration_pendoTech\lowerLimit.csv']; 
dataCalibration_upper_MEMS = csvread(direc_PendoTech_cal_upper); 
dataCalibration_lower_MEMS = csvread(direc_PendoTech_cal_lower); 
calibratedPSI_upperLimit = dataCalibration_upper_MEMS(:,1); 
sensorPSI_upperLimit = dataCalibration_upper_MEMS(:,2); 
calibratedPSI_lowerLimit = dataCalibration_lower_MEMS(:,1); 
sensorPSI_lowerLimit = dataCalibration_lower_MEMS(:,2); 
theoryCalibrated_PSI = 0:1:50; 
theorySensor_PSI = 0:1:50; 

for i=1:length(theorySensor_PSI)
    if 0<i<7
        rejection_upperLimit = theorySensor_PSI * 1.02; 
        rejection_lowerLimit = theoryCalibrated_PSI * 0.98; 
    elseif 8<i<31
        rejection_upperLimit = theorySensor_PSI * 1.03; 
        rejection_lowerLimit = theoryCalibrated_PSI * 0.97; 
    else
        rejection_upperLimit = theorySensor_PSI * 1.05; 
        rejection_lowerLimit = theoryCalibrated_PSI * 0.95; 
    end
end
        
figure; 
plot(theoryCalibrated_PSI,theorySensor_PSI,'k','LineWidth',1.5); hold on;
plot(calibratedPSI_upperLimit,sensorPSI_upperLimit,'Color',redColor,'LineWidth',1.5); hold on; 
plot(theoryCalibrated_PSI,rejection_upperLimit,'Color',greenColor,'LineWidth',1.5); hold on; 
plot(calibratedPSI_lowerLimit,sensorPSI_lowerLimit,'Color',redColor,'LineWidth',1.5); hold on; 
plot(theoryCalibrated_PSI,rejection_lowerLimit,'Color',greenColor,'LineWidth',1.5); hold off; 
legend('Ideal Pressure Sensor Output','Manufacturing Specification','PendoTech Quality Standard','Interpreter','latex','Location','northwest','FontSize',fontSize-1); 
ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
xlabel('Calibrated Pressure (PSI)','Interpreter','latex','FontSize',fontSize); 
ylabel('Sensor Reading (PSI)','Interpreter','latex','FontSize',fontSize);  
xlim([0 50]); ylim([0 50]); 
