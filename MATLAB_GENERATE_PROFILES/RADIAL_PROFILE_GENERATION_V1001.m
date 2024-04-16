% Brendan M Unikewicz, PhD Student
% MIT, Dept. Mechanical Engineering
% Date of Creation: 06162023
% Code Purpose: Generating profiles for microPump, i.e. controlPump.py

clear; close all; format bank; 

% User Inputs
aStart_cycling = input('Enter starting value for a_start (mm): ');
aStop_cycling = input('Enter stopping value for a_stop (mm): ');
nInterpolating = input('Enter the number of interpolating points (N): ');

ratesInfusion = [0.01 0.02 0.04 0.08 0.16 0.32]; % mm/s

%10uLSyringe_maxRate = 657; % nL/s
%25uLSyringe_maxRate = 1482; % nL/s
%50uLSyringe_maxRate = 2963; % nL/s

samplingRate = 500; % Hz, this is non-important... just DSP-chosen value
dt = 1/samplingRate; 

mm3_nL = 1000; % mm^3 to nL

figure; 
% Infusion Rates Generation
for i = 1:length(ratesInfusion)
    radialInfusion_profile{i,:} = aStart_cycling:(dt*ratesInfusion(i)):aStop_cycling; 
    timeVector{i,:} = (1:1:length(radialInfusion_profile{i,:}))*dt; 
    volumeInfusion_profile{i,:} = (4/3) * pi * (radialInfusion_profile{i,1}).^3 * mm3_nL; 
                         
    timeInfuse_interpolated{i,:} = linspace(min(timeVector{i,:}),max(timeVector{i,:}),nInterpolating);
    volumeInfuse_interpolated{i,:} = interp1(timeVector{i,:},volumeInfusion_profile{i,:},timeInfuse_interpolated{i,:}) ;
    subplot(2,1,2)
    plot(timeInfuse_interpolated{i,:},volumeInfuse_interpolated{i,:},'ok'); hold on; 

    volInfuse_rate{i,:} = diff(volumeInfuse_interpolated{i,:}) ./ diff(timeInfuse_interpolated{i,:}); 
    volInfuse_mag{i,:} = diff(volumeInfuse_interpolated{i,:}); 

    % Integration check
    interpInfuse(i) = trapz(timeInfuse_interpolated{i,:},volumeInfuse_interpolated{i,:});
    actualInfuse(i) = trapz(timeVector{i,:},volumeInfusion_profile{i,:});
    pDiff_infuse(i) = ((interpInfuse(i) - actualInfuse(i)) / interpInfuse(i)) * 100; 

    subplot(2,1,1); title('Infusion Radial (top) & Volumetric (bottom) Rates'); 
    plot(timeVector{i,:},radialInfusion_profile{i,:},'LineWidth',1.5); hold on; 
    ylabel('Effective Cavity Radii (mm)'); xlabel('Time (s)')
    grid on; grid minor; % legend(['Rate: ',num2str(ratesInfusion(i)),'mm/s']); hold on; 
    subplot(2,1,2); 
    plot(timeVector{i,:},volumeInfusion_profile{i,:},'LineWidth',1.5); hold on; 
    ylabel('Infused Fluid (nL)'); xlabel('Time (s)')
    grid on; grid minor; % legend(['Rate: ',num2str(ratesInfusion(i)),'mm/s']); hold on;

    infusionTime(i) = max(timeVector{i,:}); 
end
hold off; 

% Writing Infusion data to CSV
filenameInfuse_rate = 'infusionRates.csv'; 
filenameInfuse_mag = 'infusionMags.csv'; 
csvwrite(filenameInfuse_mag, cell2mat(volInfuse_mag));
csvwrite(filenameInfuse_rate, cell2mat(volInfuse_rate));

% Writing Infusion time data to CSV
filenameInfuse_time = 'infusionTimes.csv'; 
csvwrite(filenameInfuse_time, infusionTime);
