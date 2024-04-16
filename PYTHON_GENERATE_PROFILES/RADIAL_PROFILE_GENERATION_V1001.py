# Brendan M Unikewicz, PhD Student
# MIT, Dept. Mechanical Engineering
# Date of Creation: 06162023
# Code Purpose: Generating profiles for microPump, i.e. controlPump.py

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# User Inputs
aStart_cycling = float(input('Enter starting value for a_start (mm): '))
aStop_cycling = float(input('Enter stopping value for a_stop (mm): '))
nInterpolating = int(input('Enter the number of interpolating points (N): '))

ratesInfusion = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32]  # mm/s

samplingRate = 10000  # Hz, this is non-important... just DSP-chosen value
dt = 1 / samplingRate

mm3_nL = 1000  # mm^3 to nL

fig, ax = plt.subplots(2, 1)
radialInfusion_profiles = []
timeVectors = []
volumeInfusion_profiles = []
timeInfuse_interpolateds = []
volumeInfuse_interpolateds = []
volInfuse_rates = []
volInfuse_mags = []
interpInfuses = []
actualInfuses = []
pDiff_infuses = []
infusionTimes = []

# Infusion Rates Generation
for rate in ratesInfusion:
    radialInfusion_profile = np.arange(aStart_cycling, aStop_cycling, dt*rate)
    timeVector = np.arange(1, len(radialInfusion_profile)+1) * dt
    volumeInfusion_profile = (4/3) * np.pi * radialInfusion_profile**3 * mm3_nL

    timeInfuse_interpolated = np.linspace(timeVector[0], timeVector[-1], nInterpolating)
    volumeInfuse_interpolated = np.interp(timeInfuse_interpolated, timeVector, volumeInfusion_profile)

    ax[1].plot(timeInfuse_interpolated, volumeInfuse_interpolated, 'ok', label=f'Rate: {rate} mm/s')

    volInfuse_rate = np.diff(volumeInfuse_interpolated) / np.diff(timeInfuse_interpolated)
    volInfuse_mag = np.diff(volumeInfuse_interpolated)

    interpInfuse = np.trapz(volumeInfuse_interpolated, timeInfuse_interpolated)
    actualInfuse = np.trapz(volumeInfusion_profile, timeVector)
    pDiff_infuse = ((interpInfuse - actualInfuse) / interpInfuse) * 100

    ax[0].plot(timeVector, radialInfusion_profile, label=f'Rate: {rate} mm/s')
    ax[1].plot(timeVector, volumeInfusion_profile, label=f'Rate: {rate} mm/s')

    infusionTime = timeVector[-1]

    radialInfusion_profiles.append(radialInfusion_profile)
    timeVectors.append(timeVector)
    volumeInfusion_profiles.append(volumeInfusion_profile)
    timeInfuse_interpolateds.append(timeInfuse_interpolated)
    volumeInfuse_interpolateds.append(volumeInfuse_interpolated)
    volInfuse_rates.append(volInfuse_rate)
    volInfuse_mags.append(volInfuse_mag)
    interpInfuses.append(interpInfuse)
    actualInfuses.append(actualInfuse)
    pDiff_infuses.append(pDiff_infuse)
    infusionTimes.append(infusionTime)

ax[0].grid(True, which='both', linestyle='--', linewidth=0.5)
ax[1].grid(True, which='both', linestyle='--', linewidth=0.5)
ax[0].set_title('Infusion Radial (top) & Volumetric (bottom) Rates')
plt.show()

# Writing Infusion data to CSV
pd.DataFrame(volInfuse_mags).to_csv('infusionMags.csv', index=False, header=False)
pd.DataFrame(volInfuse_rates).to_csv('infusionRates.csv', index=False, header=False)
pd.DataFrame(infusionTimes).to_csv('infusionTimes.csv', index=False, header=False)
