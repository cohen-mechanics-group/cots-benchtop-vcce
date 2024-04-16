#! /usr/bin/env python
##########################################
## Code Purpose: 
## 1. Test Soft Materials
##
## Contact Information:                            
## Massachusetts Institute of Technology
## Dept. of Mechanical Engineering
## Cohen's Nonlinear Solid Mechanics Group
##                                      
## Development: 
## Brendan M Unikewicz, PhD Student, MIT
## <bmu@mit.edu> 
##
## Please CC my advisor as well:
## Tal Cohen, Professor, MIT
## <talco@mit.edu>
##########################################

import serial 
import sys
import time
import glob
import csv
from os import name, system
import tkinter as tk
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import os
import hid 
import csv
import keyboard
import datetime 

# Define Environment Variable, Begin
os.environ["BLINKA_U2IF"] = "1"
import board 

# Open Device
device=hid.device()
device.open(0x239A, 0x0109)

from cedargrove_nau7802 import NAU7802 # downloaded zip for archival in GitHub

global zcom, LogDump, current_time

# Config: Command parameters
Config = ["F5", "?V", "?C", "?B", "?E", "?M", "?S", "?D", "?U", "?G"]

# ConfigOut1: CURRENT CONFIGURATION FOR PUMP #1
## ConfigOut1 = ["", "", "", "Motor Drive Option:",
##              "Motor Counter Mode:", "Pump ", "Syringe ", "", "", "", ""]

# chgConfig[10]# RESERVED: NOT USED #
chgConfig = [""]*10

#RunPSet: INITIALIZING A DEFAULT FOR PUMP #1
RunPSet = ["1", "INFUSE", "50", "50", "SEC", "0", "0"]

# LogDump: DUMP FILE FOR ERROR EXCEPTION
LogDump = [" "]*20
LogDump[0] = "**********************************************"

# Defining the 24-bit ADC I2C, address, numChannels, etc
numChannels = 1
nau7802 = NAU7802(board.I2C(), address=0x2A, active_channels=numChannels)
samplingRate = 80 # default set in cedargrove_nau7802
gain = 1
softWareDelay = []

def enableADC_temperatureSensor(nau7802, gain, samplingRate):
    nau7802.gain = gain
    # nau7802._c2_conv_rate
    nau7802._ts = 0x02
    enabled_tempSense = nau7802.enable(True)

    return enabled_tempSense

def enableADC_wheatBridge(nau7802, gain, samplingRate):
    nau7802.gain = gain
    #nau7802._c2_conv_rate = samplingRate
    nau7802._ts = 0x00
    enabled_wheatBridge = nau7802.enable(True)

    return enabled_wheatBridge
    
def zero_channel():
    """Initiate internal calibration for current channel.Use when scale is started,
    a new channel is selected, or to adjust for measurement drift. Remove weight
    and tare from load cell before executing."""
    #print(
    #    "channel %1d calibrate.INTERNAL: %5s"
    #    % (nau7802.channel, nau7802.calibrate("INTERNAL"))
    #)
    print(
        "channel %1d calibrate.OFFSET:   %5s"
        % (nau7802.channel, nau7802.calibrate("OFFSET"))
    )
    #print(
    #    "channel %1d calibrate.GAIN:   %5s"
    #    % (nau7802.channel, nau7802.calibrate("GAIN"))
    #)
    print("...channel %1d zeroed" % nau7802.channel)
# END DEF

# Read ADC and convert to temperature
def read_and_convert_temperature(sampling_rate, gain, num_channels):
    """Read and average raw sample values from both channels, convert them to temperature and return the values."""
    
    def read_raw_value(sampling_rate):
        """Read and average consecutive raw sample values. Return average raw value."""
        samples = 3 * sampling_rate  # number of samples is three times the sampling rate
        sample_sum = 0
        sample_count = samples
        while sample_count > 0:
            if nau7802.available:
                sample_sum = sample_sum + nau7802.read()
                sample_count -= 1
        return int(sample_sum / samples)

    def convert_to_temp(value, gain):
        V_ref = 3.0  # Reference voltage -- basically 2.972V
        nBits = 2**24  # 24-bit ADC, i.e. 2^24 bits
        adcVoltage = (value / nBits) * V_ref  # This will be in Volts
        # Considering gain
        measuredVoltage = adcVoltage / gain  # Adjust for gain, this is the real voltage
        # Convert to mV for easier calculation
        measuredVoltage_mV = measuredVoltage * 1000  # Convert voltage to mV
        # Values for the temperature sensor
        baseVoltage = 109  # Base voltage for 25C
        sensitivity = 0.36  # Sensitivity in mV per degree C
        # Calculate deviation from base voltage
        voltageDeviation = measuredVoltage_mV - baseVoltage
        # Calculate temperature deviation
        tempDeviation = voltageDeviation / sensitivity
        # Base temperature + temperature deviation
        temperature = 25 + tempDeviation
        return temperature

    temp_ch1, temp_ch2 = None, None

    for i in range(1, num_channels + 1):
        nau7802.channel = i
        value = read_raw_value(sampling_rate)
        temp = convert_to_temp(value, gain)
        if i == 1:
            temp_ch1 = temp
        elif i == 2:
            temp_ch2 = temp

    return temp_ch1, temp_ch2
# END DEF

# Initialize and record temperature measurement
enabled_tempSense = enableADC_temperatureSensor(nau7802, gain, samplingRate)
temp_ch1, temp_ch2 = read_and_convert_temperature(samplingRate, gain, numChannels)

# Initialize, zero, and calibrate wheatstoneBride arrangement
enabled_wheatBridge = enableADC_wheatBridge(nau7802, gain, samplingRate)
# zero_channel()

def clear_line():
    sys.stdout.write('\r')
    sys.stdout.write('\033[K')
    sys.stdout.write('\r')
    sys.stdout.flush()
# END DEF

def read_csv_file(filename, num_rows, header_row=False):
    # Open the CSV file for reading
    with open(filename, 'r') as csv_file:
        # Create a CSV reader object
        csv_reader = csv.reader(csv_file)

        # Read the header row if it exists and skip it
        if header_row:
            next(csv_reader)

        # Get the number of columns in the CSV file
        first_row_values = next(csv_reader)
        num_columns = len(first_row_values)

        # Create a list of lists to store the values for each row and column
        data = [[None for _ in range(num_columns)] for _ in range(num_rows)]

        # Store the values of the first row in the data list
        for j, value in enumerate(first_row_values):
            variable_name = f"row1_column{j+1}"
            data[0][j] = value

        # Loop through each remaining row in the CSV file
        for i in range(1, num_rows):
            # Extract the values from the current row and store them in separate variables
            row_values = next(csv_reader)
            # Loop through each value in the current row and store it in the corresponding variable
            for j, value in enumerate(row_values):
                variable_name = f"row{i+1}_column{j+1}"
                data[i][j] = value

    return data
# END DEF

# All VCCE Profile Information (Radial Growth Constant)
TIME_FORMAT='%Y-%m-%d %H:%M:%S'

numRows_profile = 6
profileRates = read_csv_file('MATLAB_GENERATE_PROFILES/infusionRates.csv', numRows_profile, header_row=False)
profileMags = read_csv_file('MATLAB_GENERATE_PROFILES/infusionMags.csv', numRows_profile, header_row=False)

profileTimes = read_csv_file('MATLAB_GENERATE_PROFILES/infusionTimes.csv',1,header_row=False)
profileTimes = [[float(item) for item in sublist] for sublist in profileTimes]

# Initialize Relaxation Time
relaxTime = 60

# Input the number of iterations
nIterations = len(profileRates[0])
nIterations_time = nIterations

# No Tocques Muchacho; need to clean up.
nIterations *= 150
nIterations_tSequence = nIterations_time * (nIterations / 1000)
nIterations = str(nIterations)

# Preload Rates for Calling Function/Profile Indices
rate001 = 0
rate002 = 1
rate004 = 2
rate008 = 3
rate016 = 4
rate032 = 5

# Variable Writing for TEST_LOGS
needleLength = 1
needleGauge = 25
pressureSensor = 1001

print("\n")
print("######################################################################")
print("###            24-BIT ADC HAS BEEN ZEROED AND CALIBRATED           ###")
print("###                                  &                             ###")
print("###        PROFILE INFORMATION HAS BEEN SUCCESSFULLY LOADED        ###")
print("######################################################################")
print("\n")
time.sleep(3)

def WriteLog(flag):
    global LogDump, current_time
    if(flag == 0):
        LogDump = [" "]*20
        LogDump[0] = "**********************************************"
    if(flag == 1):
        td = time.strftime("%Y-%b-%d %H:%M:%S", time.localtime())
        try:
            log = open("error.log", "a")
            log.write(LogDump[0]+"\n")
            log.write(td)
            for i in range(10):
                log.write(LogDump[i]+"\n")
            log.write(LogDump[0]+"\n")
            log.close()
        except:
            print(" #! ERROR LOG EXCEPTION !#")
    return
# END DEF

def ClrScrn():
    try:
        if name == 'posix':
            _ = system('clear')  # for mac and linux (os.name is'posix')
        else:
            _ = system('cls')  # for windows (os.name is typically 'nt')
    except:
        print("\n")
    return
# END DEF

def serial_ports():
    # Find local active serial ports
    if sys.platform.startswith('win'):
        ports = ['COM%s' % (i + 1) for i in range(256)]
    elif sys.platform.startswith('linux') or sys.platform.startswith('cygwin'):
        # this excludes your current terminal "/dev/tty"
        ports = glob.glob('/dev/tty[A-Za-z]*')
    elif sys.platform.startswith('darwin'):
        ports = glob.glob('/dev/tty.*')
    else:
        raise EnvironmentError('Unsupported platform')

    result = []
    for port in ports:
        try:
            s = serial.Serial(port)
            s.close()
            result.append(port)
        except (OSError, serial.SerialException):
            pass
    return result
# END DEF

def GetCnx(results):
    # Test USB to MC2T Setting
    global zcom, current_time #, str_val
    rx = False
    print("SERIAL PORTS FOUND:", results)
    print("TESTING PORTS FOR MICROTOUCH-2T ...", end='')
    # print(str_val[0])
    for i in range(256):
        try:
            if(results[i] == "COM1" or results[i] == "COM2"):
                pass
            else:
                portstr = str(results[i])
                ser = serial.Serial(
                    portstr, baudrate=9600, bytesize=8, parity='N', stopbits=1, timeout=None)
                # this causes MC2T to Beep if connection is established.
                wxx = "F5 \r\n"
                wx3 = wxx.encode()
                q = ser.write(wx3)
                time.sleep(0.5)

                qq = ser.read(ser.in_waiting).decode()
                if len(qq) > 5:
                    print(" FOUND MC2T PORT AS: ", portstr)
                    rx = True
                    zcom = portstr
                    ser.close()
                    return(rx)

                # q=ser.write(wx3)
                ser.close()
        except:
            pass
    return(rx)
# END DEF

def SendCmd(Tag2):
    # Singular serial pump communication
    global zcom, LogDump, current_time
    # open serial port xxx
    with serial.Serial(zcom, baudrate=9600, bytesize=8, parity='N', stopbits=1, timeout=None) as ser:
        try:
            ser.flush()
            Tag = Tag2+"\r\n"
            tagx = Tag.encode()
            q = ser.write(tagx)
            time.sleep(0.002) # formerly 0.2
            ChkOK = ser.read(ser.in_waiting).decode().split('\n')
            ser.close()
        except:
            print("!# ERRROR: PUMP COMMUNICATION ERROR #!")
            ser.close()
            ChkOK = ["!#ERROR#!"]
            WriteLog(0)
            LogDump[1] = "SendCmd"
            LogDump[2] = "ZCOM:"+str(zcom)
            LogDump[3] = "tagx:"+tagx
            WriteLog()
    return ChkOK
# END DEF

# AUTOMATIC MOTOR CONTROLLER FUNCTION
def control_motor(magnitude, rates, rateChoice):
    global chkOK, chkOK1, ChkOK, RunPSet, tag_sG, tag_qG, rate001, rate002, rate004, rate008, rate016, rate032
    tag_sG = "*G"
    tag_qG = "?G" 

    for i in range(1, len(magnitude[0])):  # Assuming profileMags and profileRates are 2D lists
        RunPSet[2] = str(magnitude[rateChoice][i])   # RunPSet[2]: AMOUNT(nL)                  
        RunPSet[3] = str(rates[rateChoice][i])  # RunPSet[3]: SPEED   
        ChkOK = SendCmd("V"+RunPSet[2])  # SET TARGET AMOUNT
        ChkOK = SendCmd("R"+RunPSet[3])  # SET DELIVERY RATE
        chkOK1 = SendCmd(tag_sG)  # SEND COMMANDS TO PUMP
        chkOK = SendCmd(tag_qG)   # SEND COMMANDS TO PUMP
# END DEF

# SAMPLING DATA FUNCTIONS #
def samplingFunction():
    global zcom, ChkOK, chkOK, chkOK1, sampling, running, current_time, row, rows1, rows11, time_to_pause, csv_file1, csv_file2, csv_file3, csv_file4, csv_file5, writer
    sample_period = 1/10 # in seconds
    rows1 = []
    rows11 = []
    total_time_elapsed = 0

    # Initialize the serial connection
    ser = serial.Serial(
    port=zcom, 
    baudrate=9600, 
    bytesize=serial.EIGHTBITS, 
    parity=serial.PARITY_NONE, 
    stopbits=serial.STOPBITS_ONE, 
    timeout=None
    )

    buffer = ""  # Initialize the buffer
    time_motor_stop_archive = []  # Initialize the list to store timestamps when the motor stops

    start_time = time.time()
    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION DATA COLLECTION UNDERWAY             ###")
    print("#####################################################################") 
    print("\n")
    while total_time_elapsed <= (time_to_pause):
        buffer += ser.read(ser.in_waiting).decode()

        # If there is at least one full line in the buffer
        while "\n" in buffer:
            line, buffer = buffer.split("\n", 1)  # Split off the first line from the buffer

            if "Motor State: Stopped" in line:
                time_motor_stop = time.time()
                time_motor_stop_archive.append(time_motor_stop)  # Store the timestamp
                print(f"Motor stopped at time {time_motor_stop}")

        current_time = time.time()
        elapsed_time = current_time - start_time
        start_time = current_time  # reset start_time
        total_time_elapsed += elapsed_time  # update total_time_elapsed
        rawBits = nau7802.read()
        kPa = convert_to_kpa(rawBits,nau7802.gain)
        row = [start_time,elapsed_time,total_time_elapsed,rawBits,kPa]
        rows1.append(row)
        rows11.append(row)
    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION DATA COLLECTION COMPLETE             ###")
    print("###           Total time elapsed: {:.2f} seconds                  ###".format(total_time_elapsed))  # print the total time elapsed after the loop
    print("#####################################################################") 
    print("\n")
    return time_motor_stop_archive
# END DEF

def environmentalConditions(samplingRate):
    sample_period = 1.0 / samplingRate  # in seconds
    base_filename = input("Type in desired filename (without CSV extension):")
    csv_filename = base_filename + ".csv"
    total_duration = float(input("Specify the duration for data collection (in seconds):"))
    total_samples = int(total_duration * samplingRate)

    data_samples = []

    zero_channel()

    # Serial connection initialization
    ser = serial.Serial(
        port=zcom,
        baudrate=9600,
        bytesize=serial.EIGHTBITS,
        parity=serial.PARITY_NONE,
        stopbits=serial.STOPBITS_ONE,
        timeout=None
    )

    start_time = time.time()

    expected_time_next_sample = 0

    for i in range(total_samples):
        current_time = time.time()
        elapsed_time = current_time - start_time
        rawBits = nau7802.read()
        kPa = convert_to_kpa(rawBits, nau7802.gain)

        # Store data samples in a list for later writing to the CSV file
        utc_timestamp = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(current_time))
        data_samples.append([utc_timestamp, elapsed_time, rawBits, kPa])

        # Update terminal with timer and pressure data every second
        if i % samplingRate == 0:
            ClrScrn()
            print(f"Time Remaining: {total_duration - elapsed_time:.2f} seconds")
            print(f"Pressure: {kPa:.2f} kPa")

        expected_time_next_sample += sample_period
        time_to_sleep = expected_time_next_sample - (time.time() - start_time)

        if time_to_sleep > 0:  # Avoid negative sleep
            time.sleep(time_to_sleep)

    # Create the directories if they don't exist for CSV
    current_date = datetime.datetime.now().strftime('%Y-%m-%d')
    directory_path_csv = os.path.join("TEST_DATA", "TROUBLESHOOTING", "CSV_FILES", current_date)
    if not os.path.exists(directory_path_csv):
        os.makedirs(directory_path_csv)
    
    csv_filename_path = os.path.join(directory_path_csv, base_filename + ".csv")
    
    # Write data to the CSV file
    with open(csv_filename_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["UTC Timestamp", "Elapsed Time", "Raw Bits", "kPa"])
        writer.writerows(data_samples)

    # Plot the data
    elapsed_times = [sample[1] for sample in data_samples]
    pressures = [sample[3] for sample in data_samples]

    plt.close()
    plt.figure(figsize=(10,6))
    plt.plot(elapsed_times, pressures, label='Pressure (kPa)', color='blue')
    plt.xlabel('Elapsed Time (s)')
    plt.ylabel('Pressure (kPa)')
    plt.title('Pressure Data vs. Elapsed Time')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # Create the directories if they don't exist for IMAGES
    directory_path_img = os.path.join("TEST_DATA", "TROUBLESHOOTING", "IMAGES", current_date)
    if not os.path.exists(directory_path_img):
        os.makedirs(directory_path_img)

    img_filename_path = os.path.join(directory_path_img, base_filename + '.png')

    # Save the plot
    plt.savefig(img_filename_path)
    plt.close(img_filename_path)

    print(f"Data collection completed. Data saved to {csv_filename_path}. Plot saved as {img_filename_path}.")
# END DEF

def get_base_filename():
    # Prompt the user for the base filename
    base_filename = input("Enter the name of the CSV file (w/o CSV extension): ")
    return base_filename
# END DEF

def csv_constRadius(df2, base_name, UTC_runStart, profileRates, profileMags, time_motor_stop_archive, piecewiseTime, rates032):
    # Creating numeric arrays for processing
    profileRates1 = np.array(profileRates).astype(float)
    profileMags1 = np.array(profileMags).astype(float)

    # Define the specific rates and magnitudes
    specific_rates = profileRates1[rates032]
    specific_mags = profileMags1[rates032]

    # Find the row index of the closest UTC_EPOCH to the knownEpoch
    # start_index = (df2['UTC_EPOCH'] - UTC_runStart).abs().idxmin()

    # Create a new column called VOLUME_CONTROL and initialize RATE and VOLUME columns
    df2['VOLUME_CONTROL'] = 0
    df2['RATE'] = 0.0
    df2['VOLUME'] = 0.0

    # Initialize counter for 'specific_mags' and 'specific_rates'
    counter = 0

    # Iterate through each stop time
    for stop_time in time_motor_stop_archive:
        # Calculate the motor start time
        bufferTime = 0.2
        start_time = stop_time - piecewiseTime - bufferTime

        # Find the indices for the motor running period
        motor_start_index = (df2['UTC_EPOCH'] - start_time).abs().idxmin()
        motor_stop_index = (df2['UTC_EPOCH'] - stop_time).abs().idxmin()

        # Set VOLUME_CONTROL to 1 for the motor running period and assign values from 'specific_rates' and 'specific_mags' to 'RATE' and 'VOLUME'
        if counter < len(specific_rates) and counter < len(specific_mags):
            df2.loc[motor_start_index:motor_stop_index-1, 'VOLUME_CONTROL'] = 1
            df2.loc[motor_start_index:motor_stop_index-1, 'RATE'] = specific_rates[counter]
            df2.loc[motor_start_index:motor_stop_index-1, 'VOLUME'] = specific_mags[counter]

        # Increment the counter
        counter += 1

    # Calculate the instantaneous volume transfer at each time step
    df2['VOLUME_INSTANT_TRANSFER'] = df2['RATE'] * df2['TIME_BETWEEN_SAMPLES']

    # Calculate the cumulative volume transfer over time
    df2['VOLUME_TOTAL'] = df2['VOLUME_INSTANT_TRANSFER'].cumsum()

    # Create the directories if they don't exist
    current_date = datetime.datetime.now().strftime('%Y-%m-%d')
    directory_path = os.path.join("TEST_DATA", "CSV_FILES", current_date)
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
    
    # Define the output filename based on the input filename
    output_filename = os.path.join(directory_path, base_name + ".csv")

    # Save DataFrame to the new CSV file
    df2.to_csv(output_filename, index=False)
    
    return output_filename
# END DEF

def csv_constVolume(df2, base_name, UTC_runStart, profileRates, profileMags, time_motor_stop_archive, piecewiseTime):
    # Define the specific rates and magnitudes
    specific_rates = np.array(profileRates).astype(float)
    specific_mags = np.array(profileMags).astype(float)

    # Find the row index of the closest UTC_EPOCH to the knownEpoch
    # start_index = (df2['UTC_EPOCH'] - UTC_runStart).abs().idxmin()

    # Create a new column called VOLUME_CONTROL and initialize RATE and VOLUME columns
    df2['VOLUME_CONTROL'] = 0
    df2['RATE'] = 0.0
    df2['VOLUME'] = 0.0

    # Initialize counter for 'specific_mags' and 'specific_rates'
    counter = 0

    # Iterate through each stop time
    for stop_time in time_motor_stop_archive:
        # Calculate the motor start time
        bufferTime = 0.2
        start_time = stop_time - piecewiseTime - bufferTime

        # Find the indices for the motor running period
        motor_start_index = (df2['UTC_EPOCH'] - start_time).abs().idxmin()
        motor_stop_index = (df2['UTC_EPOCH'] - stop_time).abs().idxmin()

        # Set VOLUME_CONTROL to 1 for the motor running period and assign values from 'specific_rates' and 'specific_mags' to 'RATE' and 'VOLUME'
        if counter < len(specific_rates) and counter < len(specific_mags):
            df2.loc[motor_start_index:motor_stop_index-1, 'VOLUME_CONTROL'] = 1
            df2.loc[motor_start_index:motor_stop_index-1, 'RATE'] = specific_rates[counter]
            df2.loc[motor_start_index:motor_stop_index-1, 'VOLUME'] = specific_mags[counter]

        # Increment the counter
        counter += 1

    # Calculate the instantaneous volume transfer at each time step
    df2['VOLUME_INSTANT_TRANSFER'] = df2['RATE'] * df2['TIME_BETWEEN_SAMPLES']

    # Calculate the cumulative volume transfer over time
    df2['VOLUME_TOTAL'] = df2['VOLUME_INSTANT_TRANSFER'].cumsum()

    # Create the directories if they don't exist
    current_date = datetime.datetime.now().strftime('%Y-%m-%d')
    directory_path = os.path.join("TEST_DATA", "CSV_FILES", current_date)
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
    
    # Define the output filename based on the input filename
    output_filename = os.path.join(directory_path, base_name + ".csv")

    # Save DataFrame to the new CSV file
    df2.to_csv(output_filename, index=False)
    
    return output_filename
# END DEF

def convert_to_kpa(data_y, gain):
    mV_psi_V = 0.2584 # 0.2584 mV/psi/V
    V_ref = 3.0 # reference voltage
    fullScale_factor = 2
    mV_psi = mV_psi_V * V_ref # mV/psi
    V_psi = mV_psi / 1000 # V/psi
    V_kPa = V_psi / 6.89476 # V/kPa
    nBits = 2**24 # 24-bit ADC, i.e. 2^24 bits
    measuredVoltage = (data_y / nBits) * V_ref
    actualVoltage = (measuredVoltage  / gain) / fullScale_factor
    data_y_converted = (actualVoltage / V_kPa)
    return data_y_converted
# END DEF

def plot_and_save1(filename, base):
    # Read the CSV file into a DataFrame
    df3 = pd.read_csv(filename)

    # Close previous plots and create subplot architecture (2x2)
    plt.clf()
    fig, axs = plt.subplots(2, 2, figsize=(15, 10))

    # Set the titles
    titles = ['Pressure vs Time', 'Rate vs Time', 'Volume vs Time', 'Pressure vs Volume']

    # 1. Pressure vs Time
    axs[0, 0].plot(df3['TOTAL_TEST_TIME'], df3['KPA'])
    axs[0, 0].set(xlabel='Time (s)', ylabel='Pressure (kPa)')

    # Filter rows where VOLUME_CONTROL is 1
    df_filtered = df3[df3['VOLUME_CONTROL'] == 1]

    # 4. Pressure vs Volume
    axs[1, 1].plot(df_filtered['VOLUME_TOTAL'], df_filtered['KPA'])
    axs[1, 1].set(xlabel='Volume', ylabel='Pressure (kPa)')

    # 3. Volume vs Time
    axs[1, 0].plot(df3['TOTAL_TEST_TIME'], df3['VOLUME_TOTAL'])
    axs[1, 0].set(xlabel='Time (s)', ylabel='Volume')

    # 2. Rate vs Time
    axs[0, 1].plot(df3['TOTAL_TEST_TIME'], df3['RATE'])
    axs[0, 1].set(xlabel='Time (s)', ylabel='Rate')

    # Set titles for each subplot
    for i, ax in enumerate(axs.flat):
        ax.set_title(titles[i])

    # Adjust layout for non-overlapping
    plt.tight_layout()

    # Generate the directory and filename
    date_str = time.strftime('%Y-%m-%d', time.gmtime(time.time()))
    directory = os.path.join('TEST_DATA', 'IMAGES', date_str)
    base = os.path.basename(base)  # Extract the filename without path
    base = os.path.splitext(base)[0]
    png_file_name = os.path.join(directory, base + '_plots.png')
    
    # Ensure the directory exists
    os.makedirs(directory, exist_ok=True)

    # Save the plot as a PNG file
    plt.savefig(png_file_name)
    plt.close(fig)

    return png_file_name  # return the filename for further use if required
# END DEF

def injectFluid_constVolume_customTime(): 
    global rate_cVolume, timeRunning, a_stop, vol_stop, RunPSet, LogDump, running, sampling, ChkOk, ChkOK, chkOK, tag_sG, tag_qG, nIterations_time, UTC_runStart
    # Defaults RunPSet: SETTINGS FOR AUTOMATED PROGRAM
    a_stop = input("Input a final radius (mm) for the cavity: ")
    vol_stop = ((4/3) * np.pi * (float(a_stop))**3) * 1000 # convert mm^3 to nL
    timeRunning = input('Input required time (sec) to perform the infusion: ')
    rate_cVolume = vol_stop / float(timeRunning)

    RunPSet[0] = "1"  # ACTIVE PUMP                            #
    RunPSet[1] = "INFUSE"  # RunPSet[1]:INFUSE or WITHDRAW       #
    RunPSet[2] = str(vol_stop)  # RunPSet[2]: AMOUNT(nL) -- need waterCal_short
    RunPSet[3] = str(rate_cVolume)  # RunPSet[3]: SPEED -- need waterCal_short   
    RunPSet[4] = "SEC"  # RunPSet[4]: UNITS :nL/SEC  or nL/MIN   #
    RunPSet[5] = nIterations  # RunPSet[5]: PAUSE (x/100) SECS         #
    RunPSet[6] = "0"  # RunPSet[6]: # TIMES TO REPEAT, 0=NONE)   #
    
    xflag = 0
    RPS6 = int(RunPSet[6])
    tag1 = "I"  # not really sure what this does. copied from og code, need to investigate
    tag2 = "EI" 
    tag3 = "S" 
    # SET UP TO RUN PROGRAM
    ChkOK = SendCmd("L"+RunPSet[0])  # SET ACTIVE PUMP
    if(RunPSet[1] == "INFUSE"):
        tag1 = "I"  # SET DIRECTION INFUSE
        tag2 = "EI"  # SET VOLUME COUNTER DELIVERED
    if(RunPSet[1] == "WITHDRAW"):
        tag1 = "W"  # SET DIRECTION WITHDRAW
        tag2 = "EN"  # SET VOLUME COUNTER REMAINING

    ChkOK = SendCmd(tag1)  # SET DIRECTION
    ChkOK = SendCmd(tag2)  # SET VOLUME COUNTER

    if(RunPSet[4] == "SEC"):
        tag3 = "S"
    if(RunPSet[4] == "MIN"):
        tag3 = "M"
    ChkOK = SendCmd(tag3)  # SET UNITS

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD UNDERWAY             ###")
    print("#####################################################################") 
    print("\n")

    tag_sG = "*G"
    tag_qG = "?G"

    syringe_and_stepping()

    # Record the start time
    time_Cmd0 = time.time()
    # Send the command
    ChkOk = SendCmd("A"+RunPSet[5])  # SEND COMMAND TO PUMP
    # Record the end time
    time_Cmd1 = time.time()
    UTC_runStart = time_Cmd1+nIterations_tSequence

    ChkOK = SendCmd("V"+RunPSet[2])  # SET TARGET AMOUNT
    ChkOK = SendCmd("R"+RunPSet[3])  # SET DELIVERY RATE
    chkOK1 = SendCmd(tag_sG) # SEND COMMANDS TO PUMP
    chkOK = SendCmd(tag_qG) # SEND COMMANDS TO PUMP

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD COMPLETE             ###")
    print("#####################################################################") 
    print("\n")    
    sampling = False
# END DEF

def injectFluid_constVolume_customRate(): 
    global rate_cVolume, timeRunning, a_stop, vol_stop, RunPSet, LogDump, running, sampling, ChkOk, ChkOK, chkOK, tag_sG, tag_qG, nIterations_time, UTC_runStart
    # Defaults RunPSet: SETTINGS FOR AUTOMATED PROGRAM
    a_stop = input("Input a final radius (mm) for the cavity: ")
    vol_stop = ((4/3) * np.pi * (float(a_stop))**3) * 1000 # mm^3 to nL
    rate_cVolume = input('Input required rate (nL/sec) to perform the infusion: ')
    timeRunning = vol_stop / float(rate_cVolume)

    RunPSet[0] = "1"  # ACTIVE PUMP                            #
    RunPSet[1] = "INFUSE"  # RunPSet[1]:INFUSE or WITHDRAW       #
    RunPSet[2] = str(vol_stop)  # RunPSet[2]: AMOUNT(nL) -- need waterCal_short
    RunPSet[3] = str(rate_cVolume)  # RunPSet[3]: SPEED -- need waterCal_short   
    RunPSet[4] = "SEC"  # RunPSet[4]: UNITS :nL/SEC  or nL/MIN   #
    RunPSet[5] = nIterations  # RunPSet[5]: PAUSE (x/100) SECS         #
    RunPSet[6] = "0"  # RunPSet[6]: # TIMES TO REPEAT, 0=NONE)   #
    
    xflag = 0
    RPS6 = int(RunPSet[6])
    tag1 = "I"  # not really sure what this does. copied from og code, need to investigate
    tag2 = "EI" 
    tag3 = "S" 
    # SET UP TO RUN PROGRAM
    ChkOK = SendCmd("L"+RunPSet[0])  # SET ACTIVE PUMP
    if(RunPSet[1] == "INFUSE"):
        tag1 = "I"  # SET DIRECTION INFUSE
        tag2 = "EI"  # SET VOLUME COUNTER DELIVERED
    if(RunPSet[1] == "WITHDRAW"):
        tag1 = "W"  # SET DIRECTION WITHDRAW
        tag2 = "EN"  # SET VOLUME COUNTER REMAINING

    ChkOK = SendCmd(tag1)  # SET DIRECTION
    ChkOK = SendCmd(tag2)  # SET VOLUME COUNTER

    if(RunPSet[4] == "SEC"):
        tag3 = "S"
    if(RunPSet[4] == "MIN"):
        tag3 = "M"
    ChkOK = SendCmd(tag3)  # SET UNITS

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD UNDERWAY             ###")
    print("#####################################################################") 
    print("\n")

    tag_sG = "*G"
    tag_qG = "?G"

    syringe_and_stepping()

    # Record the start time
    time_Cmd0 = time.time()
    # Send the command
    ChkOk = SendCmd("A"+RunPSet[5])  # SEND COMMAND TO PUMP
    # Record the end time
    time_Cmd1 = time.time()
    UTC_runStart = time_Cmd1+nIterations_tSequence

    ChkOK = SendCmd("V"+RunPSet[2])  # SET TARGET AMOUNT
    ChkOK = SendCmd("R"+RunPSet[3])  # SET DELIVERY RATE
    chkOK1 = SendCmd(tag_sG) # SEND COMMANDS TO PUMP
    chkOK = SendCmd(tag_qG) # SEND COMMANDS TO PUMP

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD COMPLETE             ###")
    print("#####################################################################") 
    print("\n")    
    sampling = False
# END DEF

# MOTOR PROFILE FUNCTIONS #
def injectFluid_constRadius_001mmps(): 
    global chkOK, RunPSet, LogDump, running, sampling, ChkOK, ChkOk, tag_sG, tag_qG, nIterations_time, UTC_runStart
    # Defaults RunPSet: SETTINGS FOR AUTOMATED PROGRAM
    RunPSet[0] = "1"  # ACTIVE PUMP                            #
    RunPSet[1] = "INFUSE"  # RunPSet[1]:INFUSE or WITHDRAW       #
    RunPSet[2] = str(profileMags[rate001][0])  # RunPSet[2]: AMOUNT(nL) 
    RunPSet[3] = str(profileRates[rate001][0])  # RunPSet[3]: SPEED 
    RunPSet[4] = "SEC"  # RunPSet[4]: UNITS :nL/SEC  or nL/MIN   #
    RunPSet[5] = nIterations  # RunPSet[5]: PAUSE (x/100) SECS         #
    RunPSet[6] = "0"  # RunPSet[6]: # TIMES TO REPEAT, 0=NONE)   #

    xflag = 0
    RPS6 = int(RunPSet[6])
    tag1 = "I"  # not really sure what this does. copied from og code, need to investigate
    tag2 = "EI" 
    tag3 = "S" 
    # SET UP TO RUN PROGRAM
    ChkOK = SendCmd("L"+RunPSet[0])  # SET ACTIVE PUMP
    if(RunPSet[1] == "INFUSE"):
        tag1 = "I"  # SET DIRECTION INFUSE
        tag2 = "EI"  # SET VOLUME COUNTER DELIVERED
    if(RunPSet[1] == "WITHDRAW"):
        tag1 = "W"  # SET DIRECTION WITHDRAW
        tag2 = "EN"  # SET VOLUME COUNTER REMAINING

    ChkOK = SendCmd(tag1)  # SET DIRECTION
    ChkOK = SendCmd(tag2)  # SET VOLUME COUNTER

    if(RunPSet[4] == "SEC"):
        tag3 = "S"
    if(RunPSet[4] == "MIN"):
        tag3 = "M"
    ChkOK = SendCmd(tag3)  # SET UNITS

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD UNDERWAY             ###")
    print("#####################################################################") 
    print("\n")

    tag_sG = "*G"
    tag_qG = "?G"

    syringe_and_stepping()

    # Record the start time
    time_Cmd0 = time.time()
    # Send the command
    ChkOk = SendCmd("A"+RunPSet[5])  # SEND COMMAND TO PUMP
    # Record the end time
    time_Cmd1 = time.time()
    UTC_runStart = time_Cmd1+nIterations_tSequence

    ChkOK = SendCmd("V"+RunPSet[2])  # SET TARGET AMOUNT
    ChkOK = SendCmd("R"+RunPSet[3])  # SET DELIVERY RATE
    chkOK1 = SendCmd(tag_sG) # SEND COMMANDS TO PUMP
    chkOK = SendCmd(tag_qG) # SEND COMMANDS TO PUMP

    control_motor(profileMags,profileRates,rate001)

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD COMPLETE             ###")
    print("#####################################################################") 
    print("\n")    
    sampling = False
# END DEF

def injectFluid_constRadius_002mmps(): 
    global chgConfig, RunPSet, LogDump, running, sampling, ChkOk, ChkOK, chkOK, tag_sG, tag_qG, nIterations_time, UTC_runStart
    # Defaults RunPSet: SETTINGS FOR AUTOMATED PROGRAM
    RunPSet[0] = "1"  # ACTIVE PUMP                            #
    RunPSet[1] = "INFUSE"  # RunPSet[1]:INFUSE or WITHDRAW       #
    RunPSet[2] = str(profileMags[rate002][0])  # RunPSet[2]: AMOUNT(nL) -- need waterCal_short
    RunPSet[3] = str(profileRates[rate002][0])  # RunPSet[3]: SPEED -- need waterCal_short   
    RunPSet[4] = "SEC"  # RunPSet[4]: UNITS :nL/SEC  or nL/MIN   #
    RunPSet[5] = nIterations  # RunPSet[5]: PAUSE (x/100) SECS         #
    RunPSet[6] = "0"  # RunPSet[6]: # TIMES TO REPEAT, 0=NONE)   #
    
    xflag = 0
    RPS6 = int(RunPSet[6])
    tag1 = "I"  # not really sure what this does. copied from og code, need to investigate
    tag2 = "EI" 
    tag3 = "S" 
    # SET UP TO RUN PROGRAM
    ChkOK = SendCmd("L"+RunPSet[0])  # SET ACTIVE PUMP
    if(RunPSet[1] == "INFUSE"):
        tag1 = "I"  # SET DIRECTION INFUSE
        tag2 = "EI"  # SET VOLUME COUNTER DELIVERED
    if(RunPSet[1] == "WITHDRAW"):
        tag1 = "W"  # SET DIRECTION WITHDRAW
        tag2 = "EN"  # SET VOLUME COUNTER REMAINING

    ChkOK = SendCmd(tag1)  # SET DIRECTION
    ChkOK = SendCmd(tag2)  # SET VOLUME COUNTER

    if(RunPSet[4] == "SEC"):
        tag3 = "S"
    if(RunPSet[4] == "MIN"):
        tag3 = "M"
    ChkOK = SendCmd(tag3)  # SET UNITS

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD UNDERWAY             ###")
    print("#####################################################################") 
    print("\n")

    tag_sG = "*G"
    tag_qG = "?G"

    chgConfig = syringe_and_stepping()

    # Record the start time
    time_Cmd0 = time.time()
    # Send the command
    ChkOk = SendCmd("A"+RunPSet[5])  # SEND COMMAND TO PUMP
    # Record the end time
    time_Cmd1 = time.time()
    UTC_runStart = time_Cmd1+nIterations_tSequence

    ChkOK = SendCmd("V"+RunPSet[2])  # SET TARGET AMOUNT
    ChkOK = SendCmd("R"+RunPSet[3])  # SET DELIVERY RATE
    chkOK1 = SendCmd(tag_sG) # SEND COMMANDS TO PUMP
    chkOK = SendCmd(tag_qG) # SEND COMMANDS TO PUMP

    control_motor(profileMags,profileRates,rate002)

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD COMPLETE             ###")
    print("#####################################################################") 
    print("\n")    
    sampling = False
# END DEF

def injectFluid_constRadius_004mmps(): 
    global RunPSet, LogDump, running, sampling, ChkOk, ChkOK, chkOK, tag_sG, tag_qG, nIterations_time, UTC_runStart
    # Defaults RunPSet: SETTINGS FOR AUTOMATED PROGRAM
    RunPSet[0] = "1"  # ACTIVE PUMP                            #
    RunPSet[1] = "INFUSE"  # RunPSet[1]:INFUSE or WITHDRAW       #
    RunPSet[2] = str(profileMags[rate004][0])  # RunPSet[2]: AMOUNT(nL) -- need waterCal_short
    RunPSet[3] = str(profileRates[rate004][0])  # RunPSet[3]: SPEED -- need waterCal_short   
    RunPSet[4] = "SEC"  # RunPSet[4]: UNITS :nL/SEC  or nL/MIN   #
    RunPSet[5] = nIterations  # RunPSet[5]: PAUSE (x/100) SECS         #
    RunPSet[6] = "0"  # RunPSet[6]: # TIMES TO REPEAT, 0=NONE)   #
    
    xflag = 0
    RPS6 = int(RunPSet[6])
    tag1 = "I" 
    tag2 = "EI" 
    tag3 = "S" 
    # SET UP TO RUN PROGRAM
    ChkOK = SendCmd("L"+RunPSet[0])  # SET ACTIVE PUMP
    if(RunPSet[1] == "INFUSE"):
        tag1 = "I"  # SET DIRECTION INFUSE
        tag2 = "EI"  # SET VOLUME COUNTER DELIVERED
    if(RunPSet[1] == "WITHDRAW"):
        tag1 = "W"  # SET DIRECTION WITHDRAW
        tag2 = "EN"  # SET VOLUME COUNTER REMAINING

    ChkOK = SendCmd(tag1)  # SET DIRECTION
    ChkOK = SendCmd(tag2)  # SET VOLUME COUNTER

    if(RunPSet[4] == "SEC"):
        tag3 = "S"
    if(RunPSet[4] == "MIN"):
        tag3 = "M"
    ChkOK = SendCmd(tag3)  # SET UNITS

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD UNDERWAY             ###")
    print("#####################################################################") 
    print("\n")

    tag_sG = "*G"
    tag_qG = "?G"

    syringe_and_stepping()

    # Record the start time
    time_Cmd0 = time.time()
    # Send the command
    ChkOk = SendCmd("A"+RunPSet[5])  # SEND COMMAND TO PUMP
    # Record the end time
    time_Cmd1 = time.time()
    UTC_runStart = time_Cmd1+nIterations_tSequence

    ChkOK = SendCmd("V"+RunPSet[2])  # SET TARGET AMOUNT
    ChkOK = SendCmd("R"+RunPSet[3])  # SET DELIVERY RATE
    chkOK1 = SendCmd(tag_sG) # SEND COMMANDS TO PUMP
    chkOK = SendCmd(tag_qG) # SEND COMMANDS TO PUMP

    control_motor(profileMags,profileRates,rate004)

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD COMPLETE             ###")
    print("#####################################################################") 
    print("\n")    
    sampling = False
# END DEF

def injectFluid_constRadius_008mmps(): 
    global RunPSet, LogDump, running, sampling, ChkOk, ChkOK, chkOK, tag_sG, tag_qG, nIterations_time, UTC_runStart
    # Defaults RunPSet: SETTINGS FOR AUTOMATED PROGRAM
    RunPSet[0] = "1"  # ACTIVE PUMP                            #
    RunPSet[1] = "INFUSE"  # RunPSet[1]:INFUSE or WITHDRAW       #
    RunPSet[2] = str(profileMags[rate008][0])  # RunPSet[2]: AMOUNT(nL) -- need waterCal_short
    RunPSet[3] = str(profileRates[rate008][0])  # RunPSet[3]: SPEED -- need waterCal_short   
    RunPSet[4] = "SEC"  # RunPSet[4]: UNITS :nL/SEC  or nL/MIN   #
    RunPSet[5] = nIterations  # RunPSet[5]: PAUSE (x/100) SECS         #
    RunPSet[6] = "0"  # RunPSet[6]: # TIMES TO REPEAT, 0=NONE)   #
    
    xflag = 0
    RPS6 = int(RunPSet[6])
    tag1 = "I"  # not really sure what this does. copied from og code, need to investigate
    tag2 = "EI" 
    tag3 = "S" 
    # SET UP TO RUN PROGRAM
    ChkOK = SendCmd("L"+RunPSet[0])  # SET ACTIVE PUMP
    if(RunPSet[1] == "INFUSE"):
        tag1 = "I"  # SET DIRECTION INFUSE
        tag2 = "EI"  # SET VOLUME COUNTER DELIVERED
    if(RunPSet[1] == "WITHDRAW"):
        tag1 = "W"  # SET DIRECTION WITHDRAW
        tag2 = "EN"  # SET VOLUME COUNTER REMAINING

    ChkOK = SendCmd(tag1)  # SET DIRECTION
    ChkOK = SendCmd(tag2)  # SET VOLUME COUNTER

    if(RunPSet[4] == "SEC"):
        tag3 = "S"
    if(RunPSet[4] == "MIN"):
        tag3 = "M"
    ChkOK = SendCmd(tag3)  # SET UNITS

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD UNDERWAY             ###")
    print("#####################################################################") 
    print("\n")

    tag_sG = "*G"
    tag_qG = "?G"

    syringe_and_stepping()

    # Record the start time
    time_Cmd0 = time.time()
    # Send the command
    ChkOk = SendCmd("A"+RunPSet[5])  # SEND COMMAND TO PUMP
    # Record the end time
    time_Cmd1 = time.time()
    UTC_runStart = time_Cmd1+nIterations_tSequence

    ChkOK = SendCmd("V"+RunPSet[2])  # SET TARGET AMOUNT
    ChkOK = SendCmd("R"+RunPSet[3])  # SET DELIVERY RATE
    chkOK1 = SendCmd(tag_sG) # SEND COMMANDS TO PUMP
    chkOK = SendCmd(tag_qG) # SEND COMMANDS TO PUMP

    control_motor(profileMags,profileRates,rate008)

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD COMPLETE             ###")
    print("#####################################################################") 
    print("\n")    
    sampling = False
# END DEF

def injectFluid_constRadius_016mmps(): 
    global RunPSet, LogDump, running, sampling, ChkOk, ChkOK, chkOK, tag_sG, tag_qG, nIterations_time, UTC_runStart
    # Defaults RunPSet: SETTINGS FOR AUTOMATED PROGRAM
    RunPSet[0] = "1"  # ACTIVE PUMP                            #
    RunPSet[1] = "INFUSE"  # RunPSet[1]:INFUSE or WITHDRAW       #
    RunPSet[2] = str(profileMags[rate016][0])  # RunPSet[2]: AMOUNT(nL) -- need waterCal_short
    RunPSet[3] = str(profileRates[rate016][0])  # RunPSet[3]: SPEED -- need waterCal_short   
    RunPSet[4] = "SEC"  # RunPSet[4]: UNITS :nL/SEC  or nL/MIN   #
    RunPSet[5] = nIterations  # RunPSet[5]: PAUSE (x/100) SECS         #
    RunPSet[6] = "0"  # RunPSet[6]: # TIMES TO REPEAT, 0=NONE)   #
    
    xflag = 0
    RPS6 = int(RunPSet[6])
    tag1 = "I"  # not really sure what this does. copied from og code, need to investigate
    tag2 = "EI" 
    tag3 = "S" 
    # SET UP TO RUN PROGRAM
    ChkOK = SendCmd("L"+RunPSet[0])  # SET ACTIVE PUMP
    if(RunPSet[1] == "INFUSE"):
        tag1 = "I"  # SET DIRECTION INFUSE
        tag2 = "EI"  # SET VOLUME COUNTER DELIVERED
    if(RunPSet[1] == "WITHDRAW"):
        tag1 = "W"  # SET DIRECTION WITHDRAW
        tag2 = "EN"  # SET VOLUME COUNTER REMAINING

    ChkOK = SendCmd(tag1)  # SET DIRECTION
    ChkOK = SendCmd(tag2)  # SET VOLUME COUNTER

    if(RunPSet[4] == "SEC"):
        tag3 = "S"
    if(RunPSet[4] == "MIN"):
        tag3 = "M"
    ChkOK = SendCmd(tag3)  # SET UNITS

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD UNDERWAY             ###")
    print("#####################################################################") 
    print("\n")

    tag_sG = "*G"
    tag_qG = "?G"

    syringe_and_stepping()

    # Record the start time
    time_Cmd0 = time.time()
    # Send the command
    ChkOk = SendCmd("A"+RunPSet[5])  # SEND COMMAND TO PUMP
    # Record the end time
    time_Cmd1 = time.time()
    UTC_runStart = time_Cmd1+nIterations_tSequence

    ChkOK = SendCmd("V"+RunPSet[2])  # SET TARGET AMOUNT
    ChkOK = SendCmd("R"+RunPSet[3])  # SET DELIVERY RATE
    chkOK1 = SendCmd(tag_sG) # SEND COMMANDS TO PUMP
    chkOK = SendCmd(tag_qG) # SEND COMMANDS TO PUMP

    control_motor(profileMags,profileRates,rate016)

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD COMPLETE             ###")
    print("#####################################################################") 
    print("\n")    
    sampling = False
# END DEF

def injectFluid_constRadius_032mmps(): 
    global RunPSet, LogDump, running, sampling, ChkOk, ChkOK, chkOK, chkOK1, tag_sG, tag_qG, nIterations_time, UTC_runStart
    # Defaults RunPSet: SETTINGS FOR AUTOMATED PROGRAM
    RunPSet[0] = "1"  # ACTIVE PUMP                            #
    RunPSet[1] = "INFUSE"  # RunPSet[1]:INFUSE or WITHDRAW       #
    RunPSet[2] = str(profileMags[rate032][0])  # RunPSet[2]: AMOUNT(nL) -- need waterCal_short
    RunPSet[3] = str(profileRates[rate032][0])  # RunPSet[3]: SPEED -- need waterCal_short   
    RunPSet[4] = "SEC"  # RunPSet[4]: UNITS :nL/SEC  or nL/MIN   #
    RunPSet[5] = nIterations  # RunPSet[5]: PAUSE (x/100) SECS         #
    RunPSet[6] = "0"  # RunPSet[6]: # TIMES TO REPEAT, 0=NONE)   #
    
    xflag = 0
    RPS6 = int(RunPSet[6])
    tag1 = "I"  # not really sure what this does. copied from og code, need to investigate
    tag2 = "EI" 
    tag3 = "S" 
    # SET UP TO RUN PROGRAM
    ChkOK = SendCmd("L"+RunPSet[0])  # SET ACTIVE PUMP
    if(RunPSet[1] == "INFUSE"):
        tag1 = "I"  # SET DIRECTION INFUSE
        tag2 = "EI"  # SET VOLUME COUNTER DELIVERED
    if(RunPSet[1] == "WITHDRAW"):
        tag1 = "W"  # SET DIRECTION WITHDRAW
        tag2 = "EN"  # SET VOLUME COUNTER REMAINING

    ChkOK = SendCmd(tag1)  # SET DIRECTION
    ChkOK = SendCmd(tag2)  # SET VOLUME COUNTER

    if(RunPSet[4] == "SEC"):
        tag3 = "S"
    if(RunPSet[4] == "MIN"):
        tag3 = "M"
    ChkOK = SendCmd(tag3)  # SET UNITS

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD UNDERWAY             ###")
    print("#####################################################################") 
    print("\n")

    tag_sG = "*G"
    tag_qG = "?G"

    syringe_and_stepping()

    # Record the start time
    time_Cmd0 = time.time()
    # Send the command
    ChkOk = SendCmd("A"+RunPSet[5])  # SEND COMMAND TO PUMP
    # Record the end time
    time_Cmd1 = time.time()
    UTC_runStart = time_Cmd1+nIterations_tSequence

    ChkOK = SendCmd("V"+RunPSet[2])  # SET TARGET AMOUNT
    ChkOK = SendCmd("R"+RunPSet[3])  # SET DELIVERY RATE
    chkOK1 = SendCmd(tag_sG) # SEND COMMANDS TO PUMP
    chkOK = SendCmd(tag_qG) # SEND COMMANDS TO PUMP

    control_motor(profileMags,profileRates,rate032)

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD COMPLETE             ###")
    print("#####################################################################") 
    print("\n")    
    sampling = False
# END DEF

'''
def withdrawFluid_001mmps(): 
    global RunPSet, LogDump, running, sampling, ChkOK, ChkOk, chkOK, tag_sG, tag_qG, nIterations_time, UTC_runStart
    # Defaults RunPSet: SETTINGS FOR AUTOMATED PROGRAM
    RunPSet[0] = "1"  # ACTIVE PUMP                            
    RunPSet[1] = "WITHDRAW"  # RunPSet[1]:INFUSE or WITHDRAW   
    RunPSet[2] = str(profileMags[rate001][0])  # RunPSet[2]: AMOUNT(nL) 
    RunPSet[3] = str(profileRates[rate001][0])  # RunPSet[3]: SPEED 
    RunPSet[4] = "SEC"  # RunPSet[4]: UNITS :nL/SEC  or nL/MIN 
    RunPSet[5] = nIterations  # RunPSet[5]: PAUSE (x/100) SECS 
    RunPSet[6] = "0"  # RunPSet[6]: # TIMES TO REPEAT, 0=NONE) 

    xflag = 0
    RPS6 = int(RunPSet[6])
    tag1 = "I"  # not really sure what this does. copied from og code, need to investigate
    tag2 = "EI" 
    tag3 = "S" 
    # SET UP TO RUN PROGRAM
    ChkOK = SendCmd("L"+RunPSet[0])  # SET ACTIVE PUMP
    if(RunPSet[1] == "INFUSE"):
        tag1 = "I"  # SET DIRECTION INFUSE
        tag2 = "EI"  # SET VOLUME COUNTER DELIVERED
    if(RunPSet[1] == "WITHDRAW"):
        tag1 = "W"  # SET DIRECTION WITHDRAW
        tag2 = "EN"  # SET VOLUME COUNTER REMAINING

    ChkOK = SendCmd(tag1)  # SET DIRECTION
    ChkOK = SendCmd(tag2)  # SET VOLUME COUNTER

    if(RunPSet[4] == "SEC"):
        tag3 = "S"
    if(RunPSet[4] == "MIN"):
        tag3 = "M"
    ChkOK = SendCmd(tag3)  # SET UNITS

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD UNDERWAY             ###")
    print("#####################################################################") 
    print("\n")

    tag_sG = "*G"
    tag_qG = "?G"

    syringe_and_stepping()

    # Record the start time
    time_Cmd0 = time.time()
    # Send the command
    ChkOk = SendCmd("A"+RunPSet[5])  # SEND COMMAND TO PUMP
    # Record the end time
    time_Cmd1 = time.time()
    UTC_runStart = time_Cmd1+nIterations_tSequence

    ChkOK = SendCmd("V"+RunPSet[2])  # SET TARGET AMOUNT
    ChkOK = SendCmd("R"+RunPSet[3])  # SET DELIVERY RATE
    chkOK1 = SendCmd(tag_sG) # SEND COMMANDS TO PUMP
    chkOK = SendCmd(tag_qG) # SEND COMMANDS TO PUMP

    control_motor(profileMags,profileRates,rate001)

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION PROTOCOL UPLOAD COMPLETE             ###")
    print("#####################################################################") 
    print("\n")    
    sampling = False
'''
# END DEF

def convert_timestamp(timestamp):
    UTC = time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime(timestamp))
    return UTC

def primeMotor(): 
    global RunPSet, LogDump, running, sampling, ChkOk, ChkOK, chkOK, tag_sG, tag_qG, nIterations_time, UTC_runStart
    # Defaults RunPSet: SETTINGS FOR AUTOMATED PROGRAM
    RunPSet[0] = "1"  # ACTIVE PUMP                            #
    RunPSet[1] = "INFUSE"  # RunPSet[1]:INFUSE or WITHDRAW       #
    #RunPSet[2] = "83.4"  # RunPSet[2]: AMOUNT(nL) -- need waterCal_short (0.083323mm^3 -- 0.1mm of travel)
    RunPSet[2] = "50.1"  # RunPSet[2]: AMOUNT(nL) -- need waterCal_short (0.083323mm^3 -- 0.1mm of travel); need to confirm
    RunPSet[3] = "100"  # RunPSet[3]: SPEED -- need waterCal_short   
    RunPSet[4] = "SEC"  # RunPSet[4]: UNITS :nL/SEC  or nL/MIN   #
    RunPSet[5] = nIterations  # RunPSet[5]: PAUSE (x/100) SECS         #
    RunPSet[6] = "0"  # RunPSet[6]: # TIMES TO REPEAT, 0=NONE)   #
    
    xflag = 0
    RPS6 = int(RunPSet[6])
    tag1 = "I"  # not really sure what this does. copied from og code, need to investigate
    tag2 = "EI" 
    tag3 = "S" 
    # SET UP TO RUN PROGRAM
    ChkOK = SendCmd("L"+RunPSet[0])  # SET ACTIVE PUMP
    if(RunPSet[1] == "INFUSE"):
        tag1 = "I"  # SET DIRECTION INFUSE
        tag2 = "EI"  # SET VOLUME COUNTER DELIVERED
    if(RunPSet[1] == "WITHDRAW"):
        tag1 = "W"  # SET DIRECTION WITHDRAW
        tag2 = "EN"  # SET VOLUME COUNTER REMAINING

    ChkOK = SendCmd(tag1)  # SET DIRECTION
    ChkOK = SendCmd(tag2)  # SET VOLUME COUNTER

    if(RunPSet[4] == "SEC"):
        tag3 = "S"
    if(RunPSet[4] == "MIN"):
        tag3 = "M"
    ChkOK = SendCmd(tag3)  # SET UNITS

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION MOTOR PRIMING UNDERWAY               ###")
    print("#####################################################################") 
    print("\n")

    tag_sG = "*G"
    tag_qG = "?G"

    syringe_and_stepping()

    ChkOK = SendCmd("V"+RunPSet[2])  # SET TARGET AMOUNT
    ChkOK = SendCmd("R"+RunPSet[3])  # SET DELIVERY RATE
    chkOK1 = SendCmd(tag_sG) # SEND COMMANDS TO PUMP
    chkOK = SendCmd(tag_qG) # SEND COMMANDS TO PUMP

    print("\n")   
    print("#####################################################################")
    print("###        DIGITAL PALPATION MOTOR PRIMING COMPLETE               ###")
    print("#####################################################################") 
    print("\n")    
# END DEF

# OPEN AND CHECK PORTS #
results = serial_ports()
rx = GetCnx(results)

def syringe_and_stepping(): 
    global chgConfig
    # CONFIRM SYRINGE TYPE
    tenMicro = "T4"

    # fiftyMicro = "T7"
    twofiftyMicro = "T9" #?
    smoothDrive = "BS"
    #maxDrive = "BT"
    #chgConfig[1] = tenMicro
    chgConfig[1] = twofiftyMicro
    chgConfig[9] = smoothDrive
    #chgConfig[9]= maxDrive
    # Change Syringe Parameters
    ChkOK=SendCmd(chgConfig[9])
    ChkOK = SendCmd(chgConfig[1])
    return chgConfig
# END DEF

def test_log_setup(variables={}):
    log_dir = "TEST_LOGS"
    log_master_dir = os.path.join(log_dir, "TEST_LOGS_MASTER")
    log_date_dir = os.path.join(log_dir, "TEST_LOGS_BY_DATE")

    os.makedirs(log_master_dir, exist_ok=True)
    os.makedirs(log_date_dir, exist_ok=True)

    master_log_file = os.path.join(log_master_dir, "microPump_masterNotes.csv")
    fields = ["DATE_UTC", "DATE_TIME.TIME()", "FILENAME", "TEST_PURPOSE", "SAMPLE NUMBER", "PRESSURE_SENSOR SERIAL NUMBER", "FABRICATION DATE", 
              "KNOWN BASE/CROSS_LINKER RATIO", "FUNCTION USED", "MOTOR OPERATION",
              "NEEDLE LENGTH", "NEEDLE GAUGE", "SYRINGE TYPE", "WAIT TIME POST INSERTION", "RELAXATION TIME POST INFUSION", "TEMPERATURE", 
              "HUMIDITY", "LEAKS", "OTHER NOTES"]

    previous_entries = {field: '' for field in fields}

    if not os.path.exists(master_log_file):
        with open(master_log_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
    else:
        with open(master_log_file, 'r', newline='') as f:
            reader = csv.DictReader(f)
            try:
                last_row = list(reader)[-1]
                for field in fields:
                    previous_entries[field] = last_row[field]
            except IndexError:
                pass

    new_entry = {}
    for field in fields:
        if field in variables:
            value = variables[field]
        else:
            value = input(f"{field}: ")
            if value == '':
                value = previous_entries[field]
        new_entry[field] = value

    with open(master_log_file, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writerow(new_entry)

    # handling date subdirectories
    utc_timestamp = float(new_entry["DATE_TIME.TIME()"])
    date = datetime.datetime.utcfromtimestamp(utc_timestamp)
    date_dir = os.path.join(log_date_dir, str(date.date()))
    os.makedirs(date_dir, exist_ok=True)

    date_log_file = os.path.join(date_dir, f"{date.strftime('%Y%m%d%H%M%S')}_log.csv")
    with open(date_log_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerow(new_entry)
# END DEF

# START PRIMING ROUTINE
ClrScrn()
print("\n#########################################")
print("###     MOTOR PRIMING OPERATION       ###")
print("###  MIT NONLINEAR SOLID MECHANICS    ###")
print("###             JAN 2024              ###")
print("#########################################\n") 

primeMotor()
time.sleep(3)

ClrScrn()
'''
print("\n#########################################")
print("###           CALIBRATE ADC           ###")
print("### CLEAN NEEDLE & TYPE 'Q' TO BEGIN  ###")
print("###  MIT NONLINEAR SOLID MECHANICS    ###")
print("###             JAN 2024              ###")
print("#########################################\n") 
'''
print("\n#########################################")
print("###     NEEDLE INSERTION PROCEDURE    ###")
print("### CLEAN NEEDLE & TYPE 'Q' TO BEGIN  ###")
print("###  MIT NONLINEAR SOLID MECHANICS    ###")
print("###             JAN 2024              ###")
print("#########################################\n")

userStringStop = input()

zero_channel()

# START ZEROING ROUTINE #
print("\n#########################################")
print("###     ZERO-PRESSURE MANUALLY        ###")
print("###  MIT NONLINEAR SOLID MECHANICS    ###")
print("###             JAN 2024              ###")
print("#########################################\n")  

# OPEN AND CHECK PORTS #
results = serial_ports()
rx = GetCnx(results)

quit = False
def quit_program(): 
    global quit 
    quit = True 
# END DEF

quit = False
keyboard.add_hotkey('q', quit_program)

print("######################################################################")
print("###                     MIT ZERO PROCEDURE                         ###")
print("### [q] ONCE SATISFIED, TYPE 'q' TO CONTINUE                       ###")
print("######################################################################\n")

response = input("Would you like to save the data to a CSV file? (Y/N): ").strip().upper()
save_to_csv = response == "Y"
data_pz = []
softWareDelay = []

# Define directory structure
pressure_zeroing_dir = "INJECTION_DATA"
csv_dir = os.path.join(pressure_zeroing_dir, "CSV_FILES")
images_dir = os.path.join(pressure_zeroing_dir, "IMAGES")

# Get current date in YYYY-MM-DD format
current_date = time.strftime("%Y-%m-%d", time.gmtime())
csv_date_dir = os.path.join(csv_dir, current_date)
images_date_dir = os.path.join(images_dir, current_date)

# Create directories if they don't exist
for directory in [pressure_zeroing_dir, csv_dir, images_dir, csv_date_dir, images_date_dir]:
    if not os.path.exists(directory):
        os.makedirs(directory)

time_start_pz = time.time()

while (rx == True):
    if quit:
        break

    current_time_pz = time.time()
    elapsed_time = current_time_pz - time_start_pz

    rawValu = nau7802.read()
    kPa = convert_to_kpa(rawValu, nau7802.gain)

    data_pz.append((elapsed_time, rawValu, kPa))

    #print('\rPRESSURE:', "{:.3f}".format(kPa), end='', flush=True) # post-remove this is approx. 75Hz... need faster(?)
    sys.stdout.write(f'\rPRESSURE: {kPa:.3f}')
    sys.stdout.flush()  # Ensure that the output is flushed to the console immediately    
    #clear_line() # reduces sampling rate to ~66hz... need to make faster or remove; edit: remove

if save_to_csv:
    clear_line()
    csv_filename = input("Please enter the filename (without .csv): ").strip() + ".csv"
    csv_filepath = os.path.join(csv_date_dir, csv_filename)
    with open(csv_filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ELAPSED_TIME", "RAW_ADC_VALUE", "KPA_VALUE"])
        writer.writerows(data_pz)

# Plotting
times, raw_vals, kPas = zip(*data_pz) if data_pz else ([], [], [])
plt.plot(times, kPas, label="Pressure (kPa)")
plt.xlabel("Time Elapsed (s)")
plt.ylabel("Pressure (kPa)")
image_filename = csv_filename.replace(".csv", ".png") if save_to_csv else "pressure_plot.png"
image_filepath = os.path.join(images_date_dir, image_filename)
plt.savefig(image_filepath)
plt.show()

keyboard.unhook_all_hotkeys()  # clear all hotkeys

print("\n#########################################")
print("###      ENTERING MAIN PROCEDURE      ###")
print("###   MIT NONLINEAR SOLID MECHANICS   ###")
print("#########################################\n")
time.sleep(3)

# START MAIN ROUTINE #
ClrScrn()
print("\n#########################################")
print("###         DIGITAL PALPATION         ###")
print("###  MIT NONLINEAR SOLID MECHANICS    ###")
print("###              2023                 ###")
print("#########################################\n")

# OPEN AND CHECK PORTS #
results = serial_ports()
rx = GetCnx(results)

# CHECK CONNECTION TO MICROPUMP #
if(rx == False):
    print("MC2T PORT NOT FOUND.. \n 1. CHECK CONNECTION \n 2. ENSURE 'REMOTE' IS ENABLED ON MC2T CONFIG OPTIONS SCREEN \n 3. ENSURE MC2T IS ON OPERATIONS SCREEN THAT SHOWS PUMPS\n")
    time.sleep(0.002) # formerly 2.0

# MAIN MENU OF MAIN PROGRAM #
while (rx == True):
    print("\n")
    print("######################################################################")
    print("###                     MIT DIGITAL PALPATION                      ###")
    print("### [1] MATLAB_injectFluid_constRadius_001   Procedure             ###")
    print("### [2] MATLAB_injectFluid_constRadius_002   Procedure             ###")
    print("### [3] MATLAB_injectFluid_constRadius_004   Procedure             ###")
    print("### [4] MATLAB_injectFluid_constRadius_008   Procedure             ###")
    print("### [5] MATLAB_injectFluid_constRadius_016   Procedure             ###")
    print("### [6] MATLAB_injectFluid_constRadius_032   Procedure             ###")
    print("### [7] injectFluid_constVolume_customRate Procedure               ###")  
    print("### [8] injectFluid_constVolume_customTime Procedure               ###")
    print("###                                                                ###")
    print("### [0] Troubleshooting Environmental Conditions Procedure         ###")
    print("###                                                                ###")
    print("### [00] CONTACT INFORMATION                                       ###")    
    print("### [X] EXIT THIS PROGRAM                                          ###")
    print("######################################################################\n")
    
    plt.close()
    select = input("SELECT COMMAND [1, 2, 3, ..., or X]:")

    if(select == "1"):
        # DISPLAY/STOP/CANCEL
        ClrScrn()
        print("\n")
        print("######################################################################")
        print("###                       MIT DIGITAL PALPATION                    ###")
        print("###                           FILE GENERATION                      ###")
        print("######################################################################\n")
        
        # Get the filenames
        base = get_base_filename()

        running = True
        time_to_pause = profileTimes[0][rate001]+relaxTime+nIterations_time

        # Establish timing variable
        piecewiseTime = np.array(profileMags[rate001][0]).astype(float) / np.array(profileRates[rate001][0]).astype(float)

        # Incorporate pausing command in withdrawFluid
        injectFluid_constRadius_001mmps()

        functionUsed = "injectFluid_constRadius_001mmps()"

        time_motor_stop_archive = samplingFunction()

        sampling = False
        
        # No need to create a CSV here
        df_raw = pd.DataFrame(rows11)

        # Create dataframe and assign column names directly
        df_with_locs = pd.DataFrame(rows1, columns=['UTC_EPOCH', 'TIME_BETWEEN_SAMPLES', 'TOTAL_TEST_TIME', 'ADC_RAW_GAINED_VALUE', 'KPA'])

        UTC_date = df_with_locs.loc[1, 'UTC_EPOCH'] 
        formatted_UTC_date = convert_timestamp(UTC_date)
        
        outputFilename = csv_constRadius(df_with_locs, base, UTC_runStart, profileRates, profileMags, time_motor_stop_archive, piecewiseTime, rate001)

        plot_and_save1(outputFilename, outputFilename)

        test_log_setup({"DATE_UTC": formatted_UTC_date, "TEMPERATURE": temp_ch1, "DATE_TIME.TIME()": UTC_date, "RELAXATION TIME POST INFUSION": relaxTime, "WAIT TIME POST INSERTION": nIterations_time, "SYRINGE TYPE": chgConfig[1], "NEEDLE GAUGE": needleGauge, "NEEDLE LENGTH": needleLength, "FILENAME": outputFilename, "PRESSURE_SENSOR SERIAL NUMBER": pressureSensor, "FUNCTION USED": functionUsed, "MOTOR OPERATION": chgConfig[9]})

    elif(select == "2"):
        # DISPLAY/STOP/CANCEL
        ClrScrn()
        print("\n")
        print("######################################################################")
        print("###                       MIT DIGITAL PALPATION                    ###")
        print("###                           FILE GENERATION                      ###")
        print("######################################################################\n")
        
        # Get the filenames
        base = get_base_filename()

        running = True
        time_to_pause = profileTimes[0][rate002]+relaxTime+nIterations_time

        # Establish timing variable
        piecewiseTime = np.array(profileMags[rate002][0]).astype(float) / np.array(profileRates[rate002][0]).astype(float)

        #zero_channel()
        # Incorporate pausing command in withdrawFluid
        injectFluid_constRadius_002mmps()

        functionUsed = "injectFluid_constRadius_002mmps()"

        time_motor_stop_archive = samplingFunction()

        sampling = False
        
        # No need to create a CSV here
        df_raw = pd.DataFrame(rows11)

        # Create dataframe and assign column names directly
        df_with_locs = pd.DataFrame(rows1, columns=['UTC_EPOCH', 'TIME_BETWEEN_SAMPLES', 'TOTAL_TEST_TIME', 'ADC_RAW_GAINED_VALUE', 'KPA'])

        UTC_date = df_with_locs.loc[1, 'UTC_EPOCH'] 
        formatted_UTC_date = convert_timestamp(UTC_date)
        
        outputFilename = csv_constRadius(df_with_locs, base, UTC_runStart, profileRates, profileMags, time_motor_stop_archive, piecewiseTime, rate002)

        plot_and_save1(outputFilename, outputFilename)

        test_log_setup({"DATE_UTC": formatted_UTC_date, "TEMPERATURE": temp_ch1, "DATE_TIME.TIME()": UTC_date, "RELAXATION TIME POST INFUSION": relaxTime, "WAIT TIME POST INSERTION": nIterations_time, "SYRINGE TYPE": chgConfig[1], "NEEDLE GAUGE": needleGauge, "NEEDLE LENGTH": needleLength, "FILENAME": outputFilename, "PRESSURE_SENSOR SERIAL NUMBER": pressureSensor, "FUNCTION USED": functionUsed, "MOTOR OPERATION": chgConfig[9]})

    elif(select == "3"):
        # DISPLAY/STOP/CANCEL
        ClrScrn()
        print("\n")
        print("######################################################################")
        print("###                       MIT DIGITAL PALPATION                    ###")
        print("###                           FILE GENERATION                      ###")
        print("######################################################################\n")
        
        # Get the filenames
        base = get_base_filename()

        running = True
        time_to_pause = profileTimes[0][rate004]+relaxTime+nIterations_time

        # Establish timing variable
        piecewiseTime = np.array(profileMags[rate004][0]).astype(float) / np.array(profileRates[rate004][0]).astype(float)

        # Incorporate pausing command in withdrawFluid
        injectFluid_constRadius_004mmps()

        functionUsed = "injectFluid_constRadius_004mmps()"

        time_motor_stop_archive = samplingFunction()

        sampling = False
        
        # No need to create a CSV here
        df_raw = pd.DataFrame(rows11)

        # Create dataframe and assign column names directly
        df_with_locs = pd.DataFrame(rows1, columns=['UTC_EPOCH', 'TIME_BETWEEN_SAMPLES', 'TOTAL_TEST_TIME', 'ADC_RAW_GAINED_VALUE', 'KPA'])

        UTC_date = df_with_locs.loc[1, 'UTC_EPOCH'] 
        formatted_UTC_date = convert_timestamp(UTC_date)
        
        outputFilename = csv_constRadius(df_with_locs, base, UTC_runStart, profileRates, profileMags, time_motor_stop_archive, piecewiseTime, rate004)

        plot_and_save1(outputFilename, outputFilename)

        test_log_setup({"DATE_UTC": formatted_UTC_date, "TEMPERATURE": temp_ch1, "DATE_TIME.TIME()": UTC_date, "RELAXATION TIME POST INFUSION": relaxTime, "WAIT TIME POST INSERTION": nIterations_time, "SYRINGE TYPE": chgConfig[1], "NEEDLE GAUGE": needleGauge, "NEEDLE LENGTH": needleLength, "FILENAME": outputFilename, "PRESSURE_SENSOR SERIAL NUMBER": pressureSensor, "FUNCTION USED": functionUsed, "MOTOR OPERATION": chgConfig[9]})


    elif(select == "4"):
        # DISPLAY/STOP/CANCEL
        ClrScrn()
        print("\n")
        print("######################################################################")
        print("###                       MIT DIGITAL PALPATION                    ###")
        print("###                           FILE GENERATION                      ###")
        print("######################################################################\n")
        
        # Get the filenames
        base = get_base_filename()

        running = True
        time_to_pause = profileTimes[0][rate008]+relaxTime+nIterations_time

        # Establish timing variable
        piecewiseTime = np.array(profileMags[rate008][0]).astype(float) / np.array(profileRates[rate008][0]).astype(float)

        # Incorporate pausing command in withdrawFluid
        injectFluid_constRadius_008mmps()

        functionUsed = "injectFluid_constRadius_008mmps()"

        time_motor_stop_archive = samplingFunction()

        sampling = False
        
        # No need to create a CSV here
        df_raw = pd.DataFrame(rows11)

        # Create dataframe and assign column names directly
        df_with_locs = pd.DataFrame(rows1, columns=['UTC_EPOCH', 'TIME_BETWEEN_SAMPLES', 'TOTAL_TEST_TIME', 'ADC_RAW_GAINED_VALUE', 'KPA'])

        UTC_date = df_with_locs.loc[1, 'UTC_EPOCH'] 
        formatted_UTC_date = convert_timestamp(UTC_date)
        
        outputFilename = csv_constRadius(df_with_locs, base, UTC_runStart, profileRates, profileMags, time_motor_stop_archive, piecewiseTime, rate008)

        plot_and_save1(outputFilename, outputFilename)

        test_log_setup({"DATE_UTC": formatted_UTC_date, "TEMPERATURE": temp_ch1, "DATE_TIME.TIME()": UTC_date, "RELAXATION TIME POST INFUSION": relaxTime, "WAIT TIME POST INSERTION": nIterations_time, "SYRINGE TYPE": chgConfig[1], "NEEDLE GAUGE": needleGauge, "NEEDLE LENGTH": needleLength, "FILENAME": outputFilename, "PRESSURE_SENSOR SERIAL NUMBER": pressureSensor, "FUNCTION USED": functionUsed, "MOTOR OPERATION": chgConfig[9]})


    elif(select == "5"):
        # DISPLAY/STOP/CANCEL
        ClrScrn()
        print("\n")
        print("######################################################################")
        print("###                       MIT DIGITAL PALPATION                    ###")
        print("###                           FILE GENERATION                      ###")
        print("######################################################################\n")
        
        # Get the filenames
        base = get_base_filename()

        running = True
        time_to_pause = profileTimes[0][rate016]+relaxTime+nIterations_time

        # Establish timing variable
        piecewiseTime = np.array(profileMags[rate016][0]).astype(float) / np.array(profileRates[rate016][0]).astype(float)

        # Incorporate pausing command in withdrawFluid
        injectFluid_constRadius_016mmps()

        functionUsed = "injectFluid_constRadius_016mmps()"

        time_motor_stop_archive = samplingFunction()

        sampling = False
        
        # No need to create a CSV here
        df_raw = pd.DataFrame(rows11)

        # Create dataframe and assign column names directly
        df_with_locs = pd.DataFrame(rows1, columns=['UTC_EPOCH', 'TIME_BETWEEN_SAMPLES', 'TOTAL_TEST_TIME', 'ADC_RAW_GAINED_VALUE', 'KPA'])

        UTC_date = df_with_locs.loc[1, 'UTC_EPOCH'] 
        formatted_UTC_date = convert_timestamp(UTC_date)
        
        outputFilename = csv_constRadius(df_with_locs, base, UTC_runStart, profileRates, profileMags, time_motor_stop_archive, piecewiseTime, rate016)

        plot_and_save1(outputFilename, outputFilename)

        test_log_setup({"DATE_UTC": formatted_UTC_date, "TEMPERATURE": temp_ch1, "DATE_TIME.TIME()": UTC_date, "RELAXATION TIME POST INFUSION": relaxTime, "WAIT TIME POST INSERTION": nIterations_time, "SYRINGE TYPE": chgConfig[1], "NEEDLE GAUGE": needleGauge, "NEEDLE LENGTH": needleLength, "FILENAME": outputFilename, "PRESSURE_SENSOR SERIAL NUMBER": pressureSensor, "FUNCTION USED": functionUsed, "MOTOR OPERATION": chgConfig[9]})

    elif(select == "6"):
        # DISPLAY/STOP/CANCEL
        ClrScrn()
        print("\n")
        print("######################################################################")
        print("###                       MIT DIGITAL PALPATION                    ###")
        print("###                           FILE GENERATION                      ###")
        print("######################################################################\n")
        
        # Get the filenames
        base = get_base_filename()

        running = True
        time_to_pause = profileTimes[0][rate032]+relaxTime+nIterations_time

        # Establish timing variable
        piecewiseTime = np.array(profileMags[rate032][0]).astype(float) / np.array(profileRates[rate032][0]).astype(float)

        # Incorporate pausing command in withdrawFluid
        injectFluid_constRadius_032mmps()

        functionUsed = "injectFluid_constRadius_032mmps()"

        time_motor_stop_archive = samplingFunction()

        sampling = False
        
        # No need to create a CSV here
        df_raw = pd.DataFrame(rows11)

        # Create dataframe and assign column names directly
        df_with_locs = pd.DataFrame(rows1, columns=['UTC_EPOCH', 'TIME_BETWEEN_SAMPLES', 'TOTAL_TEST_TIME', 'ADC_RAW_GAINED_VALUE', 'KPA'])

        UTC_date = df_with_locs.loc[1, 'UTC_EPOCH'] 
        formatted_UTC_date = convert_timestamp(UTC_date)
        
        outputFilename = csv_constRadius(df_with_locs, base, UTC_runStart, profileRates, profileMags, time_motor_stop_archive, piecewiseTime, rate032)

        plot_and_save1(outputFilename, outputFilename)

        test_log_setup({"DATE_UTC": formatted_UTC_date, "TEMPERATURE": temp_ch1, "DATE_TIME.TIME()": UTC_date, "RELAXATION TIME POST INFUSION": relaxTime, "WAIT TIME POST INSERTION": nIterations_time, "SYRINGE TYPE": chgConfig[1], "NEEDLE GAUGE": needleGauge, "NEEDLE LENGTH": needleLength, "FILENAME": outputFilename, "PRESSURE_SENSOR SERIAL NUMBER": pressureSensor, "FUNCTION USED": functionUsed, "MOTOR OPERATION": chgConfig[9]})

    elif(select == "7"):
        # DISPLAY/STOP/CANCEL
        ClrScrn()
        print("\n")
        print("######################################################################")
        print("###                       MIT DIGITAL PALPATION                    ###")
        print("###                           FILE GENERATION                      ###")
        print("######################################################################\n")
        
        sampling = True

        # Get the filenames
        base = get_base_filename()

        running = True

        # Incorporate pausing command in withdrawFluid
        injectFluid_constVolume_customRate()

        # Establish timing variable
        piecewiseTime = vol_stop / float(rate_cVolume)
        time_to_pause = piecewiseTime+relaxTime+nIterations_time

        functionUsed = "injectFluid_constVolume_customRate()"

        time_motor_stop_archive = samplingFunction()

        sampling = False
        
        # No need to create a CSV here
        df_raw = pd.DataFrame(rows11)

        # Create dataframe and assign column names directly
        df_with_locs = pd.DataFrame(rows1, columns=['UTC_EPOCH', 'TIME_BETWEEN_SAMPLES', 'TOTAL_TEST_TIME', 'ADC_RAW_GAINED_VALUE', 'KPA'])

        UTC_date = df_with_locs.loc[1, 'UTC_EPOCH'] 
        formatted_UTC_date = convert_timestamp(UTC_date)

        outputFilename = csv_constVolume(df_with_locs, base, UTC_runStart, [rate_cVolume, rate_cVolume], [vol_stop, vol_stop], time_motor_stop_archive, piecewiseTime)

        test_log_setup({"DATE_UTC": formatted_UTC_date, "TEMPERATURE": temp_ch1, "DATE_TIME.TIME()": UTC_date, "RELAXATION TIME POST INFUSION": relaxTime, "WAIT TIME POST INSERTION": nIterations_time, "SYRINGE TYPE": chgConfig[1], "NEEDLE GAUGE": needleGauge, "NEEDLE LENGTH": needleLength, "FILENAME": outputFilename, "PRESSURE_SENSOR SERIAL NUMBER": pressureSensor, "FUNCTION USED": functionUsed, "MOTOR OPERATION": chgConfig[9]})

        plot_and_save1(outputFilename, outputFilename)

    elif(select == "8"):
                # DISPLAY/STOP/CANCEL
        ClrScrn()
        print("\n")
        print("######################################################################")
        print("###                       MIT DIGITAL PALPATION                    ###")
        print("###                           FILE GENERATION                      ###")
        print("######################################################################\n")
        
        sampling = True

        # Get the filenames
        base = get_base_filename()

        running = True

        # Incorporate pausing command in withdrawFluid
        injectFluid_constVolume_customTime()

        # Establish timing variable
        piecewiseTime = vol_stop / rate_cVolume
        time_to_pause = piecewiseTime+relaxTime+nIterations_time

        functionUsed = "injectFluid_constVolume_customTime()"

        time_motor_stop_archive = samplingFunction()

        sampling = False
        
        # No need to create a CSV here
        df_raw = pd.DataFrame(rows11)

        # Create dataframe and assign column names directly
        df_with_locs = pd.DataFrame(rows1, columns=['UTC_EPOCH', 'TIME_BETWEEN_SAMPLES', 'TOTAL_TEST_TIME', 'ADC_RAW_GAINED_VALUE', 'KPA'])

        UTC_date = df_with_locs.loc[1, 'UTC_EPOCH'] 
        formatted_UTC_date = convert_timestamp(UTC_date)

        outputFilename = csv_constVolume(df_with_locs, base, UTC_runStart, [rate_cVolume, rate_cVolume], [vol_stop, vol_stop], time_motor_stop_archive, piecewiseTime)

        test_log_setup({"DATE_UTC": formatted_UTC_date, "TEMPERATURE": temp_ch1, "DATE_TIME.TIME()": UTC_date, "RELAXATION TIME POST INFUSION": relaxTime, "WAIT TIME POST INSERTION": nIterations_time, "SYRINGE TYPE": chgConfig[1], "NEEDLE GAUGE": needleGauge, "NEEDLE LENGTH": needleLength, "FILENAME": outputFilename, "PRESSURE_SENSOR SERIAL NUMBER": pressureSensor, "FUNCTION USED": functionUsed, "MOTOR OPERATION": chgConfig[9]})

        plot_and_save1(outputFilename, outputFilename)
    
    elif(select == "0"):
        print("\n")
        print("######################################################################")
        print("###              Procedure in Process; Please Wait                 ###")
        print("######################################################################\n")
        
        environmentalConditions(samplingRate)

        print("\n")
        print("######################################################################")
        print("###              Procedure is Complete; Thank you                  ###")
        print("######################################################################\n")

    elif(select == "00"):
        print("\n")
        print("######################################################################")
        print("###                        CONTACT INFORMATION                     ###")
        print("###                                                                ###") 
        print("### Brendan M. Unikewicz, <bmu@mit.edu>                            ###")
        print("### Tal Cohen, <talco@mit.edu>                                     ###")
        print("###                                                                ###")
        print("### MIT, Dept. of Mechanical Engineering                           ###")
        print("### MIT, Dept. of Civil & Environmental Engineering                ###")  
        print("###                                                                ###") 
        print("###             AUTO-RETURN TO PREVIOUS SCREEN IN 10 SECONDS.      ###")
        print("######################################################################\n")
        time.sleep(10)

    elif(select == "X"):
        break

# END MAIN ROUTINE
