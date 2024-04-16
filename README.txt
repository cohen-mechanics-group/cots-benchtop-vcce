## Code Purpose

# Designed for use with the World Precision Instruments UMP3, also known as UltraMicroPump II for VCCE.

---

## Hardware Setup

To utilize this code, ensure you have the following:

- **MICROTOUCH**
- **UMP3**
- **Adafruit Trinkey QT2040 - RP2040 USB Key** with Stemma QT
- **Adafruit NAU7802 24-Bit ADC - STEMMA QT / Qwiic**

Additionally, access to two USB-A ports is mandatory.

# Software Setup

- Ensure the updated 'ADC_SITE_PACKAGE/cedargrove_nau7802.py' package is loaded in your local Python site-packages. For example, on my local machine, my ADC package is loaded to: 'C:\Python370\Lib\site-packages'

### Steps:

1. Connect the UMP3 to the MICROTOUCH in channel one.
2. Connect the USB-B to USB-A to calibrate with the MICROTOUCH controller. 
3. Connect the USB-A port to your computer and power on the MICROTOUCH.
4. Navigate: `CONFIGURE >> RESET POS >> ENABLE REMOTE ACCESS`. This completes the World Precision Instruments hardware interface.
5. Link the 24-bit ADC to the PendoTech pressure sensor. Ensure wiring matches colors on the terminal block (e.g., redWire-to-redLabeled_port).
6. Connect the Qwiic wire to the ADC channel one and then to the Trinkey.
7. Plug the Trinkey into another USB-A port on your computer. If correctly set up, a green light should turn on on the ADC.

---

# Code & Profile Setup

There are two main scripts for executing tests with the microPump:

- **VCCE_CONTROL_CODE.py**: Interface with the hardware through serial commands and collect inline pressure data.
  
- **profileGeneration_controlPump.py**: Generates .csv files for pump control.

In the MATLAB script, adjust parameters `astop` and `astart` for final and starting radii of the intended bubble growth. Execute the script to produce the .csv files, then refer to the directory containing `controlPump.py` and run the script. Please adhere to prompts in the `controlPump.py` user interface.

---

# Contact Information

**Massachusetts Institute of Technology**  
Dept. of Mechanical Engineering  
Cohen's Nonlinear Solid Mechanics Group  

- **Development**:  
  - Brendan M Unikewicz, PhD Student, MIT  
    Email: [bmu@mit.edu](mailto:bmu@mit.edu)
    
- Please CC:  
  - Tal Cohen, Professor, MIT  
    Email: [talco@mit.edu](mailto:talco@mit.edu)

---

# Update Log

- **07/06/2023**: 
  1. First commit to the GitHub branch within MIT GitHub Enterprise for lab & research application.
  2. README update with brief instructions, an update log, and datasheets.

