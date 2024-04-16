##########################################
## Code Purpose: 
## 1. install packages that bmu@mit.edu
## has on their 3.7.0 version of Python
## for VCCE operations
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

import os
import sys
import subprocess

def install(package, version):
    subprocess.call([sys.executable, "-m", "pip", "install", f"{package}=={version}"])

packages = {
    "pyserial": "3.5",
    "numpy": "1.21.6",
    "matplotlib": "3.5.3",
    "pandas": "1.1.5",
    "hidapi": "0.14.0",
    "keyboard": "0.13.5",
    "cedargrove-nau7802": "2.0.0"
}

for package, version in packages.items():
    install(package, version)

