# Brendan M Unikewicz, PhD Student
# Andre M Pincot, MSc Student
# Tal Cohen, Asc. Professor
# MIT, Dept. Mechanical Engineering
# MIT, Dept. Civil & Environmental Engineering
# Date of Creation: 03/26/2024
# Code Purpose: Plotting any data of interest against one another in 
# Python, quickly, from MATLAB_DATA_PROCESSING data folders and allowing python users
# to do further math operations / processing to data

import os
from scipy.io import loadmat
import matplotlib.pyplot as plt
import numpy as np

def load_mat_files_from_directory(directory):
    """ Load all .mat files from a given directory and classify data into simple matrix or cell array. """
    simple_data = {}  # Dictionary for non-cell data
    cell_data = {}    # Dictionary for cell array data

    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".mat"):
                file_path = os.path.join(root, file)
                contents = loadmat(file_path, simplify_cells=False)  # Do not simplify cells to manage manually

                for key in contents:
                    if key not in ['__header__', '__version__', '__globals__']:
                        item = contents[key]
                        # Check if item is a cell array
                        if isinstance(item, np.ndarray) and item.dtype == 'object' and item.size > 0:
                            # Assuming 1xN cell arrays, extract all elements
                            if item.shape[0] == 1 or item.shape[1] == 1:  # Check if it is a 1xN or Nx1 cell array
                                cell_list = [item.flat[i] for i in range(item.size)]
                                cell_data[key] = cell_list
                            else:
                                print(f"Skipping {key} as it is not a 1xN or Nx1 cell array.")
                        else:
                            # Handle simple matrices and scalars
                            simple_data[key] = item

    return simple_data, cell_data

def quickPlot_one(data, variable_name):
    """ Plot specified variable from the data, checking for both simple and cell array data """
    if variable_name in data:
        plt.figure()
        variable_data = data[variable_name]
        if isinstance(variable_data, list):
            for idx, item in enumerate(variable_data):
                if isinstance(item, np.ndarray):
                    plt.plot(item, label=f'Case {idx + 1}')
                else:
                    plt.plot(item, label=f'Single Case {idx + 1}')
        else:
            plt.plot(variable_data, label='Data')
        
        plt.title(f'{variable_name} Plot')
        plt.xlabel('Index')
        plt.ylabel('Value')
        plt.legend()
        plt.grid(True)
        plt.show()
    else:
        print(f"Data for {variable_name} not found.")

def quickPlot_two(data_x, data_y, label_x, label_y):
    """ Plot two variables against each other from the dataset. """
    if label_x in data_x and label_y in data_y:
        plt.figure()
        # Ensure the lengths of datasets are compatible
        length = min(len(data_x[label_x]), len(data_y[label_y]))
        for idx in range(length):
            x = data_x[label_x][idx]
            y = data_y[label_y][idx]
            if isinstance(x, np.ndarray) and isinstance(y, np.ndarray) and x.shape == y.shape:
                plt.plot(x, y, linestyle='-', label=f'Case {idx + 1}')
            else:
                print(f"Data shapes mismatch or not array for case {idx + 1}")

        plt.title(f'{label_x} vs {label_y}')
        plt.xlabel(label_x)
        plt.ylabel(label_y)
        plt.legend()
        plt.grid(True)
        plt.show()
    else:
        if label_x not in data_x:
            print(f"Data for {label_x} not found.")
        if label_y not in data_y:
            print(f"Data for {label_y} not found.")

def doMath(data, operation):
    """ Apply a mathematical operation to the data. """
    processed_data = []
    if isinstance(data, list):
        for item in data:
            if isinstance(item, np.ndarray):
                processed_data.append(operation(item))
            else:
                print("Item is not a NumPy array.")
    elif isinstance(data, np.ndarray):
        processed_data = operation(data)
    else:
        print("Data format not supported.")
    return processed_data

def main():
    parent_dir = os.path.dirname(os.getcwd())
    mat_data_folder = os.path.join(parent_dir, 'MATLAB_DATA_PROCESSING')

    simple_data_post, cell_data_post = load_mat_files_from_directory(os.path.join(mat_data_folder, 'data_postProcessing'))
    simple_data_pre, cell_data_pre = load_mat_files_from_directory(os.path.join(mat_data_folder, 'data_preProcessing'))

    # Merge simple and cell data from pre-processing and post-processing
    simple_data = {**simple_data_pre, **simple_data_post}
    cell_data = {**cell_data_pre, **cell_data_post}

    print("Simple data keys:", list(simple_data.keys()))
    print("Cell data keys:", list(cell_data.keys()))

    # Quick-Plotting a variable, p_43_1
    if 'p_43_1' in cell_data:
        quickPlot_one(cell_data, 'p_43_1')

    # Quick-Plotting two variables, V_T_43_1 & p_43_1
    if 'p_43_1' in cell_data and 'V_T_43_1' in cell_data:
        quickPlot_two(cell_data, cell_data, 'V_T_43_1', 'p_43_1')

    # Pulling data from dictionary, applying ~math~ operation, and renaming
    norm_p_43_1 = doMath(cell_data['p_43_1'], lambda x: x / np.max(x))

    # Plotting normalized p_43_1 variable: norm_p_43_1
    plt.figure()
    for idx, item in enumerate(norm_p_43_1):
        plt.plot(item, label=f'Normalized Case {idx + 1}')
    plt.title('Normalized p_43_1')
    plt.xlabel('Index')
    plt.ylabel('Normalized Value')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    main()
