#!/usr/bin/env python3
import os
import sys
import platform
import argparse

# add dynamic link library or shared object to python env
# attention: match current python version with the version exposing the cpp code
sys_str = platform.system()
path_1 = os.path.abspath(os.path.join(os.getcwd(), '../..'))
if sys_str == 'Windows':
    path_2 = 'lib'
elif sys_str == 'Linux':
    path_2 = 'lib'
else:
    # depend on the system
    path_2 = 'lib'
path = os.path.join(path_1, path_2)
sys.path.append(path)
# change import depending on the project name
import dw_3d_pinched_semi_cylinder_python as test_3d


def run_case(value):
    parser = argparse.ArgumentParser()
    # set case parameters
    parser.add_argument("--restart_step", default=0, type=int)
    parser.add_argument("--end_time", default=20, type=int)
    parser.add_argument("--loading_factor", default=100, type=float)
    case = parser.parse_args()
    # project = test_3d.pinched_semi_cylinder_from_sph_cpp(case.loading_factor)
    project = test_3d.pinched_semi_cylinder_from_sph_cpp(value)
    if project.CmakeTest() == 1:
        project.RunCase()
    else:
        print("check path: ", path)

def net_displacement(file_path, output_file):
    with open(file_path, 'r') as infile:
        lines = infile.readlines()
    
    first_line_values = lines[1].strip().split()
    first_value = float(first_line_values[3])
    
    last_line_values = lines[-1].strip().split()
    last_time = float(last_line_values[0])
    last_value = float(last_line_values[3])
    displacement_z = first_value - last_value
    #displacement_y = float(first_line_values[2]) - float(last_line_values[2]) 
    
    with open(output_file, 'a') as outfile:
        outfile.write("run_time = " + str(last_time) + ", disp = " + str(displacement_z) + '\n')
        #outfile.write("run_time = " + str(last_time) + ", dispy = " + str(displacement_y) + '\n')

    return last_time, displacement_z
    #return last_time, displacement_z, displacement_y

def copy_files(output_folder):
    source_files = ['output/CylinderObserver_Position.dat', 'output/SPHBody_CylinderBody_0000000000.plt', 'output/SPHBody_CylinderBody_0000000001.plt']
    destination_files = [output_folder + 'CylinderObserver_Position.dat', output_folder + 'SPHBody_CylinderBody_0000000000.plt', output_folder + 'SPHBody_CylinderBody_0000000001.plt']

    for index, source_file in enumerate(source_files):
        # Read the content of the source file
        with open(source_file, 'rb') as source:
            content = source.read()
        # Write the content to the destination file
        with open(destination_files[index], 'wb') as destination:
            destination.write(content)

def rename_files(value, output_folder):
    old_file_name = [output_folder + 'CylinderObserver_Position.dat', output_folder + 'SPHBody_CylinderBody_0000000000.plt', output_folder + 'SPHBody_CylinderBody_0000000001.plt']
    new_file_name = [(output_folder + 'CylinderObserver_Position' + value + '.dat'), (output_folder + 'SPHBody_CylinderBody_0000000000_' + value + '.plt'), (output_folder + 'SPHBody_CylinderBody_0000000001_' + value + '.plt')]

    for index, file in enumerate(old_file_name):
        try:
            os.rename(file, new_file_name[index])
            print('File renamed successfully.')
        except OSError:
            print('Error: File rename operation failed.')
     

if __name__ == "__main__":

    file_path = 'output/CylinderObserver_Position.dat'
    output_folder = 'multiple_runs_output/'
    output_file = output_folder + 'Displacements.dat'
    values = [100, 200, 300, 400, 500, 600, 700, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600 ]
    displacement = list(range(len(values)))
    run_time = list(range(len(values)))

    if not os.path.exists(output_folder):
    # Create the folder
        os.makedirs(output_folder)
    else:
        print(f"The folder '{output_folder}' already exists. Skipping creation.")

    with open(output_file, 'w') as outfile:
        outfile.truncate()

    for index,value in enumerate(values):
        with open(output_file, 'a') as outfile:
            outfile.write("loading_factor = " + str(value) + '\n')
        run_case(value)
        run_time[index], displacement = net_displacement(file_path, output_file)
        with open(output_file, 'a') as outfile:
            outfile.write('\n')
        copy_files(output_folder)
        rename_files(str(index), output_folder)
    
    print("------------------------------------------------------")
    print(f"All the cases finished! Files are saved as follows: ")
    for index,value in enumerate(values):
        print(f"For loading_factor = {value} :")
        new_file_name = [(output_folder + 'CylinderObserver_Position' + str(index) + '.dat'), (output_folder + 'SPHBody_CylinderBody_0000000000_' + str(index) + '.plt'), (output_folder + 'SPHBody_CylinderBody_0000000001_' + str(index) + '.plt')]
        print("\t" + new_file_name[0])
        print("\t" + new_file_name[1])
        print("\t" + new_file_name[2])
        
    print(f"Summary of displacements saved in : " + output_file)
    print("\n")
