######################################################################
from coffea.util import save, load
from tqdm import tqdm
import os

import config_merge_data as config

def split_list(path, n_size=config.n_size):

    '''
    This function is used to split a list of files in order to merge its files more efficienctly

    Parameters:
    
    path (str): Path to the directories named RunX, where for 2018 data X = [A, B, C, D]

    n_size (int): Number of chunks we want to split the list. 

    The function return a list with n_size size containing a list with splitted paths and the
    total number of files on the path.

    '''
    import numpy as np

    # List to store final paths
    files = []
    # With statement to open scan more efficiently
    with os.scandir(path) as aux:
        for file in aux:
            # Avoid storing already merged files.
            if file.name.endswith('_merged.coffea'):
                pass
            # Stores all other files finished with .coffea
            if file.name.endswith('.coffea') and (file.stat().st_size != 0):
                files.append(file.path)
    # Saves the total number of files found on the given path.
    n_files = len(files)
    # Used numpy array_split to split the list of filez in a list with n_size length.
    splits = np.array_split(files, n_size)
    
    return splits, n_files
    

def merge(era, file_list, counter):

    '''
    This function is used to sum the accumulator from a given list of coffea files.

    Parameters:
    
    era (str): Name of the era one wants to run the code, where for 2018 data X = [A, B, C, D]

    file_list (list): List of coffea files containing the accumulator one wants to sum.

    counter (int): An interger counter to save the files with a different name,
                   EX: file_0.coffea, file_1.coffea.

    It returns the merged file name.

    '''
    # Loads the first file on the list.
    acc = load(file_list[0])
    # Takes the length of the list
    le_files = len(file_list) 

    # Loop over the list (starts from 1 because we have already called the first file!) using tqdm
    # in order to monitor the loop progress.
    for i in tqdm(range(1, le_files), desc="Processing " + era + '_chunk_' + str(counter)  + '_merged.coffea', unit="files"):
        # Suns the accumulator for each file
        acc += load(file_list[i])
    # Name of the files to be saved
    merged_file = path + era + '/merged_data/' + era + '_chunk_' + str(counter)  + '_merged.coffea'
    # Saves the summed accumulator into a file.
    save(acc, merged_file)    
    
    print(f'File merged and saved as: {merged_file}')
    return merged_file

if __name__ == '__main__':

    '''

    Main function. In the end it will produce n_zise merged .coffea files
    
    '''

    # Takes a list with eras to be runned
    era_list = config.era_list
    # Path to the files
    path = config.path

    # Loop over the era list defined on config_merge_data
    for era in era_list:

        # Complete the main path with the era path
        paths = path + era
       
        # Calls the split_list method
        splited_paths, n_fl = split_list(paths)
        print(f'{n_fl} files to be merged')
        print(f'List divided into {len(splited_paths)} chunks')

        print(f'Reading files from: {path + era}')

        os.system("rm -rf " + paths + '/merged_data/')
        os.system('mkdir -p ' + paths + '/merged_data/')
        print(f'Creating files on:')
        print(paths + '/merged_data/')

        # Loop over the splited_paths. It calls merge function to merge all files
        # for each splite path. It stores the merged file names into a list (If one wants
        # to merge the merged files all it has to do is calling merge function again over this list)
        merged_files = []
        ct = 0
        for array in splited_paths:
            merged_files.append(merge(era, array, ct))
            ct = ct + 1