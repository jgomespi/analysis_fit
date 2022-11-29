from coffea.util import load, save
import numpy as np
import os

import uproot3

from tqdm import tqdm

import config_save_ttree as config

def load_accumulator(main_path, era, condition):

    '''
    This function is used load the accumulator for file and sum everything together. There are two conditions,
    if the condition is 'no_trigger' then it will sum all chunks produced with merge_data script. If the condition
    is 'trigger' it will simply take the accumulator produced on TriggerProcessor script.
    
    main_path (str): Path to the directories named RunX, where for 2018 data X = [A, B, C, D]

    era (str): Name of the era one wants to run the code, where for 2018 data X = [A, B, C, D]

    condition (str): It is two options: no_trigger and trigger.

    It returns the desired accumulator and the name of the trigger used (if it exists).

    '''

    # If the condition is without trigger it will sum the accumulators produced with merge_coffea.py
    if condition == 'no_trigger':

        print('Running no_trigger option')
        
        # List to store final paths
        files = []        
        # With statement to open scan more efficiently
        with os.scandir(main_path) as it:
            for file in it:
                # Stores all other files finished with .coffea
                if file.name.endswith('.coffea') and (file.stat().st_size != 0):
                    files.append(file.path)

        # Loads the first file on the list.
        acc = load(files[0])

        # Loop over the list (starts from 1 because we have already called the first file!) using tqdm
        # in order to monitor the loop progress.
        for i in tqdm(files[1:], desc="Loading " + era, unit=" files", total=len(files)-1):
            acc += load(i)
        trigger_name = ''

        return acc, trigger_name

    if condition == 'trigger':

        print('Running trigger option')
        # List to store final paths
        files = []
        # With statement to open scan more efficiently
        with os.scandir(main_path) as it:
            for file in it:
                # Stores all other files finished with .coffea
                if file.name.endswith('.coffea') and (file.stat().st_size != 0):
                    files.append(file.path)

        # Strategy to take the trigger name and store on trigger_name variable.
        # Ex: trigger_name = 'HLT_Dimuon25'
        trigger_name = files[0]
        trigger_name = trigger_name.split('HLT', 1)
        trigger_name = trigger_name[1]
        trigger_name = 'HLT' + trigger_name.replace('.coffea', '')
        trigger_name = trigger_name[:12]    

        # Loads the first file on the list.
        acc = load(files[0])

        # Loop over the list (starts from 1 because we have already called the first file!) using tqdm
        # in order to monitor the loop progress.
        for i in tqdm(files[1:], desc="Loading " + era, unit=" files", total=len(files)-1):
            acc += load(i)
        
        return acc, trigger_name

    

def create_root(accumulator, path_output, era, condition=None):

    '''
    This function is used to produce root ntuple from coffea files.
    
    accumulator (dict_accumulator): Accumulator with the particle information.

    era (str): Name of the era one wants to run the code, where for 2018 data X = [A, B, C, D]

    condition (str): It is two options: no_trigger and trigger.

    It creates a root ntuple with a branch called asso with the following variables:
    Dstar_deltamr, jpsi_pt and jpsi_mass -> All from the associated object.

    '''

    # Taks the accumulator and the trigger name
    acc, trigger_name = accumulator

    # If the option is 'no_trigger' it will take the original accumulator that comes from the
    # files produced via condor
    if condition == 'no_trigger':

        ## Dstar

        # Takes all Dstar data
        all_asso_dstar = acc['DimuDstar']['Dstar']['deltamr'].value
        
        # Takes wrong charge flags
        wrg_chg = acc['DimuDstar']['Dstar']['wrg_chg'].value
        # Takes wrong and right charge Dstars
        dstar_wrong_charge = acc['DimuDstar']['Dstar']['deltamr'].value[wrg_chg]
        dstar_right_charge = acc['DimuDstar']['Dstar']['deltamr'].value[~wrg_chg]
        # Jpsi
        all_asso_jpsi_mass = acc['DimuDstar']['Dimu']['mass'].value
        all_asso_jpsi_pt = acc['DimuDstar']['Dimu']['pt'].value
        
    
    # If the option is 'trigger' it will take the accumulator saved on the TriggerProcessor.py
    if condition == 'trigger':

        ## Dstar

        # Takes all Dstar data
        all_asso_dstar = acc['DimuDstar']['dstar_deltamr'].value
        
        # Takes wrong charge flags
        wrg_chg = acc['DimuDstar']['wrg_chg'].value
        # Takes wrong and right charge Dstars
        dstar_wrong_charge = acc['DimuDstar']['dstar_deltamr'].value[wrg_chg]
        dstar_right_charge = acc['DimuDstar']['dstar_deltamr'].value[~wrg_chg]
        # Jpsi
        all_asso_jpsi_mass = acc['DimuDstar']['jpsi_mass'].value[~wrg_chg]
        all_asso_jpsi_pt = acc['DimuDstar']['jpsi_pt'].value[~wrg_chg]
        all_asso_jpsi_dl = acc['DimuDstar']['jpsi_dl'].value[~wrg_chg]
        all_asso_jpsi_dl_err = acc['DimuDstar']['jpsi_dlErr'].value[~wrg_chg]
        # DimuDstar
        jpsi_dstar_mass = acc['DimuDstar']['dimu_dstar_mass'].value
        jpsi_dstar_deltarap = acc['DimuDstar']['dimu_dstar_deltarap'].value

    # If the variable condition is not give it will raise an exception
    else:
        raise Exception(f' The variable condition is {condition}, which is not valid! You should provide a valide one, either no_trigger or trigger!')
    
    # Creates the root file for reach input.
    if config.cate == '':
        root_name = path_output + era + '_' + trigger_name  + config.cate  + '.root'
    else:
        root_name = path_output + era + '_' + trigger_name  + '_' + config.cate  + '.root'
    with uproot3.recreate(root_name) as ds:
        ds['asso'] = uproot3.newtree({"dstar_mass": "float32",
                                      "jpsi_mass": "float32", 
                                      "jpsi_pt": "float32",
                                      "jpsi_dl": "float32",
                                      "jpsi_dlErr": "float32"})
        ds['asso'].extend({"dstar_mass": dstar_right_charge, 
                           "jpsi_mass": all_asso_jpsi_mass, 
                           "jpsi_pt": all_asso_jpsi_pt,
                           "jpsi_dl": all_asso_jpsi_dl,
                           "jpsi_dlErr": all_asso_jpsi_dl_err})

        ds['obje'] = uproot3.newtree({"jpsi_dstar_mass": "float32",
                                      "jpsi_dstar_deltarap": "float32", })
        ds['obje'].extend({"jpsi_dstar_mass": jpsi_dstar_mass, 
                           "jpsi_dstar_deltarap": jpsi_dstar_deltarap,})

if __name__ == '__main__':

    '''

    Main function. In the end it will create root files with the needed information to
    perform the fits
    
    '''

    # Takes a list with eras to be runned
    era_list= config.era_list
    # Path to the files
    main_path = config.main_path

    # Loop over the era list defined on config_trigger_procesor    
    for era in era_list:

        # Built the path to files with trigger applied
        path_merged_data = main_path  + '/' + era + '/' + 'merged_data/trigger' + '/' + config.cate 
        print('Reading files from:')
        print(main_path  + '/' + era + '/' + 'merged_data/trigger' + '/' + config.cate )

        # Analysis condition: trigger or no_trigger
        condition = config.condition
        
        # Calls load_accumulator function to load the accumulator
        acc = load_accumulator(path_merged_data, era, condition)

        # Path to store the root file
        path_output = config.path_output
        print(f'Creating files on:')
        print(path_output)

        
        # Create the root files
        create_root(acc, path_output, era, condition)

""" import os
main_path = '/home/mabarros/Analysis_2018/OniaOpenCharmRun2ULAna/output/Charmonium_2018/RunA/merged_data'
files = []
with os.scandir(main_path) as it:
    for file in it:
        if file.name.endswith('.coffea') and (file.stat().st_size != 0):
            files.append(file.path)

print(files) """