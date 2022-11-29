''' 
Config filed used to help merge_data script

The folders must be organized in the following way:

Folders with data (root files) must be in the main path, i.e., if the variable path is '/Charm', then
it should be,

/Charm/RunA; /Charm/RunB, etc.

Where RunX are part of the list stored in the variable era_list

In the RunX folder you must have coffea files produced by your condor script.

To run you simply do: python3 merge_data.py

'''

## next version: era_list is going to be a dict with key being Erax and value being n_size!!!!!!!

# List with eras to run
era_list=['RunB', 'RunC', 'RunD', 'RunE', 'RunF',]
era_list=['RunE']

# Number of chunks to divide the file list
n_size = 10

# Path where the files produced by condor are stored.
path = '/afs/cern.ch/work/m/mabarros/public/CMSSW_10_6_12/src/analysis_data/analysis_fit/data/' 

