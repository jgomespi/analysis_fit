''' 
Config filed used to help merge_data script

The folders must be organized in the following way:

Folders with files must be inside the main path, ie, if the variable path is '/Charm', then
it should be,

/Charm/RunA; /Charm/RunB, etc.

inside RunX folder you must have the coffea files produced by your condor script.

To run you simply do: python3 merge_data

'''

# List with eras to be runned
era_list=['RunB', 'RunC', 'RunD', 'RunE', 'RunF',]
#era_list=['RunE']

# Special name for save it (e.g: cate = 'prompt_jpsi')
cate = ''
# Path where the files produced by condor are stored.
main_path = '/afs/cern.ch/work/m/mabarros/public/CMSSW_10_6_12/src/analysis_data/analysis_fit/data/' 

# Analysis condition: trigger or no_trigger
condition = 'trigger'

# Path to store the root file
path_output =  "/data_root_files/"
