from coffea import processor, hist
from tqdm import tqdm

import awkward as ak
from coffea.util import load, save

from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

#import numpy as np
import config_trigger_processor as config

import os

import numpy as np

def build_p4(acc):

    '''
    This function is used to build the lorentzvector for a given particle accumulator
    
    acc (dict_accumulator): Accumulator with the particle information.

    It returns the four vector for the particle (x, y, z, t)

    '''
    # Uses awkward array zip method to build the four vector    
    p4 = ak.zip({'x': acc['x'].value, 
                 'y': acc['y'].value,
                 'z': acc['z'].value,
                 't': acc['t'].value}, with_name="LorentzVector")

    return p4

def load_accumulator(main_path, era):

    '''
    This function is used load the accumulator for each merged file and sum everything together.
    
    main_path (str): Path to the directories named RunX, where for 2017 data X = [A, B, C, D]

    era (str): Name of the era one wants to run the code, where for 2017 data X = [A, B, C, D]

    It returns the accumulator that comes from the sum of all merged_files)

    '''
    
    files = []
    # It converts the main path to the path where the merged data per era is stored
    era_path = main_path + '/' + era + '/merged_data'

    # With statement to open scan more efficiently
    with os.scandir(era_path) as it:
        for file in it:
            # Stores all files finished with .coffea
            if file.name.endswith('.coffea') and (file.stat().st_size != 0):
                files.append(file.path)
   
   # Loads the first file on the list.
    acc = load(files[0])
    
    # Loop over the list (starts from 1 because we have already called the first file!) using tqdm
    # in order to monitor the loop progress.
    for i in tqdm(files[1:], desc="Loading " + era, unit=" files", total=len(files)-1):
        acc += load(i)
    
    return acc

class TriggerProcessor(processor.ProcessorABC):
   
    '''
        
        A coffea processor class
        
     '''

    def __init__(self, analyzer_name):

        '''
        
        Initialize the processor with a name (can be any name) and an accumulator.
        
        '''
        self.analyzer_name = analyzer_name

        self._accumulator = processor.dict_accumulator({
            'JpsiDstar': processor.dict_accumulator({
                'Jpsi_mass': hist.Hist("Events", hist.Bin("mass", "$M_{\mu^+\mu^-}$ [GeV/$c^2$]", 100, 2.95, 3.25)), 
                'Jpsi_p': hist.Hist("Events", 
                                    hist.Bin("pt", "$p_{T,\mu^+\mu^-}$ [GeV/c]", 100, 0, 100),
                                    hist.Bin("eta", "$\eta_{\mu^+\mu^-}$", 60, -2.5, 2.5),
                                    hist.Bin("phi", "$\phi_{\mu^+\mu^-}$", 70, -3.5, 3.5)),
                'Jpsi_rap': hist.Hist("Events", hist.Bin("rap", "y", 60, -2.5, 2.5)),
                'Jpsi_dlSig': hist.Hist("Events", hist.Bin("dlSig", "dl Significance", 100, -50, 200)),
                'JpsiDstar_deltarap': hist.Hist("Events", hist.Bin("deltarap", "$\Delta y$", 50, -5, 5)),
                'JpsiDstar_mass': hist.Hist("Events", hist.Bin("mass", "$m_{J/\psi D*}$ [$GeV/c^2$]", 50, 0, 100)),
                'Dstar_p': hist.Hist("Events",
                                 hist.Cat("chg", "charge"), 
                                 hist.Bin("pt", "$p_{T,D*}$ [GeV/c]", 100, 0, 50),
                                 hist.Bin("eta", "$\eta_{D*}$", 80, -2.5, 2.5),
                                 hist.Bin("phi", "$\phi_{D*}$", 70, -3.5, 3.5)),
                'Dstar_rap': hist.Hist("Events", 
                                    hist.Cat("chg", "charge"), 
                                    hist.Bin("rap", "y", 60, -2.5, 2.5)),
                'Dstar_deltam': hist.Hist("Events", 
                                        hist.Cat("chg", "charge"), 
                                        hist.Bin("deltam", "$\Delta m$ [$GeV/c^2$]", 50, 0.138, 0.162)),
                'Dstar_deltamr': hist.Hist("Events", 
                                        hist.Cat("chg", "charge"), 
                                        hist.Bin("deltamr", "$\Delta m_{refit}$ [$GeV/c^2$]", 50, 0.138, 0.162)),
            }),
        })
        
    @property
    def accumulator(self):
        return self._accumulator
     
    def process(self, acc, path_output, era, mode, counter):

        '''

        This function is used load the summed accumulator for the desired particles. Then for each particles it creates
        an awkward vector, unflatten it and apply the desired trigger.
        After it save all accumulators with the trigger applied.
        
        acc (dict_accumulator): Summed accumulator with the particle information.

        path_output (str): Name of the output path, by default: main_path/ERA/merged_data/trigger

        It returns the final accumulator.

        '''
        # Creates the output
        output = self.accumulator.identity()
    
        # Opens the accumulator for each object (particles, vertices, trigger...)
        #Muon_lead_acc = acc['Muon_lead']
        #Muon_trail_acc = acc['Muon_trail']
        #Dimu_acc = acc['Dimu']
        #Dstar_acc = acc['Dstar']
        #Dstar_D0_acc = acc['Dstar_D0']
        #Dstar_trk_acc = acc['Dstar_trk']
        DimuDstar_acc = acc['DimuDstar']
        #Primary_vertex_acc = acc['Primary_vertex'] 
        HLT_acc = acc[config.hlt_year]
        DimuDstar_p4 = build_p4(DimuDstar_acc)     

        ## Muon lead collection

        # Creates the pt, eta, phi, m lorentz vector.
        """ Muon_lead = ak.zip({
            'pt' : Muon_lead_acc['pt'].value,
            'eta' : Muon_lead_acc['eta'].value,
            'phi' : Muon_lead_acc['phi'].value,}, with_name='PtEtaPhiMCandidate')
        # Uses unflatten with the number of Dimuon in order to apply trigger correction
        Muon_lead = ak.unflatten(Muon_lead, Muon_lead_acc['nMuon'].value) """

        ## Muon trail collection

        # Creates the pt, eta, phi, m lorentz vector.
        """ Muon_trail = ak.zip({
            'pt' : Muon_trail_acc['pt'].value,
            'eta' : Muon_trail_acc['eta'].value,
            'phi' : Muon_trail_acc['phi'].value,}, with_name='PtEtaPhiMCandidate')
        # Uses unflatten with the number of Dimuon in order to apply trigger correction
        Muon_trail = ak.unflatten(Muon_trail, Muon_trail_acc['nMuon'].value) """

        ## Dimuon collection

        # Creates the pt, eta, phi, m lorentz vector.
        """ Dimu = ak.zip({
            'pt': Dimu_acc['pt'].value,
            'eta': Dimu_acc['eta'].value,
            'phi': Dimu_acc['phi'].value,
            'mass': Dimu_acc['mass'].value,
            'rap': Dimu_acc['rap'].value,
            'dl': Dimu_acc['dl'].value,
            'dlSig': Dimu_acc['dlSig'].value,
            'chi2': Dimu_acc['chi2'].value,
            'cosphi': Dimu_acc['cosphi'].value,
            'is_jpsi' : Dimu_acc['is_jpsi'].value,}, with_name="PtEtaPhiMLorentzVector") 

        # Uses unflatten with the number of Dimuon in order to apply trigger correction
        Dimu = ak.unflatten(Dimu, Dimu_acc['nDimu'].value)"""
 
        ## Dstar collection

        # Creates the pt, eta, phi, m lorentz vector.
        """ Dstar = ak.zip({
            'pt' : Dstar_acc['pt'].value,
            'eta' : Dstar_acc['eta'].value,
            'phi' : Dstar_acc['phi'].value,
            'mass' : Dstar_acc['mass'].value,
            'rap': Dstar_acc['rap'].value,
            'charge' : Dstar_acc['charge'].value,
            'deltam' : Dstar_acc['deltam'].value,
            'deltamr' : Dstar_acc['deltamr'].value,
            'D0cosphi' : Dstar_D0_acc['D0cosphi'].value,
            'D0dlSig' : Dstar_D0_acc['D0dlSig'].value,
            'D0pt': Dstar_D0_acc['D0pt'].value,
            'wrg_chg' : Dstar_acc['wrg_chg'].value}, with_name='PtEtaPhiMCandidate') 

        # Uses unflatten with the number of Dimuon in order to apply trigger correction
        Dstar = ak.unflatten(Dstar, Dstar_acc['nDstar'].value) """

        ## DimuDstar collection

        # Creates the pt, eta, phi, m lorentz vector.
        DimuDstar = ak.zip({
            'jpsi_mass' : DimuDstar_acc['Dimu']['mass'].value,
            'jpsi_pt' : DimuDstar_acc['Dimu']['pt'].value,
            'jpsi_eta' : DimuDstar_acc['Dimu']['eta'].value,
            'jpsi_phi' : DimuDstar_acc['Dimu']['phi'].value,
            'jpsi_rap' : DimuDstar_acc['Dimu']['rap'].value,
            'jpsi_dl' : DimuDstar_acc['Dimu']['dl'].value,
            'jpsi_dlErr' : DimuDstar_acc['Dimu']['dlErr'].value,
            'jpsi_dlsig' : DimuDstar_acc['Dimu']['dlSig'].value,
            'dstar_deltam' : DimuDstar_acc['Dstar']['deltam'].value,
            'dstar_deltamr' : DimuDstar_acc['Dstar']['deltamr'].value,
            'dstar_pt' : DimuDstar_acc['Dstar']['pt'].value,
            'dstar_eta' : DimuDstar_acc['Dstar']['eta'].value,
            'dstar_phi' : DimuDstar_acc['Dstar']['phi'].value,
            'dstar_rap' : DimuDstar_acc['Dstar']['rap'].value,
            'associationProb' : DimuDstar_acc['Dstar']['associationProb'].value,            
            'dimu_dstar_deltarap' : DimuDstar_acc['deltarap'].value,
            'dimu_dstar_mass' : DimuDstar_p4.mass, #is_jpsi & ~wrg_chg & dlSig & dlSig_D0Dstar
            'is_jpsi' : DimuDstar_acc['Dimu']['is_jpsi'].value,
            'wrg_chg': DimuDstar_acc['Dstar']['wrg_chg'].value,}, with_name='PtEtaPhiMCandidate')  
        
        DimuDstar = ak.unflatten(DimuDstar, DimuDstar_acc['nDimuDstar'].value)

        hlt_filter = config.hlt_filter

        print(f"You are running with the trigger(s): {hlt_filter}")
        
        # Loop over trigger list. Applies OR condition to the trigger selection.
        trigger_cut = HLT_acc[hlt_filter[0]].value
        for i in range(0, len(hlt_filter)):
            trigger_cut |= HLT_acc[hlt_filter[i]].value

        # Muon lead collection
        #Muon_lead = Muon_lead[trigger_cut]
            
        # Muon trail collection
        #Muon_trail = Muon_trail[trigger_cut]

        # Jpsi collection
        #Dimu = Dimu[trigger_cut]

       #Dimu = Dimu[(Dimu.pt > 30) & (Dimu.pt < 50)]
        
        # Dstar collection
        #Dstar = Dstar[trigger_cut]

        ## DimuDstar collection

        # Trigger cut
        #DimuDstar = DimuDstar[trigger_cut]

        DimuDstar = DimuDstar[(DimuDstar.jpsi_pt > 25.0) & (DimuDstar.jpsi_pt < 150.0)]
        #DimuDstar = DimuDstar[np.absolute(DimuDstar.jpsi_rap) < 1.2]
        #DimuDstar = DimuDstar[np.absolute(DimuDstar.dstar_rap) < 2.1]
        #print(DimuDstar.jpsi_eta)

        # vtx prob cut 
        #DimuDstar = DimuDstar[DimuDstar.associationProb > 0.1]

        # To fill histograms

        jpsi_mass = ak.flatten(DimuDstar.jpsi_mass)
        jpsi_pt = ak.flatten(DimuDstar.jpsi_pt)
        jpsi_eta = ak.flatten(DimuDstar.jpsi_eta)
        jpsi_phi = ak.flatten(DimuDstar.jpsi_phi)
        jpsi_rap = ak.flatten(DimuDstar.jpsi_rap)
        jpsi_dlsig = ak.flatten(DimuDstar.jpsi_dlsig)
        dstar_deltamr = ak.flatten(DimuDstar.dstar_deltamr)
        dstar_pt = ak.flatten(DimuDstar.dstar_pt)
        dstar_eta = ak.flatten(DimuDstar.dstar_eta)
        dstar_phi = ak.flatten(DimuDstar.dstar_phi)
        dstar_rap = ak.flatten(DimuDstar.dstar_rap)
        jpsi_dstar_deltarap = ak.flatten(DimuDstar.dimu_dstar_deltarap)
        jpsi_dstar_mass = ak.flatten(DimuDstar.dimu_dstar_mass) 

        """muon_lead_acc = processor.dict_accumulator({})
        for var in Muon_lead.fields:
            muon_lead_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(Muon_lead[var])))
        muon_lead_acc["nMuon"] = processor.column_accumulator(ak.to_numpy(ak.num(Muon_lead)))
        output["Muon_lead"] = muon_lead_acc

        muon_trail_acc = processor.dict_accumulator({})
        for var in Muon_trail.fields:
            muon_trail_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(Muon_trail[var])))
        muon_trail_acc["nMuon"] = processor.column_accumulator(ak.to_numpy(ak.num(Muon_trail)))
        output["Muon_trail"] = muon_trail_acc

        dimu_acc = processor.dict_accumulator({})
        for var in Dimu.fields:
            if (var.startswith('t')): continue
            dimu_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(Dimu[var])))
        dimu_acc["nDimu"] = processor.column_accumulator(ak.to_numpy(ak.num(Dimu)))
        output["Dimu"] = dimu_acc"""

        """Dstar_acc = processor.dict_accumulator({})
        Dstar_D0_acc = processor.dict_accumulator({})
        Dstar_trk_acc = processor.dict_accumulator({})
        for var in Dstar.fields:
            if var.startswith('D0'):
                Dstar_D0_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(Dstar[var])))
            elif (var.startswith('K') or var.startswith('pi')):
                Dstar_trk_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(Dstar[var])))
            else:
                Dstar_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(Dstar[var])))
        Dstar_acc["nDstar"] = processor.column_accumulator(ak.to_numpy(ak.num(Dstar)))
        output["Dstar"] = Dstar_acc
        output["Dstar_D0"] = Dstar_D0_acc
        output["Dstar_trk"] = Dstar_trk_acc"""

        DimuDstar_acc = processor.dict_accumulator({})
        DimuDstar_acc['Dimu'] = processor.dict_accumulator({})
        DimuDstar_acc['Dstar'] = processor.dict_accumulator({})
        for var in DimuDstar.fields:
            if (var == '0') or (var =='1'):
                continue
            elif var == 'cand':
                for i0 in DimuDstar[var].fields:
                    DimuDstar_acc[i0] = processor.column_accumulator(ak.to_numpy(ak.flatten(DimuDstar[var][i0])))
            else:
                DimuDstar_acc[var] = processor.column_accumulator(ak.to_numpy(ak.flatten(DimuDstar[var])))

        for var in DimuDstar.slot0.fields:
            DimuDstar_acc['Dimu'][var] = processor.column_accumulator(ak.to_numpy(ak.flatten(DimuDstar.slot0[var])))

        for var in DimuDstar.slot1.fields:
            DimuDstar_acc['Dstar'][var] = processor.column_accumulator(ak.to_numpy(ak.flatten(DimuDstar.slot1[var])))
        DimuDstar_acc['nDimuDstar'] = processor.column_accumulator(ak.to_numpy(ak.num(DimuDstar)))
        output['DimuDstar'] = DimuDstar_acc 

        # Name to the file based on the trigger.
        hlt_name = hlt_filter[0][0:12]

        ## Histograms

        # Jpsi
        output['JpsiDstar']['Jpsi_mass'].fill(mass=jpsi_mass)
        output['JpsiDstar']['Jpsi_p'].fill(pt=jpsi_pt,
                                 eta=jpsi_eta,
                                 phi=jpsi_phi)
        output['JpsiDstar']['Jpsi_rap'].fill(rap=jpsi_rap)
        output['JpsiDstar']['Jpsi_dlSig'].fill(dlSig=jpsi_dlsig)

        # JpsiDstar
        output['JpsiDstar']['Jpsi_mass'].fill(mass=jpsi_mass)
        output['JpsiDstar']['Jpsi_p'].fill(pt=jpsi_pt,
                                           eta=jpsi_eta,
                                           phi=jpsi_phi)
        output['JpsiDstar']['Jpsi_rap'].fill(rap=jpsi_rap)

        output['JpsiDstar']['Dstar_deltamr'].fill(chg='right charge', deltamr=dstar_deltamr)
        output['JpsiDstar']['Dstar_p'].fill(chg='right charge',
                                            pt=dstar_pt,
                                            eta=dstar_eta,
                                            phi=dstar_phi)
        
        output['JpsiDstar']['Dstar_rap'].fill(chg='right charge', rap=dstar_rap)

        output['JpsiDstar']['JpsiDstar_deltarap'].fill(deltarap=jpsi_dstar_deltarap)
        output['JpsiDstar']['JpsiDstar_mass'].fill(mass=jpsi_dstar_mass)

        if mode == 'sum':
            print('Saving accumulator...')
            save(output, path_output + '/' + era + '_' + hlt_name + '.coffea')

        elif mode == 'several':
            print(f'Saving file: ')
            if config.cate == '':
                print(era + '_' + hlt_name + '_' + config.cate + str(counter) + '.coffea')
                save(output, path_output + '/' + era + '_' + hlt_name + config.cate + '_' + str(counter) + '.coffea')
            else:
                print(era + '_' + hlt_name + '_' + config.cate + '_' + str(counter) + '.coffea')
                save(output, path_output + '/' + era + '_' + hlt_name + '_' + config.cate + '_' + str(counter) + '.coffea')

        return output

    def postprocess(self, accumulator):
        return accumulator     

if __name__ == '__main__':

    '''

    Main function. In the end it will apply the trigger and save the accumulator.
    
    '''

    # Mode of running: This is a special feature created for when the sum of the accumulators is very large.
    # If mode = 'several' it will apply the trigger for each file separatedely
    # If mode = 'sum' it will apply the trigger for the summed accumulator.
    mode = config.mode

    # Takes a list with eras to be runned
    era_list = config.era_list
    # Path to the files
    main_path = config.main_path

    # Instantiates the object Trigger Processor.
    p = TriggerProcessor('ex')

    # Loop over the era list defined on config_trigger_procesor
    for era in era_list:

        print(f'Reading files from: ')
        print(main_path + '/' + era + '/merged_data')
        
        # Path to the output .coffea files
        path_output = main_path + '/' + era + '/merged_data/trigger' + '/' + config.cate 
        os.system("rm -rf " + path_output)
        os.system("mkdir -p " + path_output)
        print(f'Creating files on:')
        print(path_output)

        if mode == 'sum':
            # Calls load_accumulator function to load the summed accumulator
            acc = load_accumulator(main_path, era)
            
            
            # Calls function process on the trigger processor object.
            p.process(acc, path_output, era, mode, None)

        elif mode == 'several':

            path_mode_several = main_path + '/' + era + '/merged_data'
            files = []

            # With statement to open scan more efficiently
            with os.scandir(path_mode_several) as it:
                for file in it:
                    # Stores all files finished with .coffea
                    if file.name.endswith('.coffea') and (file.stat().st_size != 0):
                        files.append(file.path)
            
            c = 0
            for f in files:
                acc = load(f)
                p.process(acc, path_output, era, mode, c)
                c = c +1
        
        

