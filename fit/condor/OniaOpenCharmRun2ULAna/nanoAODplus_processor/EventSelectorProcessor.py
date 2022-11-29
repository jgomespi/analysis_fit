import awkward as ak
import numpy as np
import coffea.processor as processor
from coffea.util import save

from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

import random

from tools.collections import *

D0_PDG_MASS = 1.864

def association(cand1, cand2):
    ''' Function for association of the particles. The cuts that operates on all of them and 
    computation of quantities can go here. individual cuts can go on the main processing'''

    asso = ak.cartesian([cand1, cand2])    

    cand1 = ak.zip({
            'pt': asso.slot0.pt,
            'eta': asso.slot0.eta,
            'phi': asso.slot0.phi,
            'mass': asso.slot0.mass,
            'charge': asso.slot0.charge}, with_name="PtEtaPhiMCandidate")

    cand2 = ak.zip({
            'pt': asso.slot1.pt,
            'eta': asso.slot1.eta,
            'phi': asso.slot1.phi,
            'mass': asso.slot1.mass,
            'charge': asso.slot1.charge}, with_name="PtEtaPhiMCandidate")

    asso['deltarap'] = asso.slot0.rap - asso.slot1.rap
    asso['cand'] = cand1 + cand2
    
    return asso

class EventSelectorProcessor(processor.ProcessorABC):
    def __init__(self, analyzer_name):
        self.analyzer_name = analyzer_name

        self._accumulator = processor.dict_accumulator({
            'cutflow': processor.defaultdict_accumulator(int),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        ############### Cuts
        # Dimu cuts: charge = 0, mass cuts and chi2...
        # test if there is any events in the file
        if len(events) == 0:
            return output

        ############### Get the primary vertices ############### 
        Primary_vertex = ak.zip({**get_vars_dict(events, primary_vertex_cols)})
              
        ############### Get All the interesting candidates from NTuples
        Dimu = ak.zip({**get_vars_dict(events, dimu_cols)}, with_name="PtEtaPhiMCandidate")
        Muon = ak.zip({**get_vars_dict(events, muon_cols)}, with_name="PtEtaPhiMCandidate")
        D0 = ak.zip({'mass': events.D0_mass12, **get_vars_dict(events, d0_cols)}, with_name="PtEtaPhiMCandidate")
        Dstar = ak.zip({'mass': (events.Dstar_D0mass + events.Dstar_deltamr),
                        'charge': events.Dstar_pischg,
                        **get_vars_dict(events, dstar_cols)}, 
                        with_name="PtEtaPhiMCandidate")
        # Triggers for 2016 charmonium
        hlt_char_2016 = ak.zip({**get_vars_dict(events, hlt_cols_charm_2016)})

        ##### Trigger cut

        # Activate trigger
        hlt = False
        
        # HLT to be used
    
        hlt_filter = ['HLT_Dimuon20_Jpsi'] #['HLT_Dimuon10_Jpsi_Barrel', 'HLT_Dimuon16_Jpsi' , 'HLT_Dimuon20_Jpsi']
        
        # Trigger choice
        if hlt:
            #print(f"You are running with the trigger: {hlt_filter}")

            trigger_cut1 = hlt_char_2016[hlt_filter[0]]
            #trigger_cut2 = hlt_char_2016[hlt_filter[1]]
            #trigger_cut3 = hlt_char_2016[hlt_filter[2]]

            Dimu = Dimu[(trigger_cut1)]
            Muon = Muon[(trigger_cut1)]
            Dstar = Dstar[(trigger_cut1)]
            hlt_char_2016 = hlt_char_2016[(trigger_cut1)]
        
        if not hlt:
            #print("You are not running with trigger")
            # Assign 1 to all events.
            trigger_cut = np.ones(len(Dimu), dtype=bool)

            Dimu = Dimu[trigger_cut]
            Muon = Muon[trigger_cut]
            Dstar = Dstar[trigger_cut]

        ############### Dimu cuts charge = 0, mass cuts and chi2...
        Dimu = ak.mask(Dimu, Dimu.charge == 0)

        Dimu = ak.mask(Dimu, ((Dimu.mass > 2.95) & (Dimu.mass < 3.25)))

        ############### Get the Muons from Dimu, for cuts in their params
        Muon = ak.zip({'0': Muon[Dimu.t1muIdx], '1': Muon[Dimu.t2muIdx]})

        # SoftId and Global Muon cuts
        soft_id = (Muon.slot0.softId > 0) & (Muon.slot1.softId > 0)
        Dimu = ak.mask(Dimu, soft_id)
        Muon = ak.mask(Muon, soft_id)

        #!!!!!!!!!!!!!!!!! We decided to remove the global cuts !!!!!!!! #

        #global_muon = (Muon.slot0.isGlobal > 0) & (Muon.slot1.isGlobal > 0)
        #Dimu = Dimu[global_muon]
        #Muon = Muon[global_muon]
        #output['cutflow']['Dimu muon global'] += ak.sum(ak.num(Dimu))

        ## pT and eta/rapidity cuts
    
        # Muon pT
        muon_pt_cut = (Muon.slot0.pt > 3) & (Muon.slot1.pt > 3)
        Dimu = ak.mask(Dimu, muon_pt_cut)
        Muon = ak.mask(Muon, muon_pt_cut)

        # Muon eta 
        muon_eta_cut = (np.absolute(Muon.slot0.eta) <= 2.4) & (np.absolute(Muon.slot1.eta) <= 2.4)
        Dimu = ak.mask(Dimu, muon_eta_cut)
        Muon = ak.mask(Muon, muon_eta_cut)

        # Dimuon pT
        dimu_pt_cut = (Dimu.pt > 20) & (Dimu.pt < 150)
        Dimu = ak.mask(Dimu, dimu_pt_cut)
        Muon = ak.mask(Muon, dimu_pt_cut)

        # Dimuon rapidity
        dimu_rap_cut = (np.absolute(Dimu.rap) <= 2.5)
        Dimu = ak.mask(Dimu, dimu_rap_cut)
        Muon = ak.mask(Muon, dimu_rap_cut)

        Dimu['is_ups'] = (Dimu.mass > 8.5) & (Dimu.mass < 11.5)
        Dimu['is_jpsi'] = (Dimu.mass > 2.9) & (Dimu.mass < 3.3)
        Dimu['is_psi'] = (Dimu.mass > 3.35) & (Dimu.mass < 4.05)

        ############### Cuts for Dstar

        # trks cuts
        #Dstar = Dstar[trigger_cut]
        Dstar = Dstar[~Dstar.hasMuon]

        Dstar = Dstar[Dstar.Kchg != Dstar.pichg]

        Dstar = Dstar[(Dstar.pt > 4) & (Dstar.pt < 100)]

        Dstar = Dstar[np.absolute(Dstar.rap) < 2.5]

        Dstar = Dstar[(Dstar.Kpt > 1.6) & (Dstar.pipt > 1.6)]

        Dstar = Dstar[(Dstar.Kchindof < 2.5) & (Dstar.pichindof < 2.5)]

        Dstar = Dstar[(Dstar.KnValid > 4) & (Dstar.pinValid > 4) & (Dstar.KnPix > 1) & (Dstar.pinPix > 1)]

        Dstar = Dstar[(Dstar.Kdxy < 0.5/(2 * np.arctan(np.exp(-Dstar.Keta)))) & (Dstar.pidxy < 0.5/(2 * np.arctan(np.exp(-Dstar.pieta))))]

        Dstar = Dstar[(Dstar.Kdz < 0.5/(2 * np.arctan(np.exp(-Dstar.Keta)))) & (Dstar.pidz < 0.5)/(2 * np.arctan(np.exp(-Dstar.pieta)))]

        # pis cuts
        Dstar = Dstar[Dstar.pispt > 0.3]

        Dstar = Dstar[Dstar.pischindof < 3]

        Dstar = Dstar[Dstar.pisnValid > 2]

        # D0 of Dstar cuts
        Dstar = Dstar[Dstar.D0cosphi > 0.99]

        Dstar = Dstar[(Dstar.D0mass < D0_PDG_MASS + 0.028) & (Dstar.D0mass > D0_PDG_MASS - 0.028)]

        Dstar = Dstar[Dstar.D0pt > 4]

        Dstar = Dstar[Dstar.D0dlSig > 2.5]

        Dstar['wrg_chg'] = (Dstar.Kchg == Dstar.pichg)

        ############### Dimu + OpenCharm associations

        DimuDstar = association(Dimu, Dstar) 
        DimuDstar = DimuDstar[DimuDstar.slot1.associationProb > 0.01]
        DimuDstar = DimuDstar[ak.fill_none(DimuDstar.slot0.pt, -1) > -1]

        ############### Leading and Trailing muon separation Gen_particles
        leading_mu = (Muon.slot0.pt > Muon.slot1.pt)
        Muon_lead = ak.where(leading_mu, Muon.slot0, Muon.slot1)
        Muon_trail = ak.where(~leading_mu, Muon.slot0, Muon.slot1)

        ############### Create the accumulators to save output

        ## Trigger accumulator

        # 2016 triggers
        trigger_2016_acc = processor.dict_accumulator({})
        for var in hlt_char_2016.fields:
            trigger_2016_acc[var] = processor.column_accumulator(ak.to_numpy(hlt_char_2016[var]))
        output["HLT_2016"] = trigger_2016_acc

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

        file_hash = str(random.getrandbits(128)) + str(len(events))
        save(output, "output/" + self.analyzer_name + "/" + self.analyzer_name + "_" + file_hash + ".coffea")

        # return dummy accumulator
        return processor.dict_accumulator({
                'cutflow': output['cutflow'],
                #'cut_studied': output['cut_studied']
        })
        
    def postprocess(self, accumulator):
        return accumulator
