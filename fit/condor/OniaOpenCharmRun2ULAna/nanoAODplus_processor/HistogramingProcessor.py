from coffea import processor, hist

import awkward as ak
from coffea.util import load

import numpy as np

def build_p4(acc):
    p4 = ak.zip({'x': acc['x'].value, 
                 'y': acc['y'].value,
                 'z': acc['z'].value,
                 't': acc['t'].value}, with_name="LorentzVector")

    return p4

class HistogramingProcessor(processor.ProcessorABC):
    def __init__(self, analyzer_name):
        self.analyzer_name = analyzer_name
        
        self._accumulator = processor.dict_accumulator({
            #'PV_npvs': hist.Hist("Events", hist.Bin("npvs", "Num of associated tracks", 50, 0, 100)), 
            'Muon_lead_p': hist.Hist("Events", 
                                   hist.Bin("pt", "$p_{T,\mu}$ [GeV]", 100, 0, 100),
                                   hist.Bin("eta", "$\eta_{\mu}$", 60, -2.5, 2.5),
                                   hist.Bin("phi", "$\phi_{\mu}$", 70, -3.5, 3.5)),
            'Muon_trail_p': hist.Hist("Events", 
                                       hist.Bin("pt", "$p_{T,\mu}$ [GeV]", 100, 0, 50),
                                       hist.Bin("eta", "$\eta_{\mu}$", 60, -2.5, 2.5),
                                       hist.Bin("phi", "$\phi_{\mu}$", 70, -3.5, 3.5)),
            'Jpsi_mass': hist.Hist("Events", hist.Bin("mass", "$M_{\mu^+\mu^-}$ [GeV]", 100, 2.95, 3.25)),
            'Jpsi_p': hist.Hist("Events", 
                                   hist.Bin("pt", "$p_{T,\mu^+\mu^-}$ [GeV]", 100, 0, 100),
                                   hist.Bin("eta", "$\eta_{\mu^+\mu^-}$", 60, -2.5, 2.5),
                                   hist.Bin("phi", "$\phi_{\mu^+\mu^-}$", 70, -3.5, 3.5)),
            'Jpsi_rap': hist.Hist("Events", hist.Bin("rap", "y", 60, -2.5, 2.5)),
            'Jpsi_dl': hist.Hist("Events", hist.Bin("dl", "dl", 100, -1.5, 1.5)),
            'Jpsi_dlSig': hist.Hist("Events", hist.Bin("dlSig", "dl Significance", 100, -20, 50)),
            'Jpsi_chi2': hist.Hist("Events", hist.Bin("chi2", r"$\chi^2$", 50, 0, 5)),
            'Jpsi_cosphi': hist.Hist("Events", hist.Bin("cosphi", r"$cos(\alpha)$", 50, -1, 1)),
            'D0_trk_dz': hist.Hist("Events", hist.Bin("dz", "dz", 100, -1, 1)),
            'Dstar_p': hist.Hist("Events",
                                 hist.Cat("chg", "charge"), 
                                 hist.Bin("pt", "$p_{T,D*}$ [GeV]", 100, 0, 50),
                                 hist.Bin("eta", "$\eta_{D*}$", 80, -2.5, 2.5),
                                 hist.Bin("phi", "$\phi_{D*}$", 70, -3.5, 3.5)),
            'Dstar_rap': hist.Hist("Events", 
                                   hist.Cat("chg", "charge"), 
                                   hist.Bin("rap", "y", 60, -2.5, 2.5)),
            'Dstar_deltam': hist.Hist("Events", 
                                      hist.Cat("chg", "charge"), 
                                      hist.Bin("deltam", "$\Delta m$ [GeV]", 50, 0.138, 0.162)),
            'Dstar_deltamr': hist.Hist("Events", 
                                       hist.Cat("chg", "charge"), 
                                       hist.Bin("deltamr", "$\Delta m_{refit}$ [GeV]", 50, 0.138, 0.162)),
            'Dstar_D0cosphi' : hist.Hist("Events",
                                         hist.Cat("chg", "charge"),
                                         hist.Bin("cosphi", r"$cos(\alpha)$", 50, -1, 1)),
            'Dstar_D0dlSig' : hist.Hist("Events",
                                         hist.Cat("chg", "charge"),
                                         hist.Bin("dlSig", r"$D^0$ from $D*$ dl Sig", 100, -20, 50)),
            'Dstar_D0pt' : hist.Hist("Events",
                                         hist.Cat("chg", "charge"),
                                         hist.Bin("pt", r"$D^0$ from $D*$ $p_T$", 100, 0, 50)),
            'Dstar_K_p': hist.Hist("Events", 
                                   hist.Bin("pt", "$p_{T,D* K}$ [GeV]", 100, 0, 30),
                                   hist.Bin("eta", "$\eta_{D* K}$", 60, -2.5, 2.5),
                                   hist.Bin("phi", "$\phi_{D* K}$", 70, -3.5, 3.5)),
            'Dstar_K_chindof': hist.Hist("Events", hist.Bin("chindof", r"$\chi^2/ndof$", 50, 0, 2.5)),
            'Dstar_K_nValid': hist.Hist("Events", hist.Bin("nValid", "# of Tracker Hits", 40, -0.5, 39.5)),
            'Dstar_K_nPix': hist.Hist("Events", hist.Bin("nPix", "# of Pixel Hits", 15, -0.5, 14.5)),
            'Dstar_K_dxy': hist.Hist("Events", hist.Bin("dxy", "dxy", 100, 0, 0.1)),
            'Dstar_K_dz': hist.Hist("Events", hist.Bin("dz", "dz", 100, -0.2, 0.2)),
            'Dstar_K_pt_eta': hist.Hist("Events",
                                        hist.Bin("pt", "$p_{T,D* K}$ [GeV]", 100, 0, 10),
                                        hist.Bin("eta", "$\eta_{D* K}$", 60, -2.5, 2.5)),
            'Dstar_pi_p': hist.Hist("Events", 
                                    hist.Bin("pt", "$p_{T,D* \pi}$ [GeV]", 100, 0, 30),
                                    hist.Bin("eta", "$\eta_{D* \pi}$", 60, -2.5, 2.5),
                                    hist.Bin("phi", "$\phi_{D* \pi}$", 70, -3.5, 3.5)),
            'Dstar_pi_chindof': hist.Hist("Events", hist.Bin("chindof", r"$\chi^2/ndof$", 50, 0, 2.5)),
            'Dstar_pi_nValid': hist.Hist("Events", hist.Bin("nValid", "# of Tracker Hits", 40, -0.5, 39.5)),
            'Dstar_pi_nPix': hist.Hist("Events", hist.Bin("nPix", "# of Pixel Hits", 15, -0.5, 14.5)),
            'Dstar_pi_dxy': hist.Hist("Events", hist.Bin("dxy", "dxy", 100, 0, 0.1)),
            'Dstar_pi_dz': hist.Hist("Events", hist.Bin("dz", "dz", 100, -0.1, 0.1)),
            'Dstar_pi_pt_eta': hist.Hist("Events",
                                          hist.Bin("pt", "$p_{T,D* \pi}$ [GeV]", 100, 0, 10),
                                          hist.Bin("eta", "$\eta_{D* \pi}$", 60, -2.5, 2.5)),
            'Dstar_pis_p': hist.Hist("Events", 
                                     hist.Bin("pt", "$p_{T,D* \pi_s}$ [GeV]", 100, 0, 20),
                                     hist.Bin("eta", "$\eta_{D* \pi_s}$", 60, -2.5, 2.5),
                                     hist.Bin("phi", "$\phi_{D* \pi_s}$", 70, -3.5, 3.5)),
            'Dstar_pis_chindof': hist.Hist("Events", hist.Bin("chindof", r"$\chi^2/ndof$", 50, 0, 5)),
            'Dstar_pis_nValid': hist.Hist("Events", hist.Bin("nValid", "# of Tracker Hits", 40, -0.5, 39.5)),
            'Dstar_pis_nPix': hist.Hist("Events", hist.Bin("nPix", "# of Pixel Hits", 15, -0.5, 14.5)),
            'Dstar_pis_dxy': hist.Hist("Events", hist.Bin("dxy", "dxy", 100, 0, 0.2)),
            'Dstar_pis_dz': hist.Hist("Events", hist.Bin("dz", "dz", 100, -2, 2)),
            
            'JpsiDstar': processor.dict_accumulator({
                'Jpsi_mass': hist.Hist("Events", hist.Bin("mass", "$M_{\mu^+\mu^-}$ [GeV]", 100, 2.95, 3.25)), 
                'Jpsi_p': hist.Hist("Events", 
                                    hist.Bin("pt", "$p_{T,\mu^+\mu^-}$ [GeV]", 100, 0, 100),
                                    hist.Bin("eta", "$\eta_{\mu^+\mu^-}$", 60, -2.5, 2.5),
                                    hist.Bin("phi", "$\phi_{\mu^+\mu^-}$", 70, -3.5, 3.5)),
                'Jpsi_rap': hist.Hist("Events", hist.Bin("rap", "y", 60, -2.5, 2.5)),
                'Dstar_p': hist.Hist("Events",
                                 hist.Cat("chg", "charge"), 
                                 hist.Bin("pt", "$p_{T,D*}$ [GeV]", 100, 0, 50),
                                 hist.Bin("eta", "$\eta_{D*}$", 80, -2.5, 2.5),
                                 hist.Bin("phi", "$\phi_{D*}$", 70, -3.5, 3.5)),
                'Dstar_rap': hist.Hist("Events", 
                                    hist.Cat("chg", "charge"), 
                                    hist.Bin("rap", "y", 60, -2.5, 2.5)),
                'Dstar_deltam': hist.Hist("Events", 
                                        hist.Cat("chg", "charge"), 
                                        hist.Bin("deltam", "$\Delta m$ [GeV]", 50, 0.138, 0.162)),
                'Dstar_deltamr': hist.Hist("Events", 
                                        hist.Cat("chg", "charge"), 
                                        hist.Bin("deltamr", "$\Delta m_{refit}$ [GeV]", 50, 0.138, 0.162)),
                'JpsiDstar_deltarap': hist.Hist("Events", hist.Bin("deltarap", "$\Delta y$", 50, -5, 5)),
                'JpsiDstar_mass': hist.Hist("Events", hist.Bin("mass", "$m_{J/\psi D*}$ [GeV]", 100, 0, 100)),
                'JpsiDstar_associationIdx': hist.Hist("Events", hist.Bin("associationIdx", "Association idx", 100, -100, 10)),
                'JpsiDstar_vtx_prob': hist.Hist("Events", hist.Bin("vtx_prob", "$vtx prob", 100, 0, 3)),
                'JpsiDstar_vtx_chi2': hist.Hist("Events", hist.Bin("vtx_chi2", "$vtx chi2", 80, 0, 10)),
            }),
            
        })
        
    @property
    def accumulator(self):
        return self._accumulator
     
    def process(self, file):
        output = self.accumulator.identity()
        acc = load(file)

    
        Muon_lead_acc = acc['Muon_lead']
        Muon_trail_acc = acc['Muon_trail']
        Dimu_acc = acc['Dimu']
        Dstar_acc = acc['Dstar']
        Dstar_D0_acc = acc['Dstar_D0']
        Dstar_trk_acc = acc['Dstar_trk']
        DimuDstar_acc = acc['DimuDstar']
        Primary_vertex_acc = acc['Primary_vertex'] 
        HLT_2016_acc = acc['HLT_2016']
        DimuDstar_p4 = build_p4(DimuDstar_acc)

        ######################## Cuts ######################## 

        is_jpsi = DimuDstar_acc['Dimu']['is_jpsi'].value
        wrg_chg = DimuDstar_acc['Dstar']['wrg_chg'].value          

        ##### Creates coffea lorentz vector to apply trigger on the data #####

        ## Muon lead collection

        # Creates the pt, eta, phi, m lorentz vector.
        Muon_lead = ak.zip({
            'pt' : Muon_lead_acc['pt'].value,
            'eta' : Muon_lead_acc['eta'].value,
            'phi' : Muon_lead_acc['phi'].value,}, with_name='PtEtaPhiMCandidate')
        # Uses unflatten with the number of Dimuon in order to apply trigger correction
        Muon_lead = ak.unflatten(Muon_lead, Muon_lead_acc['nMuon'].value)

        ## Muon trail collection

        # Creates the pt, eta, phi, m lorentz vector.
        Muon_trail = ak.zip({
            'pt' : Muon_trail_acc['pt'].value,
            'eta' : Muon_trail_acc['eta'].value,
            'phi' : Muon_trail_acc['phi'].value,}, with_name='PtEtaPhiMCandidate')
        # Uses unflatten with the number of Dimuon in order to apply trigger correction
        Muon_trail = ak.unflatten(Muon_trail, Muon_trail_acc['nMuon'].value)

        ## Dimuon collection

        # Creates the pt, eta, phi, m lorentz vector.
        Dimu = ak.zip({
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
        Dimu = ak.unflatten(Dimu, Dimu_acc['nDimu'].value)

        ## Dstar collection

        # Creates the pt, eta, phi, m lorentz vector.
        Dstar = ak.zip({
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
        Dstar = ak.unflatten(Dstar, Dstar_acc['nDstar'].value)

        ## DimuDstar collection

        # Creates the pt, eta, phi, m lorentz vector.
        DimuDstar = ak.zip({
            'jpsi_mass' : DimuDstar_acc['Dimu']['mass'].value,
            'jpsi_pt' : DimuDstar_acc['Dimu']['pt'].value,
            'jpsi_eta' : DimuDstar_acc['Dimu']['eta'].value,
            'jpsi_phi' : DimuDstar_acc['Dimu']['phi'].value,
            'jpsi_rap' : DimuDstar_acc['Dimu']['rap'].value,
            'dstar_deltam' : DimuDstar_acc['Dstar']['deltam'].value,
            'dstar_deltamr' : DimuDstar_acc['Dstar']['deltamr'].value,
            'dstar_pt' : DimuDstar_acc['Dstar']['pt'].value,
            'dstar_eta' : DimuDstar_acc['Dstar']['eta'].value,
            'dstar_phi' : DimuDstar_acc['Dstar']['phi'].value,
            'dstar_rap' : DimuDstar_acc['Dstar']['rap'].value,
            'dimu_dstar_deltarap' : DimuDstar_acc['deltarap'].value,
            'dimu_dstar_associationIdx' : DimuDstar_acc['Dstar']['associationIdx'].value,
            'dimu_dstar_vtx_prob' : DimuDstar_acc['Dstar']['associationProb'].value,
            'dimu_dstar_vtx_chi2' : DimuDstar_acc['Dstar']['associationchi2'].value,
            'dimu_dstar_mass' : DimuDstar_p4.mass, #is_jpsi & ~wrg_chg
            'is_jpsi' : DimuDstar_acc['Dimu']['is_jpsi'].value,
            'wrg_chg': DimuDstar_acc['Dstar']['wrg_chg'].value,}, with_name='PtEtaPhiMCandidate')  
        

        DimuDstar = ak.unflatten(DimuDstar, DimuDstar_acc['nDimuDstar'].value)

        # Trigger cut
        hlt = False
        hlt_filter_2016 = ['HLT_Dimuon20_Jpsi']        
        hlt_filter = hlt_filter_2016

        HLT_acc = HLT_2016_acc

        if hlt:
            print(f"You are running with the trigger(s): {hlt_filter}")
            
            trigger_cut = HLT_acc[hlt_filter[0]].value
            for i in range(0, len(hlt_filter)):
                trigger_cut |= HLT_acc[hlt_filter[i]].value

            # Muon lead collection
            Muon_lead = Muon_lead[trigger_cut]
              
            muon_lead_pt = ak.flatten(Muon_lead.pt)
            muon_lead_eta = ak.flatten(Muon_lead.eta)
            muon_lead_phi = ak.flatten(Muon_lead.phi)
            
            # Muon trail collection
            Muon_trail = Muon_trail[trigger_cut]

            muon_trail_pt = ak.flatten(Muon_trail.pt)
            muon_trail_eta = ak.flatten(Muon_trail.eta)
            muon_trail_phi = ak.flatten(Muon_trail.phi)

            # Jpsi collection
            Dimu = Dimu[trigger_cut]
            Dimu = Dimu[Dimu.is_jpsi]

            jpsi_mass = ak.flatten(Dimu.mass)
            jpsi_pt = ak.flatten(Dimu.pt)
            jpsi_eta = ak.flatten(Dimu.eta)
            jpsi_phi = ak.flatten(Dimu.phi)
            jpsi_rap = ak.flatten(Dimu.rap)

            jpsi_dl = ak.flatten(Dimu.dl)
            jpsi_dlSig = ak.flatten(Dimu.dlSig)
            jpsi_chi2 = ak.flatten(Dimu.chi2)
            jpsi_cosphi =ak.flatten(Dimu.cosphi)
           
            # Dstar collection
            Dstar = Dstar[trigger_cut]

            dstar_right_charge = Dstar[~Dstar.wrg_chg]
            dstar_wrong_charge = Dstar[Dstar.wrg_chg]

            dstar_right_charge_pt = ak.flatten(dstar_right_charge.pt)
            dstar_right_charge_eta = ak.flatten(dstar_right_charge.eta)
            dstar_right_charge_phi = ak.flatten(dstar_right_charge.phi)

            dstar_wrong_charge_pt = ak.flatten(dstar_wrong_charge.pt)
            dstar_wrong_charge_eta = ak.flatten(dstar_wrong_charge.eta)
            dstar_wrong_charge_phi = ak.flatten(dstar_wrong_charge.phi)
            
            dstar_right_charge_rap = ak.flatten(dstar_right_charge.rap)
            dstar_wrong_charge_rap = ak.flatten(dstar_wrong_charge.rap)

            dstar_right_charge_deltam = ak.flatten(dstar_right_charge.deltam)
            dstar_wrong_charge_deltam = ak.flatten(dstar_wrong_charge.deltam)

            dstar_right_charge_deltamr = ak.flatten(dstar_right_charge.deltamr)
            dstar_wrong_charge_deltamr = ak.flatten(dstar_wrong_charge.deltamr)

            # Vertex alignment 
            dstar_right_charge_cosphi = ak.flatten(dstar_right_charge.D0cosphi)
            dstar_wrong_charge_cosphi = ak.flatten(dstar_wrong_charge.D0cosphi)

            # Decay length significance
            dstar_right_charge_D0dlsig = ak.flatten(dstar_right_charge.D0dlSig)
            dstar_wrong_charge_D0dlsig = ak.flatten(dstar_wrong_charge.D0dlSig)

            dstar_right_charge_D0pt = ak.flatten(dstar_right_charge.D0pt)
            dstar_wrong_charge_D0pt = ak.flatten(dstar_wrong_charge.D0pt)

            ## DimuDstar collection

            # Trigger cut
            DimuDstar = DimuDstar[trigger_cut]

            # Takes Jpsi range
            DimuDstar = DimuDstar[DimuDstar.is_jpsi]

            # Takes only right charge D*
            DimuDstar = DimuDstar[~DimuDstar.wrg_chg]
            # Takes only wrong charge D*
            DimuDstar_wrg = DimuDstar[DimuDstar.wrg_chg]

            # Associated jpsi
            jpsi_asso_mass = ak.flatten(DimuDstar.jpsi_mass)
            jpsi_asso_pt = ak.flatten(DimuDstar.jpsi_pt)
            jpsi_asso_eta = ak.flatten(DimuDstar.jpsi_eta)
            jpsi_asso_phi = ak.flatten(DimuDstar.jpsi_phi)
            jpsi_asso_rap = ak.flatten(DimuDstar.jpsi_rap)

            ## Associated dstar

            # Cuts
            dstar_rgt_chg_cuts = is_jpsi & ~wrg_chg
            dstar_wrg_chg_cuts = is_jpsi & wrg_chg

            dimu_dstar_right_charge = DimuDstar[~DimuDstar.wrg_chg]
            dimu_dstar_wrong_charge = DimuDstar_wrg


            dstar_asso_right_charge_deltamr = ak.flatten(dimu_dstar_right_charge.dstar_deltamr)
            dstar_asso_wrong_charge_deltamr = ak.flatten(dimu_dstar_wrong_charge.dstar_deltamr)
            
            dstar_asso_right_charge_deltam = ak.flatten(dimu_dstar_right_charge.dstar_deltam)
            dstar_asso_wrong_charge_deltam = ak.flatten(dimu_dstar_wrong_charge.dstar_deltam)

            dstar_asso_right_charge_pt = ak.flatten(dimu_dstar_right_charge.dstar_pt)
            dstar_asso_right_charge_eta = ak.flatten(dimu_dstar_right_charge.dstar_eta)
            dstar_asso_right_charge_phi = ak.flatten(dimu_dstar_right_charge.dstar_phi)

            dstar_asso_wrong_charge_pt = ak.flatten(dimu_dstar_wrong_charge.dstar_pt)
            dstar_asso_wrong_charge_eta = ak.flatten(dimu_dstar_wrong_charge.dstar_eta)
            dstar_asso_wrong_charge_phi = ak.flatten(dimu_dstar_wrong_charge.dstar_phi)

            dstar_asso_right_charge_rap = ak.flatten(dimu_dstar_right_charge.dstar_rap)
            dstar_asso_wrong_charge_rap = ak.flatten(dimu_dstar_wrong_charge.dstar_rap)

            # Associated object (DimuDstar)
            dimuon_dstar_deltarap = ak.flatten(DimuDstar.dimu_dstar_deltarap)
            dimuon_dstar_associationIdx = ak.flatten(DimuDstar.dimu_dstar_associationIdx)
            dimuon_dstar_vtx_prob = ak.flatten(DimuDstar.dimu_dstar_vtx_prob)
            dimuon_dstar_vtx_chi2 = ak.flatten(DimuDstar.dimu_dstar_vtx_chi2)
            dimuon_dstar_mass = ak.flatten(DimuDstar.dimu_dstar_mass)
            
        if not hlt:
            #print("You are not running with trigger")
            trigger_cut = np.ones(len(Dimu), dtype=bool)
            
            # Muon lead collection
            muon_lead_pt = Muon_lead_acc['pt'].value
            muon_lead_eta = Muon_lead_acc['eta'].value
            muon_lead_phi = Muon_lead_acc['phi'].value
            
            # Muon trail collection
            muon_trail_pt = Muon_trail_acc['pt'].value
            muon_trail_eta = Muon_trail_acc['eta'].value
            muon_trail_phi = Muon_trail_acc['phi'].value
           
            # Jpsi collection
            jpsi_mass = Dimu_acc['mass'].value[Dimu_acc['is_jpsi'].value]
            jpsi_pt = Dimu_acc['pt'].value[Dimu_acc['is_jpsi'].value]
            jpsi_eta = Dimu_acc['eta'].value[Dimu_acc['is_jpsi'].value]
            jpsi_phi = Dimu_acc['phi'].value[Dimu_acc['is_jpsi'].value]
            jpsi_rap = Dimu_acc['rap'].value[Dimu_acc['is_jpsi'].value]
            jpsi_dl = Dimu_acc['dl'].value[Dimu_acc['is_jpsi'].value]
            jpsi_dlSig = Dimu_acc['dlSig'].value[Dimu_acc['is_jpsi'].value]
            jpsi_chi2 = Dimu_acc['chi2'].value[Dimu_acc['is_jpsi'].value]
            jpsi_cosphi = Dimu_acc['cosphi'].value[Dimu_acc['is_jpsi'].value]
           
            # Dstar collection
            dstar_right_charge_pt = Dstar_acc['pt'].value[~Dstar_acc['wrg_chg'].value]
            dstar_right_charge_eta = Dstar_acc['eta'].value[~Dstar_acc['wrg_chg'].value]
            dstar_right_charge_phi = Dstar_acc['phi'].value[~Dstar_acc['wrg_chg'].value]
            
            dstar_wrong_charge_pt = Dstar_acc['pt'].value[Dstar_acc['wrg_chg'].value]
            dstar_wrong_charge_eta = Dstar_acc['eta'].value[Dstar_acc['wrg_chg'].value]
            dstar_wrong_charge_phi = Dstar_acc['phi'].value[Dstar_acc['wrg_chg'].value]
            
            dstar_right_charge_rap = Dstar_acc['rap'].value[~Dstar_acc['wrg_chg'].value]
            dstar_wrong_charge_rap = Dstar_acc['rap'].value[Dstar_acc['wrg_chg'].value]

            dstar_right_charge_deltam = Dstar_acc['deltamr'].value[~Dstar_acc['wrg_chg'].value]
            dstar_wrong_charge_deltam = Dstar_acc['deltamr'].value[Dstar_acc['wrg_chg'].value]

            dstar_right_charge_deltamr = Dstar_acc['deltam'].value[~Dstar_acc['wrg_chg'].value]
            dstar_wrong_charge_deltamr = Dstar_acc['deltam'].value[Dstar_acc['wrg_chg'].value]

            # Vertex alignment 
            dstar_right_charge_cosphi = Dstar_D0_acc['D0cosphi'].value[~Dstar_acc['wrg_chg'].value] 
            dstar_wrong_charge_cosphi = Dstar_D0_acc['D0cosphi'].value[Dstar_acc['wrg_chg'].value]

            # Decay length significance
            dstar_right_charge_D0dlsig = Dstar_D0_acc['D0dlSig'].value[~Dstar_acc['wrg_chg'].value] 
            dstar_wrong_charge_D0dlsig = Dstar_D0_acc['D0dlSig'].value[Dstar_acc['wrg_chg'].value]

            # pT
            dstar_right_charge_D0pt = Dstar_D0_acc['D0pt'].value[~Dstar_acc['wrg_chg'].value] 
            dstar_wrong_charge_D0pt = Dstar_D0_acc['D0pt'].value[Dstar_acc['wrg_chg'].value]

            ## DimuonDstar

            # Filters for jpsi and dstar
            is_jpsi = DimuDstar_acc['Dimu']['is_jpsi'].value
            wrg_chg = DimuDstar_acc['Dstar']['wrg_chg'].value

            # Associated jpsi
            jpsi_asso_mass = DimuDstar_acc['Dimu']['mass'].value[is_jpsi & ~wrg_chg]
            jpsi_asso_pt = DimuDstar_acc['Dimu']['pt'].value[is_jpsi & ~wrg_chg]
            jpsi_asso_eta = DimuDstar_acc['Dimu']['eta'].value[is_jpsi & ~wrg_chg]
            jpsi_asso_phi = DimuDstar_acc['Dimu']['phi'].value[is_jpsi & ~wrg_chg]
            jpsi_asso_rap = DimuDstar_acc['Dimu']['rap'].value[is_jpsi & ~wrg_chg]

            # Associated dstar
            dstar_rgt_chg_cuts = is_jpsi & ~wrg_chg
            dstar_wrg_chg_cuts = is_jpsi & wrg_chg

            dstar_asso_right_charge_deltamr = DimuDstar_acc['Dstar']['deltamr'].value[dstar_rgt_chg_cuts]
            dstar_asso_wrong_charge_deltamr = DimuDstar_acc['Dstar']['deltamr'].value[dstar_wrg_chg_cuts]
            
            dstar_asso_right_charge_deltam = DimuDstar_acc['Dstar']['deltam'].value[dstar_rgt_chg_cuts]
            dstar_asso_wrong_charge_deltam = DimuDstar_acc['Dstar']['deltam'].value[dstar_wrg_chg_cuts]

            dstar_asso_right_charge_pt = DimuDstar_acc['Dstar']['pt'].value[dstar_rgt_chg_cuts]
            dstar_asso_right_charge_eta = DimuDstar_acc['Dstar']['eta'].value[dstar_rgt_chg_cuts]
            dstar_asso_right_charge_phi = DimuDstar_acc['Dstar']['phi'].value[dstar_rgt_chg_cuts]

            dstar_asso_wrong_charge_pt = DimuDstar_acc['Dstar']['pt'].value[dstar_wrg_chg_cuts]
            dstar_asso_wrong_charge_eta = DimuDstar_acc['Dstar']['eta'].value[dstar_wrg_chg_cuts]
            dstar_asso_wrong_charge_phi = DimuDstar_acc['Dstar']['phi'].value[dstar_wrg_chg_cuts]

            dstar_asso_right_charge_rap = DimuDstar_acc['Dstar']['rap'].value[dstar_rgt_chg_cuts]
            dstar_asso_wrong_charge_rap = DimuDstar_acc['Dstar']['rap'].value[dstar_wrg_chg_cuts]

            # Associated object
            dimuon_dstar_deltarap = DimuDstar_acc['deltarap'].value[is_jpsi & ~wrg_chg]
            dimuon_dstar_associationIdx = DimuDstar_acc['Dstar']['associationIdx'].value[is_jpsi & ~wrg_chg]
            dimuon_dstar_vtx_prob = DimuDstar_acc['Dstar']['associationProb'].value[is_jpsi & ~wrg_chg]
            dimuon_dstar_vtx_chi2 = DimuDstar_acc['Dstar']['associationchi2'].value[is_jpsi & ~wrg_chg]
            dimuon_dstar_mass = DimuDstar_p4.mass[is_jpsi & ~wrg_chg]

        ##Muon
        output['Muon_lead_p'].fill(pt=muon_lead_pt,
                                   eta=muon_lead_eta,
                                   phi=muon_lead_phi,)
        output['Muon_trail_p'].fill(pt=muon_trail_pt,
                                   eta=muon_trail_eta,
                                   phi=muon_trail_phi,)
     
        # Jpsi
        output['Jpsi_mass'].fill(mass=jpsi_mass)
        output['Jpsi_p'].fill(pt=jpsi_pt,
                                 eta=jpsi_eta,
                                 phi=jpsi_phi)
        output['Jpsi_rap'].fill(rap=jpsi_rap)
        output['Jpsi_dl'].fill(dl=jpsi_dl)
        output['Jpsi_dlSig'].fill(dlSig=jpsi_dlSig)
        output['Jpsi_chi2'].fill(chi2=jpsi_chi2)
        output['Jpsi_cosphi'].fill(cosphi=jpsi_cosphi)
        
        # Dstar
        output['Dstar_p'].fill(chg='right charge', 
                               pt=dstar_right_charge_pt,
                               eta=dstar_right_charge_eta,
                               phi=dstar_right_charge_phi)
        output['Dstar_p'].fill(chg='wrong charge', 
                               pt=dstar_wrong_charge_pt,
                               eta=dstar_wrong_charge_eta,
                               phi=dstar_wrong_charge_phi)
        output['Dstar_rap'].fill(chg='right charge', rap=dstar_right_charge_rap)
        output['Dstar_rap'].fill(chg='wrong charge', rap=dstar_wrong_charge_rap)
        output['Dstar_deltamr'].fill(chg='right charge', deltamr=dstar_right_charge_deltam)
        output['Dstar_deltamr'].fill(chg='wrong charge', deltamr=dstar_wrong_charge_deltam)
        output['Dstar_deltam'].fill(chg='right charge', deltam=dstar_right_charge_deltamr)
        output['Dstar_deltam'].fill(chg='wrong charge', deltam=dstar_wrong_charge_deltamr)
        output['Dstar_D0cosphi'].fill(chg='right charge', cosphi=dstar_right_charge_cosphi)
        output['Dstar_D0cosphi'].fill(chg='wrong charge', cosphi=dstar_wrong_charge_cosphi)
        output['Dstar_D0dlSig'].fill(chg='right charge', dlSig=dstar_right_charge_D0dlsig)
        output['Dstar_D0dlSig'].fill(chg='wrong charge', dlSig=dstar_wrong_charge_D0dlsig)
        output['Dstar_D0pt'].fill(chg='right charge', pt=dstar_right_charge_D0pt)
        output['Dstar_D0pt'].fill(chg='wrong charge', pt=dstar_wrong_charge_D0pt)
        # Dstar trks
        output['Dstar_K_p'].fill(pt=Dstar_trk_acc['Kpt'].value[~Dstar_acc['wrg_chg'].value],
                                 eta=Dstar_trk_acc['Keta'].value[~Dstar_acc['wrg_chg'].value],
                                 phi=Dstar_trk_acc['Kphi'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_K_chindof'].fill(chindof=Dstar_trk_acc['Kchindof'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_K_nValid'].fill(nValid=Dstar_trk_acc['KnValid'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_K_nPix'].fill(nPix=Dstar_trk_acc['KnPix'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_K_dxy'].fill(dxy=Dstar_trk_acc['Kdxy'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_K_dz'].fill(dz=Dstar_trk_acc['Kdz'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_K_pt_eta'].fill(pt=Dstar_trk_acc['Kpt'].value[~Dstar_acc['wrg_chg'].value],
                                      eta=Dstar_trk_acc['Keta'].value[~Dstar_acc['wrg_chg'].value])

        output['Dstar_pi_p'].fill(pt=Dstar_trk_acc['pipt'].value[~Dstar_acc['wrg_chg'].value],
                                  eta=Dstar_trk_acc['pieta'].value[~Dstar_acc['wrg_chg'].value],
                                  phi=Dstar_trk_acc['piphi'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_pi_chindof'].fill(chindof=Dstar_trk_acc['pichindof'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_pi_nValid'].fill(nValid=Dstar_trk_acc['pinValid'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_pi_nPix'].fill(nPix=Dstar_trk_acc['pinPix'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_pi_dxy'].fill(dxy=Dstar_trk_acc['pidxy'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_pi_dz'].fill(dz=Dstar_trk_acc['pidz'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_pi_pt_eta'].fill(pt=Dstar_trk_acc['pipt'].value[~Dstar_acc['wrg_chg'].value],
                                       eta=Dstar_trk_acc['pieta'].value[~Dstar_acc['wrg_chg'].value])

        output['Dstar_pis_p'].fill(pt=Dstar_trk_acc['pispt'].value[~Dstar_acc['wrg_chg'].value],
                                   eta=Dstar_trk_acc['piseta'].value[~Dstar_acc['wrg_chg'].value],
                                   phi=Dstar_trk_acc['pisphi'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_pis_chindof'].fill(chindof=Dstar_trk_acc['pischindof'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_pis_nValid'].fill(nValid=Dstar_trk_acc['pisnValid'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_pis_nPix'].fill(nPix=Dstar_trk_acc['pisnPix'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_pis_dxy'].fill(dxy=Dstar_trk_acc['pisdxy'].value[~Dstar_acc['wrg_chg'].value])
        output['Dstar_pis_dz'].fill(dz=Dstar_trk_acc['pisdz'].value[~Dstar_acc['wrg_chg'].value])

        ############# DimuDstar
     

        # JpsiDstar
        output['JpsiDstar']['Jpsi_mass'].fill(mass=jpsi_asso_mass)
        output['JpsiDstar']['Jpsi_p'].fill(pt=jpsi_asso_pt,
                                           eta=jpsi_asso_eta,
                                           phi=jpsi_asso_phi)
        output['JpsiDstar']['Jpsi_rap'].fill(rap=jpsi_asso_rap)

        output['JpsiDstar']['Dstar_deltamr'].fill(chg='right charge', deltamr=dstar_asso_right_charge_deltamr)
        output['JpsiDstar']['Dstar_deltamr'].fill(chg='wrong charge', deltamr=dstar_asso_wrong_charge_deltamr)
        output['JpsiDstar']['Dstar_deltam'].fill(chg='right charge', deltam=dstar_asso_right_charge_deltam)
        output['JpsiDstar']['Dstar_deltam'].fill(chg='wrong charge', deltam=dstar_asso_wrong_charge_deltam)
        output['JpsiDstar']['Dstar_p'].fill(chg='right charge',
                                            pt=dstar_asso_right_charge_pt,
                                            eta=dstar_asso_right_charge_eta,
                                            phi=dstar_asso_right_charge_phi)
        output['JpsiDstar']['Dstar_p'].fill(chg='wrong charge',
                                            pt=dstar_asso_wrong_charge_pt,
                                            eta=dstar_asso_wrong_charge_eta,
                                            phi=dstar_asso_wrong_charge_phi)
        output['JpsiDstar']['Dstar_rap'].fill(chg='right charge', rap=dstar_asso_right_charge_rap)
        output['JpsiDstar']['Dstar_rap'].fill(chg='wrong charge', rap=dstar_asso_wrong_charge_rap)

        output['JpsiDstar']['JpsiDstar_deltarap'].fill(deltarap=dimuon_dstar_deltarap)
        output['JpsiDstar']['JpsiDstar_associationIdx'].fill(associationIdx=dimuon_dstar_associationIdx)
        output['JpsiDstar']['JpsiDstar_vtx_prob'].fill(vtx_prob=dimuon_dstar_vtx_prob)
        output['JpsiDstar']['JpsiDstar_vtx_chi2'].fill(vtx_chi2=dimuon_dstar_vtx_chi2)
        output['JpsiDstar']['JpsiDstar_mass'].fill(mass=dimuon_dstar_mass)

       
        return output

    def postprocess(self, accumulator):
        return accumulator      
