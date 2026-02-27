# coding: utf-8

"""
Definition of triggers
"""

import order as od

from zttpol.config.trigger_util import Trigger, TriggerLeg



# ----------------------------------------------- #
#                   Run2 2016 UL                  #
# ----------------------------------------------- #
def add_triggers_2016(config: od.Config, postfix: str) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger, [
        # ===>>> single muon
        Trigger(
            name="HLT_IsoMu22",
            id=131,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=23.0,
                    # filter names:
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu22_eta2p1",
            id=132,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=23.0,
                    # filter names:
                    # hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        # ===>>> mu-tauh
        Trigger(
            name="HLT_IsoMu19_eta2p1_LooseIsoPFTau20",
            id=13151,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=20.0,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=21.0,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=1024 + 512,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1",
            id=13152,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=20.0,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=21.0,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=1024 + 512,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
    ])    



# ----------------------------------------------- #
#                   Run2 2017 UL                  #
# ----------------------------------------------- #
def add_triggers_UL2017(config: od.Config) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger, [
        # ===>>> single muon
        Trigger(
            name="HLT_IsoMu24",
            id=131,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    # filter names:
                    # hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu27",
            id=132,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    # filter names:
                    # hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        # ===>>> mu-tauh
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1",
            id=13151,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=21.0,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=1024 + 512,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
    ])    



# ----------------------------------------------- #
#                   Run2 2018 UL                  #
# ----------------------------------------------- #
def add_triggers_2018(config: od.Config) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger, [
        # ===>>> single muon
        Trigger(
            name="HLT_IsoMu24",
            id=131000,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    max_abseta=None,
                    # filter names:
                    # hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu27",
            id=132000,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=25.0,
                    max_abseta=None,
                    # filter names:
                    # hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        # ===>>> mu-tauh
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1",
            id=13151,
            run_range=(None,317509), #315974),
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=22.0, # 21.0
                    max_abseta=2.1, # None
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=64, #2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15, 
                    min_pt=32.0,
                    max_abseta=2.1,  #None,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=64, #1024 + 512,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1",
            id=13153,
            run_range=(317508,None),
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=22.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=64, #2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=512, #1024 + 512,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era >= "B"),
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),        
    ])



# ----------------------------------------------- #
#                   Run3 2022                     #
# ----------------------------------------------- #
def add_triggers_2022(config: od.Config, postfix: str) -> None:
    """
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger#Tau_Triggers_in_NanoAOD
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    ** Tau Trigger: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger#Trigger_Table_for_2022
    ** Electron Trigger: https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIIISummary
    ** Muon Trigger: https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2022
    ** HiggsDNA : https://gitlab.cern.ch/dwinterb/HiggsDNA/-/blob/NewProduction_Run3-2022/higgs_dna/metaconditions/Era2022_ditau.json?ref_type=heads
    PreEE:
      /afs/cern.ch/work/g/gsaha/public/IPHC/Work/ColumnFlowAnalyses/CPinHToTauTau/yamls/HLTlog_2022PreEE.log
    PostEE:
      /afs/cern.ch/work/g/gsaha/public/IPHC/Work/ColumnFlowAnalyses/CPinHToTauTau/yamls/HLTlog_2022PostEE.log
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger,[
        ## ===>>> single muon
        # https://cmshltinfo.app.cern.ch/path/HLT_IsoMu24_v#state=53332ee4-249a-4fd3-8152-5ddbcec178c6&session_state=4d32aa9a-e3a3-4947-8c1a-6d96d0a49833&code=90d18d3f-5685-4f64-899f-02aa42d33b2c.4d32aa9a-e3a3-4947-8c1a-6d96d0a49833.1363e04b-e180-4d83-92b3-3aca653d1d8d
        Trigger(
            name="HLT_IsoMu24",
            id=131000,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=26.0,
                    min_pt_online=24.0,
                    max_abseta=2.4,
                    max_abseta_online=2.4,
                    # filter names:
                    # hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p08
                    trigger_bits=2**1 + 2**3,  # Iso (bit 1) + 1mu (bit 3)
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        ## ===>>> mu-tauh
        # https://cmshltinfo.app.cern.ch/path/HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1_v
        # mu : https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L118C66-L118C74
        # tau : https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L143
        # https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv14/2024Prompt/doc_EGamma1_Run2024D-PromptReco-v1.html#TrigObj
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1",
            id=13151,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=21.0,
                    min_pt_online=20.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # 
                    trigger_bits=2**6,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    min_pt_online=27.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsOverlapFilterIsoMu20LooseMuTauWPDeepTauPFTau27L1Seeded
                    trigger_bits=2**13,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
    ])


# ----------------------------------------------- #
#                   Run3 2023                     #
# ----------------------------------------------- #
def add_triggers_2023(config: od.Config, postfix: str) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    PreBPix:
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger,[
        # ===>>> single muon
        Trigger(
            name="HLT_IsoMu24",
            id=131000,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=26.0,
                    min_pt_online=24.0,
                    max_abseta=2.4,
                    max_abseta_online=2.4,
                    # filter
                    # hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered
                    # https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv14/2024Prompt/doc_TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_RunIII2024Summer24NanoAOD-140X_mcRun3_2024_realistic_v26-v2.html#TrigObj
                    # https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L115
                    trigger_bits=2**1 + 2**3, # Iso (bit 1) + 1mu (bit 3)
                ),
            ],
            tags={"single_trigger", "single_mu", "channel_mu_tau"},
        ),
        # ===>>> mu-tauh
        Trigger(
            name="HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1",
            id=13151,
            legs=[
                TriggerLeg(
                    pdg_id=13,
                    min_pt=21.0,
                    min_pt_online=20.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2**6,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    min_pt_online=27.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched or
                    # hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded
                    trigger_bits=2**13,
                ),
            ],
            tags={"cross_trigger", "cross_mu_tau", "channel_mu_tau"},
        ),
    ])
