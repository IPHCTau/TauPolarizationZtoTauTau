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
        # ===>>> single electron
        Trigger(
            name="HLT_Ele25_eta2p1_WPTight_Gsf",
            id=111,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=26.0,
                    # filter names:
                    # hltEle25erWPTightGsfTrackIsoFilter
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
    ])    



# ----------------------------------------------- #
#                   Run2 2017 UL                  #
# ----------------------------------------------- #
def add_triggers_2017(config: od.Config) -> None:
    """
    Adds all triggers to a *config*. For the conversion from filter names to trigger bits, see
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
    """
    config.x.triggers = od.UniqueObjectIndex(Trigger, [
        # ===>>> single electron
        Trigger(
            name="HLT_Ele27_WPTight_Gsf",
            id=111,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=28.0,
                    # filter names:
                    # hltEle32L1DoubleEGWPTightGsfTrackIsoFilter
                    # hltEGL1SingleEGOrFilter
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        Trigger(
            name="HLT_Ele32_WPTight_Gsf",
            id=112,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=28.0,
                    # filter names:
                    # hltEle32WPTightGsfTrackIsoFilter
                    trigger_bits=2 + 1024,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era >= "D"),
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        # ===>>> e-tauh
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1",
            id=11151,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=25.0,
                    # filter names:
                    # hltEle24erWPTightGsfTrackIsoFilterForTau
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=32.0,
                    # filter names:
                    # hltSelectedPFTau30LooseChargedIsolationL1HLTMatched
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=1024 + 256,
                ),
            ],
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
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
        # ===>>> single electron
        Trigger(
            name="HLT_Ele32_WPTight_Gsf",
            id=111000,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=33.0,
                    max_abseta=None,
                    # filter names:
                    # hltEle32WPTightGsfTrackIsoFilter
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        Trigger(
            name="HLT_Ele35_WPTight_Gsf",
            id=112000,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=33.0,
                    max_abseta=None,
                    # filter names:
                    # hltEle35noerWPTightGsfTrackIsoFilter
                    trigger_bits=2,
                ),
            ],
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        # ===>>> e-tauh
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1",
            id=11151,
            run_range=(None,317509),
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=26.0, #25.0,
                    max_abseta=2.1, #None,
                    # filter names:
                    # hltEle24erWPTightGsfTrackIsoFilterForTau
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=64, #2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    max_abseta=2.1, #None,
                    # filter names:
                    # hltSelectedPFTau30LooseChargedIsolationL1HLTMatched
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=256, #1024 + 256,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
        ),
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1",
            id=11152,
            run_range=(317508,None),
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=26.0, #25.0,
                    max_abseta=2.1, #None,
                    # filter names:
                    # hltEle24erWPTightGsfTrackIsoFilterForTau
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=2 + 64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    max_abseta=2.1, #None,
                    # filter names:
                    # hltSelectedPFTau30LooseChargedIsolationL1HLTMatched
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=1024 + 256,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era >= "B"),
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
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
        ## ===>>> single electron
        Trigger(
            name="HLT_Ele30_WPTight_Gsf",
            id=111000,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=31.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.3,
                    # filter names:
                    # hltEle30WPTightGsfTrackIsoFilter
                    trigger_bits=2**1,  # 1e (WPTight) (bit 1)
                ),
            ],
            tags=["single_trigger", "single_e", "channel_e_tau"],
        ),
        ## ===>>> e-tauh
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1",
            id=11151,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=25.0, # 26 -> Imperial
                    min_pt_online=24.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltOverlapFilterIsoEle24IsoTau30WPTightGsfCaloJet5
                    # hltHpsOverlapFilterIsoEle24WPTightGsfLooseETauWPDeepTauPFTau30
                    trigger_bits=2**7,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltOverlapFilterIsoEle24IsoTau30WPTightGsfCaloJet5
                    # hltHpsOverlapFilterIsoEle24WPTightGsfLooseETauWPDeepTauPFTau30
                    trigger_bits=2**12,
                ),
            ],
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
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
        # ===>>> single electron
        Trigger(
            name="HLT_Ele30_WPTight_Gsf",
            id=111000,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=31.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.3,
                    # filter names:
                    # 
                    trigger_bits=2**1,
                ),
            ],
            #applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era >= "D"),
            tags={"single_trigger", "single_e", "channel_e_tau"},
        ),
        # ===>>> e-tauh
        Trigger(
            name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1",
            id=11151,
            legs=[
                TriggerLeg(
                    pdg_id=11,
                    min_pt=25.0,
                    min_pt_online=24.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltEle24erWPTightGsfTrackIsoFilterForTau
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=2**7,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltSelectedPFTau30LooseChargedIsolationL1HLTMatched
                    # hltOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30
                    trigger_bits=2**12,
                ),
            ],
            tags={"cross_trigger", "cross_e_tau", "channel_e_tau"},
        ),
    ])
