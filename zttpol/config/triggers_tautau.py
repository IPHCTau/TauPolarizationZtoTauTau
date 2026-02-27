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
        # ===>>> tauh-tauh
        Trigger(
            name="HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg",
            id=15151,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=36.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=36.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc or dataset_inst.is_mc) if postfix == "preVFP" else (lambda dataset_inst: dataset_inst.x.era < "H"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg",
            id=15152,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=36.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=36.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_mc) if postfix == "preVFP" else (lambda dataset_inst: dataset_inst.is_mc or dataset_inst.x.era == "H"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
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
        # ===>>> tauh-tauh
        Trigger(
            name="HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg",
            id=15151,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
            id=15152,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=42.0,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            #applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
            id=15153,
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=37.0,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=37.0,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            #applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
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
        # ===>>> tauh-tauh 
        Trigger(
            name="HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
            id=15151,
            run_range=(None,317509), # cover up to 317509
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg",
            id=15152,
            run_range=(None,317509), # cover up to 317509
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=45.0,
                    max_abseta=None,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        Trigger(
            name="HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
            id=15153,
            run_range=(None,317509), # cover up to 317509
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        #Trigger(
        #    name="HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
        #    id=15152,
        #    run_range=(None,317510), # cover up to 317509
        #    legs=[
        #        TriggerLeg(
        #            pdg_id=15,
        #            min_pt=45.0,
        #            max_abseta=None,
        #            # filter names:
        #            # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
        #            trigger_bits=64,
        #        ),
        #        TriggerLeg(
        #            pdg_id=15,
        #            min_pt=45.0,
        #            max_abseta=None,
        #            # filter names:
        #            # hltDoublePFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
        #            trigger_bits=64,
        #        ),
        #    ],
        #    applies_to_dataset=(lambda dataset_inst: dataset_inst.is_data and dataset_inst.x.era <= "B"),
        #    tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        #),
        Trigger(
            name="HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg",
            id=15154,
            run_range=(317508,None), # cover 317509
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    max_abseta=2.1,
                    # filter names:
                    # hltDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg
                    trigger_bits=64,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
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
        ## ===>>> tauh-tauh
        # https://cmshltinfo.app.cern.ch/path/HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1_v
        # https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L141C73-L141C87
        # https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv14/2024Prompt/doc_EGamma1_Run2024D-PromptReco-v1.html#TrigObj
        Trigger(
            name="HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1",
            id=15151,
            #run_range=[355862,362760],
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    min_pt_online=35.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsSelectedPFTausMediumDitauWPDeepTau
                    # hltHpsDoublePFTau35MediumDitauWPDeepTauL1HLTMatched
                    # 3 => DeepTau no spec WP, 11 => di-tau
                    trigger_bits=2**3 + 2**11,
               ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    min_pt_online=35.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsDoublePFTau35MediumDitauWPDeepTauDz02 (Deeptau + HPS)
                    # 3 => DeepTau no spec WP, 11 => di-tau 
                    trigger_bits=2**3 + 2**11,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        # ==> diTau + Jet
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger#Trigger_Table_for_2022
        # https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L148
        # https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv14/2024Prompt/doc_EGamma1_Run2024D-PromptReco-v1.html#TrigObj
        # same as https://gitlab.cern.ch/dwinterb/HiggsDNA/-/blob/NewProduction_Run3-2022/higgs_dna/metaconditions/Era2022_ditau.json?ref_type=heads#L146
        Trigger(
            name="HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60",
            id=15152,
            #run_range=[355862,362760],
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60
                    # 14 => di-tau + PFJet
                    trigger_bits=2**3 + 2**14,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60
                    # 14 => di-tau + PFJet
                    trigger_bits=2**3 + 2**14,
                ),
                TriggerLeg(
                    pdg_id=1,
                    min_pt=60.0,
                    min_pt_online=60.0,
                    max_abseta=4.7,
                    max_abseta_online=4.9,
                    # filter names:
                    # hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60
                    # 14 => di-tau + PFJet
                    trigger_bits=2**17,
                ),                
            ],
            tags={"cross_trigger", "cross_tau_tau_jet", "channel_tau_tau"},
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
        # ===>>> tauh-tauh
        Trigger(
            name="HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1",
            id=15151,
            #run_range=[355862,362760],
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    min_pt_online=35.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsSelectedPFTausMediumDitauWPDeepTau
                    # hltHpsDoublePFTau35MediumDitauWPDeepTauL1HLTMatched
                    # 3 => DeepTau no spec WP, 11 => di-tau
                    trigger_bits=2**3 + 2**11,
               ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=40.0,
                    min_pt_online=35.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsDoublePFTau35MediumDitauWPDeepTauDz02 (Deeptau + HPS)
                    # 3 => DeepTau no spec WP, 11 => di-tau 
                    trigger_bits=2**3 + 2**11,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau", "channel_tau_tau"},
        ),
        # ==> diTau + Jet
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger#Trigger_Table_for_2022
        # https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L148
        # https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv14/2024Prompt/doc_EGamma1_Run2024D-PromptReco-v1.html#TrigObj
        # same as https://gitlab.cern.ch/dwinterb/HiggsDNA/-/blob/NewProduction_Run3-2022/higgs_dna/metaconditions/Era2022_ditau.json?ref_type=heads#L146
        Trigger(
            name="HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60",
            id=15152,
            #run_range=[355862,362760],
            legs=[
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60
                    # 14 => di-tau + PFJet
                    trigger_bits=2**3 + 2**14,
                ),
                TriggerLeg(
                    pdg_id=15,
                    min_pt=35.0,
                    min_pt_online=30.0,
                    max_abseta=2.1,
                    max_abseta_online=2.1,
                    # filter names:
                    # hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60
                    # 14 => di-tau + PFJet
                    trigger_bits=2**3 + 2**14,
                ),
                TriggerLeg(
                    pdg_id=1,
                    min_pt=60.0,
                    min_pt_online=60.0,
                    max_abseta=4.7,
                    max_abseta_online=4.9,
                    # filter names:
                    # hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60
                    # 14 => di-tau + PFJet
                    trigger_bits=2**17,
                ),
            ],
            tags={"cross_trigger", "cross_tau_tau_jet", "channel_tau_tau"},
        ),
    ])
