# coding: utf-8

"""
Example inference model.
"""

from columnflow.inference import inference_model, ParameterType, ParameterTransformation


@inference_model
def main(self):
    # name of the (single) config this model is built for
    config_inst = self.config_insts[0]
    config_name = config_inst.name

    #
    # categories
    #

    self.add_category(
        "emu",
        config_data={
            config_name: self.category_config_spec(
                category="emu__real_1__real_2__D",
                variable="zcand_invm",
                data_datasets=["data_muoneg_A", "data_muoneg_B", "data_muoneg_C", "data_muoneg_D"],
            ),
        },
        #data_from_processes=["TT", "ST", "VV", "VVV", "W", "EWKW", "EWKZ", "Higgs", "DY_ll", "DY_2tau_minus", "DY_2tau_plus"],
        mc_stats=True,
    )

    #
    # processes
    #

    def add_proc(name, config_process, datasets, **kwargs):
        self.add_process(
            name,
            config_data={
                config_name: self.process_config_spec(
                    process=config_process,
                    mc_datasets=datasets,
                ),
            },
            **kwargs,
        )

    # backgrounds
    add_proc("TT", "tt", ["tt_sl_powheg", "tt_dl_powheg", "tt_fh_powheg",])
    add_proc("ST", "st", [
        "st_tchannel_t_powheg", "st_tchannel_tbar_powheg",
        "st_tw_t_powheg", "st_tw_tbar_powheg",
    ])
    add_proc("VV", "diboson", [
        "ww_dl_powheg", "ww_sl_powheg",
        "wz_wlnu_zll_powheg", "wz_wqq_zll_powheg", "wz_wlnu_zqq_powheg",
        "zz_zqq_zll_powheg", "zz_zll_zll_powheg",
    ])
    add_proc("VVV", "triboson", [
        "www_amcatnlo", "wwz_amcatnlo", "wzz_amcatnlo", "zzz_amcatnlo",
    ])
    add_proc("EWKW", "ewk_w", ["ewk_wp_lnu_madgraph", "ewk_wm_lnu_madgraph",])
    add_proc("EWKZ", "ewk_z_ll_m50toinf", ["ewk_z_ll_madgraph",])
    add_proc("Higgs", "higgs", [
        "h_ggf_htt_powheg", "zh_htt_powheg", "wph_htt_powheg", "wmh_htt_powheg",
    ])
    add_proc("W", "w_lnu", ["wj_incl_madgraph"])
    add_proc("DY_ll", "dy_2e_or_2mu_m50toinf_nj", ["dy_lep_m50_amcatnlo"])
    add_proc("qcd", "qcd", [], is_dynamic=True)

    # signals
    add_proc("DY_2tau_minus", "dy_2tau_m50toinf_nj_LHEspin_minus",
             ["dy_2tau_m50_amcatnlo"], is_signal=True)
    add_proc("DY_2tau_plus", "dy_2tau_m50toinf_nj_LHEspin_plus",
             ["dy_2tau_m50_amcatnlo"], is_signal=True)

    #
    # parameters
    #

    # lumi
    self.add_parameter(
        "lumi_13TeV",
        type=ParameterType.rate_gauss,
        effect=1.0084,
        process="!qcd",
    )

    """ lumi = config_inst.x.luminosity
    for unc_name in lumi.uncertainties:
        self.add_parameter(
            unc_name,
            type=ParameterType.rate_gauss,
            effect=lumi.get(names=unc_name, direction=("down", "up"), factor=True),
            #transformations=[ParameterTransformation.symmetrize],
        ) """

    # trigger
    """ self.add_parameter(
        "trigger_emu",
        type=ParameterType.shape,
        config_data={config_name: self.parameter_config_spec(shift_source="trigger_emu")}
    ) """

    # trigger
    self.add_parameter(
        "trigger_emu",
        type=ParameterType.rate_gauss,
        effect=1.02,
        process="!qcd",
    )

    # cross sections / normalisations (process-specific)
    self.add_parameter(
        "xs_DY",
        type=ParameterType.rate_gauss,
        effect=1.056,
        process=["DY_ll", "DY_2tau_minus", "DY_2tau_plus"],
    )

    self.add_parameter(
        "xs_TT",
        type=ParameterType.rate_gauss,
        effect=1.042,
        process="TT",
    )

    self.add_parameter(
        "xs_VV_ST",
        type=ParameterType.rate_gauss,
        effect=1.05,
        process=["VV", "ST"],
    )

    self.add_parameter(
        "xs_EWKZ",
        type=ParameterType.rate_gauss,
        effect=1.04,
        process="EWKZ",
    )

    self.add_parameter(
        "norm_W",
        type=ParameterType.rate_gauss,
        effect=1.20,
        process="W",
    )

    self.add_parameter(
        "xs_W",
        type=ParameterType.rate_gauss,
        effect=1.04,
        process="W",
    )

    # qcd normalisation (the only parameter acting on qcd)
    self.add_parameter(
        "norm_qcd",
        type=ParameterType.rate_gauss,
        effect=1.10,
        process="qcd",
    )

    # PU
    self.add_parameter(
        "PU",
        type=ParameterType.shape,
        config_data={config_name: self.parameter_config_spec(shift_source="minbias_xs")},
        process="!qcd",
    )

    # Jet resolution
    self.add_parameter(
        "JER",
        type=ParameterType.shape,
        config_data={config_name: self.parameter_config_spec(shift_source="jer")},
        process="!qcd",
    )

    # Electron corrections
    self.add_parameter(
        "ele_id",
        type=ParameterType.shape,
        config_data={config_name: self.parameter_config_spec(shift_source="e_id")},
        process="!qcd",
    )

    self.add_parameter(
        "ele_reco",
        type=ParameterType.shape,
        config_data={config_name: self.parameter_config_spec(shift_source="e_reco")},
        process="!qcd",
    )

    """ self.add_parameter(
        "ele_scale",
        type=ParameterType.shape,
        config_data={config_name: self.parameter_config_spec(shift_source="eec")}
    ) """

    """ self.add_parameter(
        "ele_smear",
        type=ParameterType.shape,
        config_data={config_name: self.parameter_config_spec(shift_source="eer")}
    ) """

    # Muon corrections
    self.add_parameter(
        "muon_id",
        type=ParameterType.shape,
        config_data={config_name: self.parameter_config_spec(shift_source="mu_id")},
        process="!qcd",
    )

    self.add_parameter(
        "muon_iso",
        type=ParameterType.shape,
        config_data={config_name: self.parameter_config_spec(shift_source="mu_iso")},
        process="!qcd",
    )

    """ self.add_parameter(
        "muon_scale",
        type=ParameterType.shape,
        config_data={config_name: self.parameter_config_spec(shift_source="mec")}
    ) """

    """ self.add_parameter(
        "muon_res",
        type=ParameterType.shape,
        config_data={config_name: self.parameter_config_spec(shift_source="mer")}
    ) """

    # Theory
    """ self.add_parameter(
        "QCD_scale_fact",
        type=ParameterType.shape,
        config_data={config_name: self.parameter_config_spec(shift_source="murmuf")}
    )

    self.add_parameter(
        "ISR",
        type=ParameterType.shape,
        config_data={config_name: self.parameter_config_spec(shift_source="isr")}
    )

    self.add_parameter(
        "FSR",
        type=ParameterType.shape,
        config_data={config_name: self.parameter_config_spec(shift_source="fsr")}
    ) """