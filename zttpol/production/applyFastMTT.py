# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

import os
import sys
from pathlib import Path

import law
import order as od
from typing import Optional
from columnflow.production import Producer, producer
from columnflow.production.util import attach_coffea_behavior

from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, remove_ak_column
from columnflow.columnar_util import optional_column as optional

from zttpol.production.FastMTT import FastMTT
from zttpol.wrappers import KinFit
from zttpol.wrappers import classicsvfit_cpp


np = maybe_import("numpy")
ak = maybe_import("awkward")
pd = maybe_import("pandas")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

logger = law.logger.get_logger(__name__)


def as_float64(array):
    """Convert a flat Awkward/NumPy column into contiguous float64."""
    arr = np.ascontiguousarray(
        ak.to_numpy(ak.fill_none(array, np.nan)).squeeze(axis=1),
        dtype=np.float64,
    )    
    return arr

def as_int32(array):
    """Convert a flat Awkward/NumPy column into contiguous float64."""
    arr = np.ascontiguousarray(
        ak.to_numpy(ak.fill_none(array, np.nan)).squeeze(axis=1),
        dtype=np.int32,
    )
    return arr
    

@producer(
    uses={
        # nano columns
        "zcand.*", "channel_id",
        "PuppiMET.{pt,phi,covXX,covXY,covYY}",
        "PVBS.*"
    },
)
def apply_fastMTT(
        self: Producer, 
        events: ak.Array,
        run_fmtt = False,
        run_fmtt_kinfit = False,
        run_classic_svfit = False,
        **kwargs
) -> ak.Array:

    events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
    zcand_ = ak.with_name(events.zcand, "PtEtaPhiMLorentzVector")
    hmass  = (zcand_[:,:1] + zcand_[:,1:2]).mass

    # dummy if fastMTT does not run
    zcand_pt_fastMTT = zcand_.pt
    zcand_eta_fastMTT = zcand_.eta
    zcand_phi_fastMTT = zcand_.phi
    zcand_mass_fastMTT = zcand_.mass
    mass_h = hmass * 1.2

    Z_mass = 91.1876 # GeV
    Z_width = {
        "constant": 20.0,
        "bw": 2.4952,
        "gaussian": 3.0,
        "exact": 2.0, # not required
    }

    etau_id   = self.config_inst.get_channel("etau").id
    mutau_id  = self.config_inst.get_channel("mutau").id
    tautau_id = self.config_inst.get_channel("tautau").id

    
    
    if run_fmtt:
        

        metpt  = ak.to_numpy(events.PuppiMET.pt[:,None])
        metphi = ak.to_numpy(events.PuppiMET.phi[:,None])
        metcovxx = ak.to_numpy(events.PuppiMET.covXX[:,None])
        metcovxy = ak.to_numpy(events.PuppiMET.covXY[:,None])
        metcovyy = ak.to_numpy(events.PuppiMET.covYY[:,None])
        
        pt1    = ak.to_numpy(events.zcand.pt[:,0:1])
        eta1   = ak.to_numpy(events.zcand.eta[:,0:1])
        phi1   = ak.to_numpy(events.zcand.phi[:,0:1])
        mass1  = ak.to_numpy(events.zcand.mass[:,0:1])
        
        pt2    = ak.to_numpy(events.zcand.pt[:,1:2])
        eta2   = ak.to_numpy(events.zcand.eta[:,1:2])
        phi2   = ak.to_numpy(events.zcand.phi[:,1:2])
        mass2  = ak.to_numpy(events.zcand.mass[:,1:2])

        dm1    = events.zcand.decayMode[:,0:1]
        dm1_dummy = ak.values_astype(-1 * ak.ones_like(dm1), np.int32)
        dm1    = ak.to_numpy(ak.where(events.channel_id < 4, dm1_dummy, dm1))
        dm2    = ak.to_numpy(events.zcand.decayMode[:,1:2])

        # decay_type:
        #1 - TauToHad
        #2 - TauToElec
        #3 - TauToMu
        type_1 = ak.ones_like(hmass)
        type_2 = 2 * type_1
        type_3 = 3 * type_1
        
        type1  = ak.to_numpy(ak.where(events.channel_id == etau_id,
                                      type_2,
                                      ak.where(events.channel_id == mutau_id,
                                               type_3,
                                               type_1)
                                      )
                             )
        type2  = ak.to_numpy(type_1)


        
        higgs_mass = ak.to_numpy(hmass)
        
        
        # prepare inputs properly
        metcov = np.concatenate((metcovxx,metcovxy,metcovxy,metcovyy), axis=1)
        metcov = np.reshape(metcov, (metcov.shape[0], 2, 2))
        measuredTauLeptons = np.concatenate((dm1,
                                             pt1,
                                             eta1,
                                             phi1,
                                             mass1,
                                             type1,
                                             dm2,
                                             pt2,
                                             eta2,
                                             phi2,
                                             mass2,
                                             type2), axis=1)
        measuredTauLeptons = np.reshape(measuredTauLeptons, (measuredTauLeptons.shape[0], 2, 6))
        metx   = (metpt * np.cos(metphi)).reshape(-1)
        mety   = (metpt * np.sin(metphi)).reshape(-1)
        
        
        # Launch FastMTT
        fMTT = FastMTT(enable_gauss=False,
                       enable_BW=False,
                       enable_window=True,
                       constrain_window=[70.0, 110.0],
                       massX=Z_mass,
                       widthX=Z_width['bw'])
        #You can choose to plot likelihood for one of the events. -1 means no plot.
        fMTT.WhichLikelihoodPlot = -1
        #You can also choose to calculate uncertainties by:
        fMTT.CalculateUncertainties = True
        
        fMTT.run(measuredTauLeptons, metx, mety, metcov)
        
        mass_h = fMTT.mass
        mass_h = ak.from_regular(ak.Array(mass_h.reshape(mass_h.shape[0],1)))
        
        p4_h1   = fMTT.tau1P4
        px_h1   = ak.from_regular(ak.Array(p4_h1[:,0:1]))
        py_h1   = ak.from_regular(ak.Array(p4_h1[:,1:2]))
        pz_h1   = ak.from_regular(ak.Array(p4_h1[:,2:3]))
        energy_h1 = ak.from_regular(ak.Array(p4_h1[:,3:4]))
        p4_h1_reg = ak.zip(
            {
                "x": px_h1, "y": py_h1, "z": pz_h1, "t": energy_h1,
            },
            with_name="LorentzVector",
            behavior=coffea.nanoevents.methods.vector.behavior,
        )
    
        p4_h2  = fMTT.tau2P4
        px_h2   = ak.from_regular(ak.Array(p4_h2[:,0:1]))
        py_h2   = ak.from_regular(ak.Array(p4_h2[:,1:2]))
        pz_h2   = ak.from_regular(ak.Array(p4_h2[:,2:3]))
        energy_h2 = ak.from_regular(ak.Array(p4_h2[:,3:4]))
        p4_h2_reg = ak.zip(
            {
                "x": px_h2, "y": py_h2, "z": pz_h2, "t": energy_h2,
            },
            with_name="LorentzVector",
            behavior=coffea.nanoevents.methods.vector.behavior,
        )

        
    elif run_fmtt_kinfit:
        is_3pr = events.zcand.decayMode >= 10
        tau_3pr = events.zcand[is_3pr]
        tau_other = events.zcand[~is_3pr]


        algo = "constant" # "gaussian", "constant", "bw", "exact"
        results = KinFit.kinfit_3pr(
            # First visible tau: must be the 3-prong tau
            as_float64(tau_3pr.pt),
            as_float64(tau_3pr.eta),
            as_float64(tau_3pr.phi),
            as_float64(tau_3pr.mass),
            # Second visible tau
            as_float64(tau_other.pt),
            as_float64(tau_other.eta),
            as_float64(tau_other.phi),
            as_float64(tau_other.mass),
            # MET
            as_float64(events.PuppiMET.pt[:,None]),
            as_float64(events.PuppiMET.phi[:,None]),
            # MET covariance
            as_float64(events.PuppiMET.covXX[:,None]),
            as_float64(events.PuppiMET.covXY[:,None]),
            as_float64(events.PuppiMET.covYY[:,None]),
            # PV-to-SV displacement vector of the 3-prong tau
            as_float64(tau_3pr.SVx - events.PVBS.x),
            as_float64(tau_3pr.SVy - events.PVBS.y),
            as_float64(tau_3pr.SVz - events.PVBS.z),
            # Secondary-vertex covariance
            as_float64(tau_3pr.SVcovxx),
            as_float64(tau_3pr.SVcovxy),
            as_float64(tau_3pr.SVcovxz),
            as_float64(tau_3pr.SVcovyy),
            as_float64(tau_3pr.SVcovyz),
            as_float64(tau_3pr.SVcovzz),
            # Steering parameters
            True,      # use full SV covariance
            mass_constraint=algo,
            mX=Z_mass,
            widthX=Z_width[algo],
        )

        px_3pr = ak.Array(results['px_1'])
        py_3pr = ak.Array(results['py_1'])
        pz_3pr = ak.Array(results['pz_1'])
        px_other = ak.Array(results['px_2'])
        py_other = ak.Array(results['py_2'])
        pz_other = ak.Array(results['pz_2'])
        chi_2  = ak.Array(results['chi2'])
        chi2_sv = ak.Array(results['chi2_sv'])
        chi2_met = ak.Array(results['chi2_met'])

        
        is_3pr_1st = ak.fill_none(ak.firsts(is_3pr, axis=1), False)

        px_1 = ak.where(is_3pr_1st, px_3pr, px_other)[:,None]
        py_1 = ak.where(is_3pr_1st, py_3pr, py_other)[:,None]
        pz_1 = ak.where(is_3pr_1st, pz_3pr, pz_other)[:,None]
        px_2 = ak.where(is_3pr_1st, px_other, px_3pr)[:,None]
        py_2 = ak.where(is_3pr_1st, py_other, py_3pr)[:,None]
        pz_2 = ak.where(is_3pr_1st, pz_other, pz_3pr)[:,None]
        
        mass_1 = events.zcand.mass[:,0:1]
        mass_2 = events.zcand.mass[:,1:2]

        energy_1 = np.sqrt(px_1**2 + py_1**2 + pz_1**2 + mass_1**2)
        energy_2 = np.sqrt(px_2**2 + py_2**2 + pz_2**2 + mass_2**2)        

        p4_h1_reg = ak.zip(
            {
                "x": px_1, "y": py_1, "z": pz_1, "t": energy_1,
            },
            with_name="LorentzVector",
            behavior=coffea.nanoevents.methods.vector.behavior,
        )
        p4_h2_reg = ak.zip(
            {
                "x": px_2, "y": py_2, "z": pz_2, "t": energy_2,
            },
            with_name="LorentzVector",
            behavior=coffea.nanoevents.methods.vector.behavior,
        )

    elif run_classic_svfit:
        """
        
        """
        
        # decay_type:
        #0 - TauToEle
        #1 - TauToMu
        #2 - TauToHad
        type_1 = ak.ones_like(hmass)
        type_0 = type_1 - 1
        type_2 = 2 * type_1
        
        
        type1  = ak.to_numpy(ak.where(events.channel_id == etau_id,
                                      type_0,
                                      ak.where(events.channel_id == mutau_id,
                                               type_1,
                                               type_2)
                                      )
                             )
        type2  = ak.to_numpy(type_2)
        
        
        lep1 = events.zcand[:,0:1]
        lep2 = events.zcand[:,1:2]
        
        results = classicsvfit_cpp.classic_svfit(
            as_float64(lep1.pt),
            as_float64(lep1.eta),
            as_float64(lep1.phi),
            as_float64(lep1.mass),
            as_int32(type1),
            as_float64(lep2.pt),
            as_float64(lep2.eta),
            as_float64(lep2.phi),
            as_float64(lep2.mass),
            as_int32(type2),
            as_float64(events.PuppiMET.pt[:,None]),
            as_float64(events.PuppiMET.phi[:,None]),
            as_float64(events.PuppiMET.covXX[:,None]),
            as_float64(events.PuppiMET.covXY[:,None]),
            as_float64(events.PuppiMET.covYY[:,None]),
            decay_mode_2=as_int32(lep2.decayMode),
            mX=-1.0, # no resonance-mass constraint
            max_calls=100000,
        )
        
        valid  = ak.Array(results['valid'])
        status = ak.Array(results['status'])
        px_1   = ak.Array(results['px_1'])[:,None]
        py_1   = ak.Array(results['py_1'])[:,None]
        pz_1   = ak.Array(results['pz_1'])[:,None]
        energy_1 = ak.Array(results['energy_1'])[:,None]
        px_2   = ak.Array(results['px_2'])[:,None]
        py_2   = ak.Array(results['py_2'])[:,None]
        pz_2   = ak.Array(results['pz_2'])[:,None]
        energy_2 = ak.Array(results['energy_2'])[:,None]
        mass   = ak.Array(results['mass'])[:,None]
        mass_err = ak.Array(results['mass_err'])[:,None]
        mass_from_taus = ak.Array(results['mass_from_taus'])[:,None]
        cpu_time = round(float(np.sum(results['cpu_time'])), 3)
        real_time = round(float(np.sum(results['real_time'])), 3)

        p4_h1_reg = ak.zip(
            {
                "x": px_1, "y": py_1, "z": pz_1, "t": energy_1,
            },
            with_name="LorentzVector",
            behavior=coffea.nanoevents.methods.vector.behavior,
        )
        p4_h2_reg = ak.zip(
            {
                "x": px_2, "y": py_2, "z": pz_2, "t": energy_2,
            },
            with_name="LorentzVector",
            behavior=coffea.nanoevents.methods.vector.behavior,
        )
        logger.info(f"Time elapsed by Classic-SVfit : CPU - {cpu_time} & Real - {real_time} seconds")
        
    else:
        logger.critical("Set [cfg.x.run_fastMTT = True] in the main config to make fastMTT run ... setting dummy variables as output")


    #from IPython import embed; embed()
    
    return p4_h1_reg,p4_h2_reg
