import os
from typing import Optional
from columnflow.util import maybe_import

from zttpol.production.PolarimetricA1 import PolarimetricA1

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")

import law
logger = law.logger.get_logger(__name__)



def remove_empty_places(array, debug=False):
    counts = ak.num(array)
    max_count = int(np.max(counts))
    if debug == True:
        logger.info(f"DEBUG -- Max count : {max_count}")
        logger.info(f"DEBUG -- {ak.sum(counts == max_count, axis=0)} valid entries out of {ak.num(counts, axis=0)}")
    flat = ak.flatten(array)
    out = ak.unflatten(flat, max_count)
    return out, counts


def back_to_original(array, counts, debug=False):
    array = ak.flatten(array)
    out = ak.unflatten(array, counts)
    if debug == True:
        counts = ak.num(out)
        yes = ak.sum(counts > 0, axis=0)
        total = ak.num(counts, axis=0)
        logger.info(f"DEBUG -- {yes} valid entries among {total}")
    return out


def clean_rearranged_dict(p4dict, keylist=[], debug=False):

    if len(keylist) == 0:
        raise RuntimeWarning("must provide some key")
    
    count = None
    out = {key:None for key in p4dict}
    for key in keylist:
        val = p4dict[key]
        val, count = remove_empty_places(val, debug=debug)
        out[key] = ak.from_regular(val)

    if count is None:
        raise RuntimeWarning('check clean_rearranged_dict')

    count = ak.values_astype(ak.where(count > 1, 1, count), np.int32)
    mask = count > 0
    
    for okey in out:
        if okey not in keylist:
            masked_val = p4dict[okey][mask]
            max_count_safe = ak.max(ak.num(masked_val, axis=1))
            min_count_safe = ak.min(ak.num(masked_val, axis=1))
            if max_count_safe == min_count_safe:
                masked_val = ak.drop_none(masked_val)
            else:
                logger.warning(f"drop_none won't probably be safe. Check again the array for {okey}")
            out[okey] = masked_val
    
    return out, count
            


def getlistofobservables():
    variables = ["costheta_1", "costheta_2",
                 "cosalpha_1", "cosalpha_2",
                 "omega_1",    "omega_2",
                 "omegavis_1", "omegavis_2",
                 "omegabar_1", "omegabar_2",
                 "OMEGA", "OMEGAVIS", "OMEGABAR",
                 "massvis"]
    return variables



def check_nan(var, val):
    n_nan = ak.sum(np.isnan(val))
    if n_nan > 0:
        val = ak.nan_to_num(val, -999.9)
        n_tot = ak.sum(np.abs(val) >= 0)
        logger.warning(f"{var} --> {n_nan} out of {n_tot} is NaN")
        
    n_none = ak.sum(ak.is_none(val, axis=-1))
    if n_none > 0:
        logger.warning(f"{var} --> {n_none} out of {n_tot} is None")
    val = ak.fill_none(val, -999.9)
        
    return val
    

def wrap(mask, src, dest):
    for key,val in src.items():
        if val is None:
            continue
        if key not in dest:
            raise RuntimeWarning(f'{key} not found in final dictionary')

        val = check_nan(key, val)
        dest[key] = ak.where(mask, val, dest[key])
    return dest



def unwrap(vardict, counts, debug=False):
    for key, val in vardict.items():
        if val is None:
            continue
        vardict[key] = back_to_original(val, counts, debug=debug)

    return vardict


"""
def to_cartesian(p4):
    pt = p4.pt
    eta = p4.eta
    phi = p4.phi
    mass = p4.mass
    charge = p4.charge if "charge" in p4.fields else ak.zeros_like(mass)
    
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    E = np.sqrt(mass**2 + pt**2 * (np.cosh(eta)**2))

    cartesian_p4 = ak.zip(
        {
            'x': ak.values_astype(ak.where(px == np.nan, -999.0, px), np.float64),
            'y': ak.values_astype(ak.where(py == np.nan, -999.0, py), np.float64),
            'z': ak.values_astype(ak.where(pz == np.nan, -999.0, pz), np.float64),
            't': ak.values_astype(ak.where(E == np.nan, -999.0, E), np.float64),
            'charge': charge,
        },
        with_name="LorentzVector",
        behavior=coffea.nanoevents.methods.vector.behavior
    )
    
    return cartesian_p4
"""

def to_cartesian(p4):
    pt = p4.pt
    eta = p4.eta
    phi = p4.phi
    mass = p4.mass
    charge = p4.charge if "charge" in p4.fields else ak.zeros_like(mass)

    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    energy = np.sqrt(mass**2 + pt**2 * np.cosh(eta)**2)

    return ak.zip(
        {
            "x": ak.nan_to_none(ak.values_astype(px, np.float64)),
            "y": ak.nan_to_none(ak.values_astype(py, np.float64)),
            "z": ak.nan_to_none(ak.values_astype(pz, np.float64)),
            "t": ak.nan_to_none(ak.values_astype(energy, np.float64)),
            "charge": charge,
        },
        with_name="LorentzVector",
        behavior=coffea.nanoevents.methods.vector.behavior
    )


#def rotate(vec, rot):
#    angle_z = 0.5 * np.pi - rot.phi
#    angle_x = rot.theta
#    rot = vec.rotateZ(angle_z).rotateX(angle_x)
#    print(rot.x)
#    print(rot.y)
#    print(rot.z)
#    return rot


def rotate(vec, rot):
    """
    Equivalent to C++:
        vec.RotateZ(0.5*pi - Rot.Phi())
        vec.RotateX(Rot.Theta())

    vec: LorentzVectorArray with x,y,z,t
    rot: ThreeVectorArray with x,y,z
    """

    angle_z = 0.5 * np.pi - rot.phi
    angle_x = rot.theta

    x = vec.x
    y = vec.y
    z = vec.z
    t = vec.t

    # RotateZ(angle_z)
    cz = np.cos(angle_z)
    sz = np.sin(angle_z)

    x1 = cz * x - sz * y
    y1 = sz * x + cz * y
    z1 = z

    # RotateX(angle_x)
    cx = np.cos(angle_x)
    sx = np.sin(angle_x)

    x2 = x1
    y2 = cx * y1 - sx * z1
    z2 = sx * y1 + cx * z1

    fields = {
        "x": x2,
        "y": y2,
        "z": z2,
        "t": t,
    }

    if "charge" in vec.fields:
        fields["charge"] = vec.charge

    rot = ak.zip(
        fields,
        with_name="LorentzVector",
        behavior=coffea.nanoevents.methods.vector.behavior,
    )

    return rot


def getPolarimetricVector(tauP4    = None,
                          pionP4   = None,
                          pizeroP4 = None,
                          #boostv   = None,
                          tauch    = None,
                          leg = ""):
    """
    Everything must be in some rest frame
    """
    h  = None
    if leg == "pi":

        #h = pionP4.boost(boostv.negative()).pvec
        h = pionP4.pvec
        
    elif leg == "rho":

        P   = tauP4 #.boost(boostv.negative())
        pi  = pionP4 #.boost(boostv.negative())
        pi0 = pizeroP4 #.boost(boostv.negative())
        q   = pi.subtract(pi0)
        N   = P.subtract(pi.add(pi0))

        # From IPHC Polarization code
        ##### R = P.M()*(2*(q*N)*q.Vect() - q.Mag2()*N.Vect()) * (1/ (2*(q*N)*(q*P) - q.Mag2()*(N*P)))
        qN = q.dot(N)
        qP = q.dot(P)
        NP = N.dot(P)
        q2 = q.mass2
        num = 2 * qN * q.pvec - q2 * N.pvec
        den = 2 * qN * qP - q2 * NP
        h = P.mass * num / den
        #h = ( 2*  ( ( q.energy*N.energy ) - (q.dot(N)) )*q.pvec)  - (q.mass2*N.pvec)
        
        
    elif leg == "a1":
        
        P       = tauP4 #.boost(boostv.negative())
        os_pi   = pionP4[:, 0:1] #.boost(boostv.negative())
        ss1_pi  = pionP4[:, 1:2] #.boost(boostv.negative()) 
        ss2_pi  = pionP4[:, 2:3] #.boost(boostv.negative()) 
        a1pol   = PolarimetricA1(P,
                                 os_pi,
                                 ss1_pi,
                                 ss2_pi,
                                 tauch)
        h = -a1pol.PVC().pvec
        

    else:
        raise RuntimeError(f"Wrong leg assignment : {leg}")
        
    return h


    

def getCombOMEGA(omega1 : ak.Array,
                 omega2 : ak.Array,
                 **kwargs):
    
    debug = kwargs.get('debug', False)

    omega = (omega1 + omega2)/(1 + omega1*omega2)
    nan_omega = np.isnan(omega)
    inf_omega = ~np.isfinite(omega)
    bad_omega = (nan_omega | inf_omega)
    if (debug == True) and ak.any(bad_omega):
        logger.warning(f'comb : {ak.sum(nan_omega)} NaNs and {ak.sum(inf_omega)} Idfinites')
    omega = ak.where(bad_omega, -999.0, omega)
    
    return omega



