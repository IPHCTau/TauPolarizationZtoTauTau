# coding: utf-8

"""
Main categories file for the Higgs CP analysis
"""

from columnflow.categorization import Categorizer, categorizer
from columnflow.util import maybe_import

ak = maybe_import("awkward")


# ---------------------------------------------------------- #
#                          Channels                          #
# ---------------------------------------------------------- #

# inclusive
@categorizer(uses={"event"})
def cat_incl(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ak.ones_like(events.event) == 1

# mutau
@categorizer(uses={"channel_id"})
def cat_mutau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("mutau")
    return events, events["channel_id"] == ch.id

# ---------------------------------------------------------- #
#                          For nJets                         #
# ---------------------------------------------------------- #
@categorizer(uses={"has_0jet"})
def cat_0j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.has_0jet

@categorizer(uses={"has_1jet"})
def cat_1j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.has_1jet

@categorizer(uses={"has_2jet"})
def cat_2j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.has_2jet

# ---------------------------------------------------------- #
#                          For PhiCP                         #
# ---------------------------------------------------------- #


# ---------- >>> for e/mu-tauh
# tau -> pi
@categorizer(uses={"is_pi_2"})
def cat_pi_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_pi_2
# tau -> rho
@categorizer(uses={"is_rho_2"})
def cat_rho_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_rho_2
# tau -> a1 (DM2)
@categorizer(uses={"is_a1_1pr_2pi0_2"})
def cat_a1dm2_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_a1_1pr_2pi0_2
# tau -> a1 (DM10)
@categorizer(uses={"is_a1_3pr_0pi0_2"})
def cat_a1dm10_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_a1_3pr_0pi0_2
# tau -> a1 (DM11)
@categorizer(uses={"is_a1_3pr_1pi0_2"})
def cat_a1dm11_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_a1_3pr_1pi0_2

# IPSig
@categorizer(uses={"is_ipsig_0to1_1"})
def cat_ipsig_0to1_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_ipsig_0to1_1
@categorizer(uses={"is_ipsig_0to1_1"})
def cat_ipsig_1toany_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ~events.is_ipsig_0to1_1




# to know true or fake tau
@categorizer(uses={"is_real_1"})
def cat_real_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_real_1
@categorizer(uses={"is_fake_1"})
def cat_fake_1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_fake_1
@categorizer(uses={"is_real_2"})
def cat_real_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_real_2
@categorizer(uses={"is_fake_2"})
def cat_fake_2(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_fake_2




## --- e/mutau --->>>

#        bveto      no-bveto                 bveto                 
# << ------------ | --------- | -------------------------------- >>
# |-------------------------------------------------------------- | ^
# |               |           |           |                       | ^
# |               |           |           |                       | | tau-2
# |      QCD      |     T     |     SR    |          W            | |  iso
# |       A       |    A0     |     D     |          A1           | |
# |               |           |           |                       | |
# |-------------------------------------------------------------- | -
# |               |           |           |                       | |
# |               |           |           |                       | |
# |      QCD      |     T     |     SR    |          W            | |  tau-2
# |       B       |    B0     |     C     |          B1           | | antiIso
# |               |           |           |                       | |
# |-------------------------------------------------------------- | v
# << ------------ | -------------------------------------------- >> v
#        SS                           OS                           
# << ------------------------------------ | -------------------- >>
#                  mT < 50                          mT > 50        
#
# A   : e/mutau [ss__iso2__nobjet__lowmt]
# B   : e/mutau [ss__noniso2__nobjet__lowmt]
# A0  : e/mutau [os__iso2__bjet__lowmt]
# B0  : e/mutau [os__noniso2__bjet__lowmt]
# A1  : e/mutau [os__iso2__nobjet__highmt]
# B1  : e/mutau [os__noniso2__nobjet__highmt]
# D   : e/mutau [os__iso2__nobjet__lowmt]
# C   : e/mutau [os__noniso2__nobjet__lowmt]

# A
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_ss_iso2_bveto_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ~events.is_os & events.is_iso_2 & events.is_b_veto & events.is_low_mt 
# B
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_ss_noniso2_bveto_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ~events.is_os & ~events.is_iso_2 & events.is_b_veto & events.is_low_mt 
# A0
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_os_iso2_nobveto_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & events.is_iso_2 & ~events.is_b_veto & events.is_low_mt 
# B0
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_os_noniso2_nobveto_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & ~events.is_iso_2 & ~events.is_b_veto & events.is_low_mt 
# A1
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_os_iso2_bveto_highmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & events.is_iso_2 & events.is_b_veto & ~events.is_low_mt 
# B1
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_os_noniso2_bveto_highmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & ~events.is_iso_2 & events.is_b_veto & ~events.is_low_mt 
# D
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_os_iso2_bveto_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & events.is_iso_2 & events.is_b_veto & events.is_low_mt 
# C
@categorizer(uses={"is_os", "is_iso_2", "is_b_veto", "is_low_mt"})
def cat_os_noniso2_bveto_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & ~events.is_iso_2 & events.is_b_veto & events.is_low_mt


# DESY (Stepan)
# SR
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_low_mt"})
def cat_os_iso1_iso2_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & events.is_iso_1 & events.is_iso_2 & events.is_low_mt 
# AR
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_low_mt"})
def cat_ss_iso1_iso2_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ~events.is_os & events.is_iso_1 & events.is_iso_2 & events.is_low_mt 
# DR_Num
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_low_mt"})
def cat_os_noniso1_iso2_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, events.is_os & ~events.is_iso_1 & events.is_iso_2 & events.is_low_mt 
# DR_Den
@categorizer(uses={"is_os", "is_iso_1", "is_iso_2", "is_low_mt"})
def cat_ss_noniso1_iso2_lowmt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    return events, ~events.is_os & ~events.is_iso_1 & events.is_iso_2 & events.is_low_mt
