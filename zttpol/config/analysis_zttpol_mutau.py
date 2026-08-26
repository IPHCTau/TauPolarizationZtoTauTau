# coding: utf-8

"""
Configuration of the CPinHToTauTau analysis.
"""
import law
import order as od
from scinum import Number

# ------------------------ #
# The main analysis object #
# ------------------------ #
analysis_zttpol_mutau = ana = od.Analysis(
    name="analysis_zttpol_mutau",
    id=3,
)

# analysis-global versions
# (see cfg.x.versions below for more info)
ana.x.versions = {}

# files of bash sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.bash_sandboxes = [
    "$CF_BASE/sandboxes/cf.sh",
    "$ZTTPOL_BASE/sandboxes/venv_zttpol.sh",
]
#default_sandbox = law.Sandbox.new(law.config.get("analysis", "default_columnar_sandbox"))
#if default_sandbox.sandbox_type == "bash" and default_sandbox.name not in ana.x.bash_sandboxes:
#    ana.x.bash_sandboxes.append(default_sandbox.name)

# files of cmssw sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.cmssw_sandboxes = [
    #"$CF_BASE/sandboxes/cmssw_default.sh",
]

# config groups for conveniently looping over certain configs
# (used in wrapper_factory)
ana.x.config_groups = {}

# ------------- #
# setup configs #
# ------------- #

from zttpol.config.build_analysis import *
build_analysis(analysis=analysis_zttpol_mutau,
               era=2018,
               postfix="",
               channel="mutau",
               islimited=True,
               isfull=False)
