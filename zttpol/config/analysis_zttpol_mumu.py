# coding: utf-8

"""
Configuration of the CPinHToTauTau analysis.
"""
import luigi
import law
import order as od
from scinum import Number

# ------------------------ #
# The main analysis object #
# ------------------------ #
analysis_zttpol_mumu = ana = od.Analysis(
    name="analysis_zttpol_mumu",
    id=6,
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
default_sandbox = law.Sandbox.new(law.config.get("analysis", "default_columnar_sandbox"))
if default_sandbox.sandbox_type == "bash" and default_sandbox.name not in ana.x.bash_sandboxes:
    ana.x.bash_sandboxes.append(default_sandbox.name)

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

islimited, isfull = get_cfg_type_from_cmd(luigi.cmdline_parser.CmdlineParser.get_instance(),
                                          law.parser.root_task_parser())

build_analysis(analysis=analysis_zttpol_mumu,
               era=2018,
               postfix="",
               channel="mumu",
               islimited=islimited,
               isfull=isfull)
