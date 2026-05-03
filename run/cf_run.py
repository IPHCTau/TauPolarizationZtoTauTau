import os
import yaml
import argparse 
import datetime
import time
import sys


#from utilities import setup_logger, runShellCmd, run_command
datetime_tag = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
ListOfMainFuncs = ['CalibrateEvents',
                   'SelectEvents',
                   'ReduceEvents',
                   'PlotCutFlow',
                   'ProduceColumns',
                   'CreateHistograms','MergeHistograms',
                   'PlotVariables1D', 'PlotVariables2D', 'PlotVariablesPerProcess2D',
                   'PlotShiftedVariables1D', 'PlotShiftedVariables2D',
                   'UniteColumns']

parser = argparse.ArgumentParser(description='Pre-Processing')
parser.add_argument('-i',
                    '--inconfig',
                    type=str,
                    required=True,
                    help="cf.SelectEvents arguments")
parser.add_argument('-f',
                    '--func',
                    type=str,
                    required=True,
                    help="e.g. SelectEvents")
parser.add_argument('-c',
                    '--channel',
                    type=str,
                    required=True,
                    help="e.g. emu / etau / mutau / tautau")
parser.add_argument('-l',
                    '--limited',
                    action='store_true',
                    required=False,
                    default=False,
                    help="full or limited")

pargs = parser.parse_args()

yml_config = None
yml_config_file = pargs.inconfig

with open(yml_config_file,'r') as conf:
    yml_config = yaml.safe_load(conf)

main_func = pargs.func
assert main_func in ListOfMainFuncs, f"Error: {main_func} is wrong"

islimited = pargs.limited
yml_config = yml_config.get("FULL") if not islimited else yml_config.get("LIMITED")


era           = yml_config.get("era")
run           = yml_config.get("run")
postfix       = yml_config.get("postfix")
limit         = "limited" if islimited else "full"
nworkers      = yml_config.get("workers")
tasks_per_job = yml_config.get("tasks_per_job")
wrapper       = yml_config.get("wrapper")

main_args = yml_config.get("args")
version_tag = main_args.get("version")
if version_tag == "dummy":
    version_tag = f"Run{run}_{era}{postfix}_{limit}_{datetime_tag}"

print(f"Yaml     : {yml_config_file}")
print(f"Execute  : cf.{main_func}")
print(f"Era      : {era}")
print(f"Run      : {run}")
print(f"Postfix  : {postfix}")
print(f"Limited? : {islimited}")
print(f"nWorkers : {nworkers}")
print(f"nTasks_per_job : {tasks_per_job}")
print(f"Wrapper  : {wrapper}")

channel = pargs.channel
print(f"Channel : {channel}")


config    = main_args.get("config")
config    = config if not islimited else f"{config}_limited" 
print(f"Config   : {config}")

dataset_list = main_args.get("datasets")
if len(dataset_list) > 1 and wrapper == False:
    print(f"\tmore than one dataset, making wrapper True")
    wrapper = True
datasets  = ",".join(dataset_list)
print(f"datasets : {datasets}")

processes = ",".join(main_args.get("processes"))
print(f"processes: {processes}")

shiftList = main_args.get("shifts")
shifts = ",".join(shiftList) if shiftList else ""
print(f"shifts: {shifts}")

workflow  = main_args.get("workflow")
print(f"worklflow: {workflow}")

branch    = main_args.get("branch")
print(f"branch   : {branch}")

version = version_tag
print(f"version  : {version}")

extras = []
if "extras" in list(main_args.keys()):
    extras = [f"--{extra}" for extra in main_args.get("extras")]

categories = ",".join(main_args.get("categories")) if "categories" in list(main_args.keys()) else None
variables  = ",".join(main_args.get("variables")) if "variables" in list(main_args.keys()) else None
variables1D = ",".join([var for var in main_args.get("variables") if len(var.split('-')) < 2])
variables2D = ",".join([var for var in main_args.get("variables") if len(var.split('-')) == 2])

shifts_val = "nominal,{"+f"{shifts}"+"}_{up,down}" if shiftList else "nominal"

skip_wrapping = main_func.startswith("PlotVariables")
if wrapper:    
    if not skip_wrapping:
        cmd_list = [
            "law", "run", f"cf.{main_func}Wrapper",
            "--config", config,
            "--datasets", datasets,
            f"--shifts", shifts_val,
            f"--cf.{main_func}-workflow", workflow,
            #f"--cf.{main_func}-branch", branch,
            f"--cf.{main_func}-tasks-per-job", tasks_per_job,
            f"--cf.{main_func}-pilot",
            "--version", version,
            "--workers", nworkers,
            #f"--cf.{main_func}-log-file", jobfile,
            f"--analysis zttpol.config.analysis_zttpol_{channel}.analysis_zttpol_{channel}"
        ]
        if main_func.startswith("CreateHistograms") or main_func.startswith("MergeHistograms"):
            cmd_list.append(f"--cf.{main_func}-variables")
            cmd_list.append(variables)
    else:
        _variables = variables1D if main_func == "PlotVariables1D" else variables2D
        print(f"variables : {_variables}")
        cmd_list = [
            "law", "run", f"cf.{main_func}",
            "--config", config,
            "--datasets", datasets,
            "--processes", processes,
            f"--workflow", workflow,
            #f"--branch", branch,
            "--version", version,
            "--categories", categories,
            "--variables", _variables,
            "--workers", nworkers,
            "--pilot",
            #f"--log-file", jobfile,
            f"--analysis zttpol.config.analysis_zttpol_{channel}.analysis_zttpol_{channel}"
        ] + extras

else:
    cmd_list = [
        "law", "run", f"cf.{main_func}",
        "--config", config,
        "--dataset", datasets,
        "--workflow", workflow,
        "--workers", nworkers,
        #"--branch", branch,
        "--version", version,
        "--pilot",
        #"&>", jobfile, '&'
        f"--analysis zttpol.config.analysis_zttpol_{channel}.analysis_zttpol_{channel}"
    ]
    if branch > 0:
        cmd_list.append("--branch")
        cmd_list.append(branch)
    if main_func in ["CreateHistograms", "MergeHistograms", "PlotVariables1D", "PlotVariables2D"]:
        cmd_list.append("--variables")
        if main_func in ["PlotVariables1D", "PlotVariables2D"]:
            cmd_list.append(variables1D) if "1D" in main_func else cmd_list.append(variables2D)
            cmd_list.append("--categories")
            cmd_list.append(categories)
        else:
            cmd_list.append(variables)

    
        
cmd_list = [str(item) for item in cmd_list]
cmd = " ".join(cmd_list)

print(f"CMD      : \n\n{cmd}\n")
print(f"Copy this line above, paste it and hit enter : BABUSHCHA !")

#print(f"writing CMD in {cmdfile}")
#cmdf = open(cmdfile, 'w')
#cmdf.write("#!/bin/sh\n\n")
#cmdf.write(cmd + '\n')
#cmdf.close()
