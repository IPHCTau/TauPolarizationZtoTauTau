# Recommended : two patches need to be used to use IPHC slurm. Easy way to apply.


## `job.py`

 - `cd ../modules/columnflow/modules/law/src/law/contrib/slurm`
 - `mv job.py job.py.default`
 - `cd -`
 - `cp job.py ../modules/columnflow/modules/law/src/law/contrib/slurm`


## `remote.py`

 - `cd ../modules/columnflow/columnflow/tasks/framework`
 - `mv remote.py remote.py.default`
 - `cd -`
 - `cp remote.py ../modules/columnflow/columnflow/tasks/framework`