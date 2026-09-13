# TauPolarizationZtoTauTau Analysis in IPHC computing cluster

Analysis is now moved from CERN lxplus to IPHC UI --> So, CERN HTCondor to IPHC Slurm

## IMPORTANT:
 - `remote.py` is now merged to `columnflow` [PR #820](https://github.com/columnflow/columnflow/pull/820)
 - `job.py` in `law` is pending [PR #212](https://github.com/riga/law/pull/212)
If the later is merged to law, columnflow with updated law submodule can be pulled safely, and copying these two files won't be neecessary.



## Resources

 - `columnflow` : [origin](https://github.com/IPHCTau/columnflow/tree/865e970d4a20b87b9020687a616be9cc3dbb879b) -- [upstream](https://github.com/columnflow/columnflow/)
 - `cmsdb` : [origin](https://github.com/IPHCTau/cmsdb/tree/04bcb41b23a3a68423026f859d2ea42bce494102) -- [upstream](https://github.com/uhh-cms/cmsdb)
 - [law](https://github.com/riga/law)
 - [order](https://github.com/riga/order)
 - [luigi](https://github.com/spotify/luigi)


## To setup the analysis

 - Clone the repo: `git clone --recurse-submodules <git-id>`
 - Before setting up the environment, check `README` inside `custom_patches`
   - copy `remote.py` and `job.py` as mentioned there. 
 - Use `slurm_setup.sh` --> `source slurm_setup.sh`
   - it will build a new environment named `zttenv`
     - Use ` CF_DATA="/slurm/shared/cms/<username>/TauPolarizationZtoTauTauData"`
   - it will also fix all path related issues
   - make sure it shows expected information at the end after it install


## use `--workflow slurm` to run jobs in `slurm`