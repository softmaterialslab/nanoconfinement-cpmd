# nanoconfinement-cpmd

## Install instructions on BigRed2
* First, git clone the project
```git clone https://github.com/softmaterialslab/nanoconfinement-cpmd```
* Then, load the required modules using following command:
```module swap PrgEnv-cray PrgEnv-gnu && module load boost/1.65.0 && module load gsl```
* Next, go to the root directory:
 ```cd nanoconfinement-cpmd```
* Then, install the project:
```make cluster-install```
* Fianlly, submit the job:
```make cluster-submit```
* All outputs from the simulation will be stored in the bin/outfiles folder when the simulation is completed.

