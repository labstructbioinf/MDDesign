# **MDDesign** #

Set of scripts that facilitate molecular dynamics aided screening of the models generated with Rosetta and selection of the structural ensembles that can be used as input in protein design simulations.

## **Requirements**

* Installed AmberTools 16 and properly set Amber environment (usually "*source /path/to/amber16/amber.sh*" is enough)
* BioPDB (for PDB manipulation)
* DSSP (for secondary structure assignment)
* Imagemagick (and its "convert" executable set in path)

## **Installation**

Clone the repository. This will create the mdscripts directory.
```
$ git clone https://jludw@bitbucket.org/jludw/mddesign.git
```
Now add the whole directory to PATH.

```
$ export PATH=$PATH:/path/to/mdscripts
```

## **MDPrep.py**
Prepares necessary files (directory tree) and inputs for the MD simulations.

Preferably a new directory should be created with the analyzed PDB files (with .pdb extension) or (solvated!) topologies and restarts in format MODEL_NAME_i_s.parm7 MODEL_NAME_i_s.rst7.

Arguments:

| Option    | Description |
|:----------:|-------------|
| **`-start`** | Required if script is used in the directory for the first time. This will process all the pdb files in directory, create Amber topologies and restarts as well as process inputs and generate job files for the first run. |
| **`-list`** | Name of the filename of the models list included in simulations / analyses. If used along with -start then the file will be created based on the pdb files in directory, otherwise models will be read from the supplied list. |
| **`-rep`** | Used only with `-start` - states the number of individual replicates started with different velocity distribution.|
| **`-i`** | Amber input for the certain step. |
| **`-one`** | Used if only one trajectory is to be calculated prior to the assignnment of different velocities. |
| **`-copy`** | Used after calculation with `-one` parameter. This will copy the files resulting from one job to the rest of subdirectories of skipped replicas. |
| **`-cmd`**    | Command to execute amber, e.g. 'pmemd.cuda', 'pmemd.MPI', 'mpirun -np 24 pmemd.MPI' etc., default 'pmemd.cuda' |
| **`-ref`**    | Used if reference structure for Amber restraints is required - currently always the restart structure from the previous step is used. |
| **`-nox`**    | Do not output trajectory. |
| **`-soll_wall_dist`** | Used only with `-start`. Determines the minimum solute wall distance for the topology preparation. Default 10. |
| **`-cluster_header`**    | Cluster queuing system header. For examples look at inputs directory. |
| **`-cluster_add_cmd`**    | Command to add batch job to the queue. Default is the slurms' sbatch. |

### Example workflow: ###
```
#!bash
$ cd /path/with/pdb/files/to/analyze

$ MDPrep.py -list MODEL_LIST -i min.in -start -one -cmd 'pmemd' -rep 50

$ ./STEP_1.sh

$ MDPrep.py -copy -step 1 -list MODEL_LIST

$ MDPrep.py -list MODEL_LIST -i heat1.in -one -ref -nox -step 2

$ ./STEP_2.sh

$ MDPrep.py -list MODEL_LIST -i heat2.in -one -ref -nox -step 3

$ ./STEP_3.sh

$ MDPrep.py -list MODEL_LIST -copy -step 3

$ MDPrep.py -list MODEL_LIST -i md1.in -step 4

$ ./STEP_4.sh
```

If pipeline is to be run on a cluster **`-cluster_header`** and **`-cluster_add_cmd`** have to be added. First one is a batch job header used by e.g. slurm queing system. All parameters have to be specified in this file except actual command to run the calculation which will be filled by the script automatically. Second command, **`-cluster_add_cmd`**, allows to specify the command which will add job to the queing system. Depending on the type of controller this can be 'sbatch', 'qsub', etc.

### Example cluster header file: ###
```
#!bash
#!/bin/bash
#SBATCH -p gpu
#SBATCH --exclusive

module load cuda7.5
source /home/amber16/amber.sh
export CUDA_HOME=/usr/local/cuda-7.5/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUDA_HOME/lib:$CUDA_HOME/lib64
export CUDA_VISIBLE_DEVICES=0
```

### Example commands to run the calculation on the cluster: ###

```
#!bash
$ cd /path/with/pdb/files/to/analyze

$ MDPrep.py -list MODEL_LIST -i min.in -start -one -cmd 'mpirun -np 24 pmemd.MPI' -rep 50 -cluster_add_cmd sbatch -cluster_header file.sh

$ ./STEP_1.sh

$ MDPrep.py -copy -step 1 -list MODEL_LIST

$ MDPrep.py -list MODEL_LIST -i heat1.in -one -ref -nox -step 2 -cluster_add_cmd sbatch -cluster_header file.sh

$ ./STEP_2.sh

$ MDPrep.py -list MODEL_LIST -i heat2.in -one -ref -nox -step 3 -cluster_add_cmd sbatch -cluster_header file.sh

$ ./STEP_3.sh

$ MDPrep.py -list MODEL_LIST -copy -step 3 

$ MDPrep.py -list MODEL_LIST -i md1.in -step 4 -cluster_add_cmd sbatch -cluster_header file.sh

$ ./STEP_4.sh
```
## **MDAnalyze.py**

Analyzes the output of the runs done with MDPrep.py. Run this only after the simulations are finished.

Arguments:

| Option    | Description |
|:----------:|-------------|
| **`-job`** | Job type. Available are: **autoimage** (autoimages all trajectories), **refstruct** (prepares reference structures for the RMSD calculation), **rmsd_bb** (calculates backbone RMSD throughout simulation), **rmsd_sse** (calculates RMSD over secondary structure elements throughout the simulations), **rmsd_custom** (custom RMSD calculation, -mask has to be specified), **secstruct** (secondary structure progression in simulation). |
| **`-list`** | Model list generated with MDPrep.py |
| **`-step`** | Steps number to perform analysis on (eg. '4' or '4 5')|
| **`-suffix`** | Suffix names for rmsd calculations to identify them and further plot with MDPlot.py. |
| **`-mask`** | Amber type mask (e.g. ':1-20') for custom rmsd analysis. |
| **`-nofit`** | Perform backbone rms first and then calculate rmsd of certain part of the protein. Useful for custom residue ranges. |
| **`-ai`** | Use the autoimaged trajectory for the analysis. Such trajectory must be first calculated as in the example workflow. |

### Example workflow: ###

```
#!bash

$ MDAnalyze.py -list MODEL_LIST -job autoimage -step 4

$ MDAnalyze.py -list MODEL_LIST -job refstruct -step 1

$ MDAnalyze.py -list MODEL_LIST -job rmsd_bb -step 4 -suffix rmsd_bb -ai

$ MDAnalyze.py -list MODEL_LIST -job rmsd_sse -step 4 -suffix rmsd_sse -ai

$ MDAnalyze.py -list MODEL_LIST -job secstruct -step 4 -ai

```

## **MDPlot.py **

Creates a visualization of the output of MDAnalyze.

| Option    | Description |
|:----------:|-------------|
| **`-list`** | Model list generated with MDPrep.py |
| **`-rmsd`** | List of suffixes (generated with MDAnalyze) to include in rmsd visulization. |
| **`-secstruct`** | Plot secondary structure time progression (must be calculated with MDAnalyze) |
| **`-out`** | Out PDF filename. |
| **`-dt`** | Determines every which picosecond frames are output from MD run. Default - 2. |
| **`-skip_first`** | Skip first n frames in analysis. |
| **`-skip_last`** | Skip last n frames in analysis. |
| **`-rmsd_min`** | Minimal RMSD values to include in the plots. |
| **`-rmsd_max`** | Maxmimal RMSD values to include in the plots. This allows to keep highly unstable variants with high RMSD values out of the plot drawing range. |
| **`-percentile`** | Adjust maximum visible RMSD values on the plots to the percentile of the distribution of all RMSD values from all simulations. Scale 0-1. |
| **`-format`** | Image format for plots drawing. Default - 'png' |
| **`-dpi`** | DPI of the resultant images. |
| **`-keep_images`** | Keep images after merging them into PDF. |

### Example workflow: ###

```
#!bash
$ MDPlot.py -list MODEL_LIST -rmsd rmsd_bb rmsd_sse -secstruct -skip_first 1500 -out FINAL.pdf

```
If problems with tiff images appear in matplotlib try:

```
#!bash

$ pip install --user pillow

```

## **MDCombClust.py**

Performs combined clustering on the trajectories to extract a diverse set of backbones and respective full atom structures visited during simulations.

| Option    | Description |
|:----------:|-------------|
| **`-list`** | Model list (can be manually adjusted to remove not stable structures) |
| **`-step`** | Step number from MDPrep, determines from which step trajectory is taken. |
| **`-stride`** | Every n-th frame from simulation is take for analysis. |
| **`-nstruct`** | Number of representative structures to output from clustering. |
| **`-start_frame`** | Start frame from each trajectory to include in analysis (otherwise first is used). 
| **`-last_frame`** | Last frame number of each simulation (required) |
| **`-out_dir`** | Output directory for PDB structures and statistics (default - COMBCLUST) |
| **`-ai`** | Use the autoimaged trajectory for the analysis. Such trajectory must be first calculated as in the example workflow. |

Warning! Clustering can take few hours and progress will be shown after run of the script. Also be careful with number of frames used in analysis. It is advisable to tweak -start_frame and -stride in such a way (depending on the number of models investigated) to have ~20k structures in clustering. Otherwise RAM problems may occur.
In the **`-out_dir`** two types of file will be available after calculation: **rep_bb** - aligned backbone representative structures and **rep_allatom** - reconstituted, based on frames with backbone representative structures, full atom conformations.

### Example workflow: ###

```
#!bash
$ MDCombClust.py -list MODEL_LIST -start_frame 1001 -last_frame 2000 -stride 40 -step 4 -nstruct 500 -ai

```