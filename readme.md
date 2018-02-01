# Objective

A distributed computing algorithm is developed to identify the crystalline and amorphous phases in a coarse-grained polymer system (coarse-grained linear polyethylene).
Distributed computing is separately implemented in both MPI (message passing interface) mode and multiprocessing mode. 
The bin-table technique is used in both modes to further improve computing efficiency.

This code was used in my project 6 and project 7, as demonstrated in [my webpage](https://liyiyang.weebly.com/).

The local p2 order parameter is computed for every particle in system, which is used to determine if a particle should be assigned to the crystalline phase.


### MPI mode

Run distributed computing on all the worker nodes in your cluster network.
The master node can also be utilized along side the worker nodes (this can be changed by modifying the 'host' file.)
 
Usage:
```
cd run
mpiexec -np <n> python ../src/main.py mpi
```
or
```
cd run
mpiexec --hostfile ./hosts -np <n> python ../src/main.py mpi
```
where <n> is the number of physical cores in your cluster.
To use the `--hostfile` option, you need to modify the 'host' file to accommodate the setting of your cluster.

### Multiprocessing mode

Run parallel computing only on the master node. Usage:
```
cd run
python ../src/main.py mp
```


# Prerequisites

##### MPI mode:

&nbsp;&nbsp;&nbsp;&nbsp;1. Open MPI or MPICH <br />
&nbsp;&nbsp;&nbsp;&nbsp;2. mpi4py

These should be installed on the master node as well as all the worker nodes in your cluster. <br />
Master node should be able to passwordlessly ssh to worker nodes. <br />

##### Multiprocessing mode: <br />

&nbsp;&nbsp;&nbsp;&nbsp;None


# Code structure

. <br />
├── src <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── main.py <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── modules <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── __init__.py <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── parser.py <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── md_system.py <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── bin_table.py <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── identify_mpi.py <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── identify_mp.py <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── p2_parameters.py <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── postprocess.py <br />
├── data <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── cg-80blk-50chain.lammps <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── cg-80blk-50chain-1atm-300K-40ns.lammpstrj <br />
├── run <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── hosts <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── input <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── image.png <br />
├── vmd_scripts <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── color_scale.vmd <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── setuserfield.vmd


##### src/main.py

Main script that initializes input parameters, and launches either MPI or multiprocessing jobs to identify crystalline particles in polymer system.

##### src/modules/parser.py

Script providing functions to read parameters from the input file 'input'.

##### src/modules/md_system.py

A class to store molecular information of a coarse-grained linear polyethylene system.

##### src/modules/bin_table.py

A class to generate a bin table for a polymer system.

##### src/modules/identify_mpi.py

Script using MPI to identify possible crystalline particles, utilizing all worker nodes in cluster.

##### src/modules/identify_mp.py

Script using multiprocessing to identify possible crystalline particles, utilizing only master node.

##### src/modules/p2_parameters.py

Script to calculate local p2 order parameter for each atom i in a list of local particles to be handled by a parallel process.

##### src/modules/postprocess.py

Script to run on master process, that completes the identification of crystalline particles.

##### data/cg-80blk-50chain.lammps

A LAMMPS data file for a reference state of a polymer system.

##### data/cg-80blk-50chain-1atm-300K-40ns.lammpstrj

A LAMMPS trajectory file containing two snap shots (or time steps) of the evolution of the polymer system defined by the LAMMPS data file.

##### run/hosts

Configuration file listing the nodes that will be used for distributed computing.

##### run/input

Input file containing required parameters.

##### run/image.png

An example randering of the polymer system with VMD. <br />
Left:  at 0.5 ns. <br />
Right: at 40.0 ns.

##### vmd_scripts/color_scale.vmd

A VMD script to define a custom color scale used for rendering crystalline particles (blue) and amorphous particles (yellow) at all time steps listed in the trajectory file.

##### vmd_scripts/setuserfield.vmd

A VMD script to set the values of 'user' field as the values in the 'vx' column of LAMMPS trajectory file, at every time step listed in the trajectory file.

### New folder created

After executing the code, an 'output' folder will be created under the 'run' folder, which contains three files:

. <br />
├── run <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── output <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── identified.lammpstrj <br />
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── identified_crystalline_atoms <br />
└&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── crystallinty

##### run/output/identified.lammpstrj

A new LAMMPS trajectory file, using a 'vx' column to indicate if a particle belongs to crystalline phase.

##### run/output/identified_crystalline_atoms

Each line lists the IDs of particles belongs to a crystalline chunk, for every time step listed in the trajectory file.
IDs are sorted according to particle connectivity.

##### run/output/crystallinty

Lists the crystallinty of polymer system at every time listed in the trajectory file.


# Example output (MPI mode)

With a cluster with a master node and two slave nodes, as listed in the 'host' file where each node has two physical cores, computation is launched by
```
cd run
mpiexec --hostfile ./hosts -np 6 python ../src/main.py mpi
```

### Screen output

> LAMMPS data file: data/cg-80blk-50chain.lammps
> LAMMPS trajectory file: data/cg-80blk-50chain-1atm-300K-40ns.lammpstrj
> 
> Time: 0.5 ns <br />
> atom_id |      p2  is_crys          proc
> ------- | -------  -------  ------------
>       0 | -0.0079    false  Master-rank0
>    2000 | -0.0321    false  Slave1-rank3
>    2750 | -0.0792    false  Slave2-rank4
>     750 | -0.0521    false  Master-rank1
>    3500 | -0.0912    false  Slave2-rank5
>    1500 |  0.0909    false  Slave1-rank2
>     250 |  0.0148    false  Master-rank0
>    2250 | -0.0451    false  Slave1-rank3
>    3000 |  0.0114    false  Slave2-rank4
>    1000 |  0.0562    false  Master-rank1
>    3750 |  0.1311    false  Slave2-rank5
>    1750 | -0.0259    false  Slave1-rank2
>     500 |  0.0109    false  Master-rank0
>    2500 |  0.0604    false  Slave1-rank3
>    3250 | -0.0114    false  Slave2-rank4
>    1250 |  0.0521    false  Master-rank1
>    3999 |  0.1746    false  Slave2-rank5
> Filtering crystalline atoms...
> Done
> Crystalline atoms:  0 / 4000
> Crystallinty:  0.0
> 
> Time: 25.0 ns
> atom_id |      p2  is_crys          proc
> ------- | -------  -------  ------------
>       0 |  0.0851    false  Master-rank0
>    2000 |  0.5543     true  Slave1-rank3
>    2750 |  0.2190    false  Slave2-rank4
>     750 |  0.1121    false  Master-rank1
>    1500 |  0.5630     true  Slave1-rank2
>    3500 |  0.8161     true  Slave2-rank5
>    2250 |  0.4310     true  Slave1-rank3
>     250 |  0.8452     true  Master-rank0
>    3000 |  0.7537     true  Slave2-rank4
>    1000 |  0.1774    false  Master-rank1
>    3750 |  0.6286     true  Slave2-rank5
>    1750 |  0.6259     true  Slave1-rank2
>    2500 | -0.0190    false  Slave1-rank3
>     500 |  0.5150     true  Master-rank0
>    3250 |  0.5011     true  Slave2-rank4
>    1250 | -0.0591    false  Master-rank1
>    3999 |  0.1601    false  Slave2-rank5
> Filtering crystalline atoms...
> Done
> Crystalline atoms:  1983 / 4000
> Crystallinty: 0.496


# Visualization

I use VMD to visualize the identified polymer system.
For this, I installed:

&nbsp;&nbsp;&nbsp;&nbsp;1. Tachyon <br />
&nbsp;&nbsp;&nbsp;&nbsp;2. VMD

Then open the VMD GUI interface with
```
cd run
vmd output/identified.lammpstrj 
```
is opened. Now modify some settings as stated below: 

'Display' -> check 'Orthographic' <br />
'Display' -> uncheck 'Depth Cueing' <br />
'Display' -> 'Axes' -> 'Off'

'Extensions' -> 'Tk Console', this opens a console, from where execute these commands:
```
  play ../vmd_scripts/color_scale.vmd
  play ../vmd_scripts/setuserfield.vmd
  pbc box -center unitcell -shiftcenter {0.168269 0.168269 0.168269}  -color gray
```

'Graphics' -> 'Representations...', this opens a new panel: <br />
&nbsp;&nbsp;&nbsp;&nbsp;Modify 'Coloring Method' to 'Trajectory' -> 'User' -> 'User' <br />
&nbsp;&nbsp;&nbsp;&nbsp;Modify 'Drawing Method' -> 'VDW' (Sphere Scale 0.6) <br />
&nbsp;&nbsp;&nbsp;&nbsp;Modify 'Material' -> 'AOShiny' <br />
&nbsp;&nbsp;&nbsp;&nbsp;Click  'Create Rep' <br />
&nbsp;&nbsp;&nbsp;&nbsp;Modify 'Drawing Method' -> 'DynamicBonds' (Distance Cutoff 2.7)

'Display' -> 'Display Settings' -> check 'Shadows' On <br />
'Display' -> 'Display Settings' -> Check 'Amb. Occl.' On

'Graphics' -> 'Colors...' -> Categories 'Display' -> Names 'Background' -> Colors '8 White'

'Extensions' -> 'Tk Console', execute following command to render an image file:
```
   render Tachyon vmdscene.dat tachyon -aasamples 24 -fullshade -res 1000 1000 %s -format PNG -o image.png
```
![Example rendering.](./run/image.png)

