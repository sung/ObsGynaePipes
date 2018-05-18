#!/bin/bash
#!
#! Example SLURM job script for CSD3 (Skylake CPUs, OPA)
#! Last updated: 26/Feb/2018
#! Modified by Sung Gong <ssg29>

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J SLX-ILM-Plasma2017.Homo_sapiens.v1.ILM-PLS.PL-35_S34.XXX.MY_LANE.MY_CHUNK.per.sample.salmon
#! Which project should be charged:
#SBATCH -A CHARNOCK-JONES-SL2-CPU 
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each
#SBATCH --ntasks=32
#! How much wallclock time will be required?
#SBATCH --time=0:30:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! long job (>7days)
##SBATCH --partition=long

#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":      
#SBATCH -p skylake

#SBATCH --error=SLX-ILM-Plasma2017.Homo_sapiens.v1.ILM-PLS.PL-35_S34.XXX.MY_LANE.MY_CHUNK.per.sample.salmon.log
#SBATCH --output=SLX-ILM-Plasma2017.Homo_sapiens.v1.ILM-PLS.PL-35_S34.XXX.MY_LANE.MY_CHUNK.per.sample.salmon.log

#! sbatch directives end here (put any additional directives above this line)

#! Notes:
#! Charging is determined by node number*walltime. Allocation is in entire nodes.
#! The --ntasks value refers to the number of tasks to be launched by SLURM only. This
#! usually equates to the number of MPI tasks launched. Reduce this from nodes*16 if
#! demanded by memory requirements, or if OMP_NUM_THREADS>1.

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
#! Modify the settings below to specify the application's environment, location 
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment            

#! Insert additional module load commands after this line if needed:
source $HOME/lib/modules.sh # defines necessary software to load 

#! Full path to application executable: 
application="$BIN_TOP/script/SLX-ILM-Plasma2017.Homo_sapiens.v1.PL-35_S34.per.sample.salmon.XXX.sh"
if [ ! -f $application ];then
	echo "$application not found"
	exit
fi

#! Run options for the application:
options=

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$[${numnodes}*${mpi_tasks_per_node}]

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.


#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
CMD="bash $application"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! (OMP_NUM_THREADS threads will be created):
#CMD="$application $options"

#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"


###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to `pwd`.\n"

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        export NODEFILE=`generate_pbs_nodefile`
        cat $NODEFILE | uniq > machine.file.$JOBID
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi

echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD 
