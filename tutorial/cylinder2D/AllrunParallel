#PBS -l walltime=240:30:00,nodes=1:ppn=12

rm -rf log
#decomposePar | tee -a log

mpirun -np 12 -machinefile $PBS_NODEFILE iceFoam -parallel >log 2>&1




