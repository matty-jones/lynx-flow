module load bwpy
source /u/eot/anderso4/projects/rhaco/rhaco-virtualenv-rhaas/bin/activate 
export PYTHONPATH="/u/eot/anderso4/software/build-hoomd-on-blue-waters/hoomd_blue/build"
echo "Running generate"
echo $1
aprun -n 1 python operations.py generate $1
#bwpy-environ python operations.py generate $1
