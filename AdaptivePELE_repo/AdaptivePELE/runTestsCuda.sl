#!/bin/bash
#SBATCH --job-name="AdaptiveTest_CUDA"
#SBATCH -D .
#SBATCH --output=AdaptiveTest_CUDA.out
#SBATCH --error=AdaptiveTest_CUDA.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=5
#SBATCH --time=00:20:00
#SBATCH --qos=debug
#SBATCH --constraint=k80
#SBATCH --gres gpu:2


srun python runAllTests.py --run MD_CUDA
