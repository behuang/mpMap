module load cuda/5.0.35
R CMD INSTALL --configure-args='--with-cuda-home=$CUDA_HOME/cuda' . --preclean
