#!/bin/bash
#BSUB -R "rusage[mem=50GB]"
#BSUB -W 128:00

module load anaconda3/2019.07
source activate /icgc/dkfzlsdf/analysis/C010/scRNA_lung_mouse/cnmf_env

cd /icgc/dkfzlsdf/analysis/C010/scRNA_lung_mouse/cNMF/

python ./cnmf.py prepare --output-dir ./integrated_no_ciliated -c /icgc/dkfzlsdf/analysis/C010/scRNA_lung_mouse/WOT/data/integrated_all_tp.h5ad -k 8 9 10 11 --n-iter 100 --total-workers 10 --seed 14 --numgenes 2000

python ./cnmf.py factorize --output-dir ./integrated_no_ciliated --worker-index 0 
python ./cnmf.py factorize --output-dir ./integrated_no_ciliated --worker-index 1 
python ./cnmf.py factorize --output-dir ./integrated_no_ciliated --worker-index 2 
python ./cnmf.py factorize --output-dir ./integrated_no_ciliated --worker-index 3 
python ./cnmf.py factorize --output-dir ./integrated_no_ciliated --worker-index 4 
python ./cnmf.py factorize --output-dir ./integrated_no_ciliated --worker-index 5 
python ./cnmf.py factorize --output-dir ./integrated_no_ciliated --worker-index 6 
python ./cnmf.py factorize --output-dir ./integrated_no_ciliated --worker-index 7 
python ./cnmf.py factorize --output-dir ./integrated_no_ciliated --worker-index 8
python ./cnmf.py factorize --output-dir ./integrated_no_ciliated --worker-index 9 

# if tehre is a problem with the K selection plot
#The error is happening on the step where the Frobenius norm is being computed

#LINE 562: prediction_error = ((norm_counts.X.todense() - rf_pred_norm_counts)**2).sum().sum()
#Maybe try replacing it with something like this?

#LINE 562: prediction_error = ((np.array(norm_counts.X.todense()) - rf_pred_norm_counts)**2).sum().sum()
 
 
python ./cnmf.py combine --output-dir ./integrated_no_ciliated 

python ./cnmf.py k_selection_plot --output-dir ./integrated_no_ciliated 

python ./cnmf.py consensus --output-dir ./integrated_no_ciliated  --components 9 --local-density-threshold 0.01 --show-clustering