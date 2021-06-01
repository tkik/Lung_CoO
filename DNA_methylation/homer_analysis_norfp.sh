module load anaconda3/2019.07
source activate /icgc/dkfzlsdf/analysis/C010/lung_mouse_cells/homer

cd /icgc/dkfzlsdf/analysis/C010/lung_mouse_cells/methylation_analysis/homer_analysis/homer/bin/
bed_dir=/icgc/dkfzlsdf/analysis/C010/lung_mouse_cells/methylation_analysis/homer_analysis/beds_new_groups_norfp
output_dir=/icgc/dkfzlsdf/analysis/C010/lung_mouse_cells/methylation_analysis/homer_analysis/homer_res_new_groups_norfp/
beds=`ls -l $bed_dir | awk '{print $9}'`

perl ../configureHomer.pl -install mm10
perl ../configureHomer.pl -install mouse

for bed in $beds
do
comp="${bed%.*}"
echo $comp
#mkdir $output_dir/$comp
/icgc/dkfzlsdf/analysis/C010/lung_mouse_cells/methylation_analysis/homer_analysis/homer/bin/findMotifsGenome.pl $bed_dir/$bed mm10  $output_dir/$comp  -len 8,10,12 -size 100 -S 8 -p 8 -cache 6921 -fdr 0
done

output_dir=/icgc/dkfzlsdf/analysis/C010/lung_mouse_cells/methylation_analysis/homer_analysis/

#/icgc/dkfzlsdf/analysis/C010/lung_mouse_cells/methylation_analysis/homer_analysis/homer/bin/findMotifs.pl /icgc/dkfzlsdf/analysis/C010/lung_mouse_cells/methylation_analysis/homer_analysis/MsCCSP_MsSPC_tumor_different_genes.txt mouse  $output_dir -start -400 -end 100 -len 8,10,12 -p 8
#-len 8,10,12 -size 100 -S 8 -p 8 -cache 6921 -fdr 0
