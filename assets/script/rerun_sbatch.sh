#!/bin/env bash

tmp_scripts_dir='/project/4180000.36/awake/tmp_rabies_scripts_rat'
complete_scripts_dir='/project/4180000.36/awake/complete_rabies_scripts_rat'
output_dir='/project/4180000.36/awake/output_rat'

cd $tmp_scripts_dir

ls *.sh | while read line
do

#line='sub-0100200_ses-1_run-1_task-rest_bold.sh' 

#keep the first 3 sections of line separated by underscores
subject=$(echo $line | cut -d'_' -f1)
session=$(echo $line | cut -d'_' -f2)
run=$(echo $line | cut -d'_' -f3)

# build path to anat reg output

anat_reg_output=$(find ${output_dir}/preprocess_QC_report/commonspace_reg_wf.Native2Atlas -type f | grep ${subject}_${session})
#build path to the func reg output
func_reg_output=$(find ${output_dir}/preprocess_QC_report/EPI2Anat/ -type f | grep ${subject}_${session}_${run})

#build path to aCompCor seed output
acompcor_seed_output=$(find ${output_dir}/aCompCor3/analysis_datasink/seed_correlation_maps/_split_name_${subject}_${session}_${run}_task-rest_bold/_seed_name_vpm_r -type f | grep combined_cleaned_vpm_r_corr_map.nii.gz)

if [ -z $anat_reg_output ]; then 
echo "${subject}_${session}_${run} incomplete"
continue
fi
if [ -z $func_reg_output ]; then 
echo "${subject}_${session}_${run} incomplete"
continue
fi
if [ -z $acompcor_seed_output ]; then
echo "${subject}_${session}_${run} incomplete"
continue 
fi


#check if all files exist, if so, move the script to the complete directory, else run `sbatch` on the script
if [ -f $anat_reg_output ] && [ -f $func_reg_output ] && [ -f $acompcor_seed_output ]; then
    echo "${subject}_${session}_${run} complete"
    mv $line $complete_scripts_dir

else
    echo "${subject}_${session}_${run} incomplete"
fi

done

 
mkdir -p $output_dir/wmcsf1
mkdir -p $output_dir/wmcsf2
mkdir -p $output_dir/wmcsf3
mkdir -p $output_dir/gsr1
mkdir -p $output_dir/gsr2
mkdir -p $output_dir/gsr3
mkdir -p $output_dir/aCompCor1
mkdir -p $output_dir/aCompCor2
mkdir -p $output_dir/aCompCor3




#find . -type f -exec sed -i 's/3dWarp /3dWarp --bold_inho_cor method=N4_reg,otsu_thresh=2,multiotsu=false --anat_inho_cor method=N4_reg,otsu_thresh=2,multiotsu=false /g' {} +
#find . -type f -exec sed -i 's/--bold_inho_cor /--anat_robust_inho_cor apply=true,masking=false,brain_extraction=false,template_registration=SyN --bold_robust_inho_cor apply=true,masking=false,brain_extraction=false,template_registration=SyN --bold_inho_cor /g' {} +

#find . -type f -exec sed -i 's/mem=24GB/mem=64GB/g' {} +

#find . -type f -exec sed -i 's/--detect_dummy //g' {} +

#find . -type f -exec sed -i 's/fast_commonspace=true /fast_commonspace=true --bold_only --bold_inho_cor method=Rigid,otsu_thresh=2,multiotsu=false --bold_robust_inho_cor apply=true,masking=true,brain_extraction=true,template_registration=Affine/g' {} +

#find . -type f -exec sed -i 's/--bold_inho_cor method=N4_reg,otsu_thresh=2,multiotsu=false --anat_inho_cor method=N4_reg,otsu_thresh=2,multiotsu=false //g' {} +

#find . -type f -exec sed -i 's/3dWarp /3dWarp --bold_inho_cor method=N4_reg,otsu_thresh=2,multiotsu=false /g' {} +


#find . -type f -exec sed -i 's/mem=64GB/mem=128GB/g' {} +

#find . -type f -exec sed -i 's/mem=64GB/mem=128GB/g' {} +


#find . -type f -exec sed -i 's/--detect_dummy/--detect_dummy --bold_autobox --anat_autobox/g' {} +

#find . -type f -exec sed -i 's/--bold_inho_cor method=N4_reg,otsu_thresh=2,multiotsu=false --anat_inho_cor method=N4_reg,otsu_thresh=2,multiotsu=false //g' {} +


#find . -type f -exec sed -i 's/--detect_dummy //g' {} +
#find . -type f -exec sed -i 's/--detect_dummy /--detect_dummy --oblique2card 3dWarp /g' {} +
#find . -type f -exec sed -i 's/--oblique2card 3dWarp //g' {} +

#find . -type f -exec sed -i 's/--time=12:00:00/--time=24:00:00/g' {} +
#find . -type f -exec sed -i 's/--bold_autobox --anat_autobox //g' {} +

#find . -type f -exec sed -i 's/MultiProc/Linear/g' {} +



