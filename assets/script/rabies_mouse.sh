#!/bin/sh bash

# author: Joanes Grandjean
# intial date: 05.05.2022
# last modified: 08.04.2026 (Jo)

# changelog
# 8.04.2026 
# add rabies_version to make rabies version selection more explicit
# 29.05.2024
# use the --inclusion_ids flag to run rabies one func at a time
# use the $TMPDIR environment variable to run rabies on /scratch and not on local project folder
# 03.08.2024
# modify it to run on the awake project. 


#define what root dir you want to use, where the bids folder is, where the tmp scripts will go, and where the output will go
root_dir="/project/4180000.36/awake/"
bids=$root_dir"bids/mouse/"
script_dir=$root_dir"/tmp_rabies_scripts_mouse"
output_dir=$root_dir"/output_mouse"
template_dir="/home/traaffneu/joagra/code/awake/assets/template/mouse"

train_fix=true
train_fix_dir=$root_dir"/train_fixrabies_mouse"
train_fix_model='awakeMouseV1'

#define what version of rabies you want to use. run `ls /opt/rabies/` to see what versions are on
rabies_version=0.6.4
rabies="/opt/rabies/${rabies_version}/rabies.sif"

#arguments for RABIES preprocessing, confound regression, analysis. see https://rabies.readthedocs.io/ for more info
prep_arg='--commonspace_resampling 0.2x0.2x0.2 --anatomical_resampling 0.2x0.2x0.2 --remove_dummy --oblique2card 3dWarp --despiking apply=true,for_reg_only=false,smoothing_sigma=None,threshold=4.0,gen_report=false,fps=10 --commonspace_reg masking=false,brain_extraction=false,template_registration=SyN,fast_commonspace=true --anat_template '${template_dir}'/template.nii.gz --brain_mask '${template_dir}'/mask.nii.gz --WM_mask '${template_dir}'/wm.nii.gz --CSF_mask '${template_dir}'/csf.nii.gz --vascular_mask '${template_dir}'/csf.nii.gz' 

conf_arg_gen=' --smoothing_filter 0.5 --highpass 0.01 --lowpass 0.1 --image_scaling voxelwise_standardization --read_datasink --TR '

conf_arg_wmcsf1=' --nuisance_regressors mot_6 WM_signal CSF_signal --frame_censoring FD_censoring=true,FD_threshold=0.1,DVARS_censoring=true,minimum_timepoint=3'
conf_arg_wmcsf2=' --nuisance_regressors mot_6 WM_signal CSF_signal --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=true,minimum_timepoint=3'
conf_arg_wmcsf3=' --nuisance_regressors mot_6 WM_signal CSF_signal --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=false,minimum_timepoint=3'

conf_arg_gsr1=' --nuisance_regressors mot_6 global_signal --frame_censoring FD_censoring=true,FD_threshold=0.1,DVARS_censoring=true,minimum_timepoint=3'
conf_arg_gsr2=' --nuisance_regressors mot_6 global_signal --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=true,minimum_timepoint=3'
conf_arg_gsr3=' --nuisance_regressors mot_6 global_signal --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=false,minimum_timepoint=3'

conf_arg_aCompCor1=' --nuisance_regressors mot_6 aCompCor_percent --frame_censoring FD_censoring=true,FD_threshold=0.1,DVARS_censoring=true,minimum_timepoint=3'
conf_arg_aCompCor2=' --nuisance_regressors mot_6 aCompCor_percent --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=true,minimum_timepoint=3'
conf_arg_aCompCor3=' --nuisance_regressors mot_6 aCompCor_percent --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=false,minimum_timepoint=3'


analysis_arg="-DR_template ~/code/awake/assets/template/mouse/ica.nii.gz -extract_threshold 7.5z -seed '{\"apply\": true, \"seeds\": [\"~/code/awake/assets/template/mouse/s1_r.nii.gz\", \"~/code/awake/assets/template/mouse/s1_l.nii.gz\",\"~/code/awake/assets/template/mouse/aca_r.nii.gz\"]}'"

#make the script directory. this is where your runnable rabies script per func scan will be run. 
mkdir -p $script_dir

mkdir -p $output_dir/wmcsf1
mkdir -p $output_dir/wmcsf2
mkdir -p $output_dir/wmcsf3
mkdir -p $output_dir/gsr1
mkdir -p $output_dir/gsr2
mkdir -p $output_dir/gsr3
mkdir -p $output_dir/aCompCor1
mkdir -p $output_dir/aCompCor2
mkdir -p $output_dir/aCompCor3
mkdir -p $output_dir/fixrabies


cd $script_dir

#this is the main loop. by default, it will loop over every func scan that you have in your bids directory and make a separate script for it. 
find $bids -name *_bold.nii.gz | while read line
do

#need to find the corresponding json to find tr, extract it from json, and div by 1000 to get val in sec.
json=$(echo $line | sed "s/.nii.gz/.json/g")
tr=$(grep "RepetitionTime" $json | sed -r 's/"RepetitionTime"://g' | sed "s/[^[:digit:]]//g")
trcor=`echo $tr / 1000 | bc -l`


#check for anatomical scans
anat_folder=$(dirname $line)
anat_scan=$(dirname $anat_folder)'/anat/*.nii.gz'
bold_only=''
if [ ! -f $anat_scan ]; then
    bold_only=' --bold_only'
fi



#edit the func file name and path for rabies
##replace the full path to the bids directory with a relative path for rabies
func_file=$line


##set the name of the script file that will be created. 
func_base=$(basename $func_file)
func_noext="$(remove_ext $func_base)"
script_file=$script_dir/$func_noext'.sh'

echo "now doing subject "$func_noext

#initialize the script with a bang and slurm header. you can edit the time and mem options if you think you need more or less resources. 
echo '#!/bin/bash' > $script_file
echo "#SBATCH --job-name="$func_noext >> $script_file
echo "#SBATCH --nodes=1" >> $script_file
echo "#SBATCH --time=12:00:00" >> $script_file
echo "#SBATCH --mail-type=FAIL" >> $script_file
echo "#SBATCH --partition=batch" >> $script_file
echo "#SBATCH --mem=24GB" >> $script_file

#create temporary folders in scratch folder so you don't clutter your project folder
echo " " >> $script_file
echo "" >> $script_file
echo " " >> $script_file
echo "#### init varibles and make tmp directories ####" >> $script_file
echo " " >> $script_file

echo "preprocess=$""TMPDIR/preprocess" >> $script_file

echo "confound_wmcsf1=$""TMPDIR/confound_wmcsf1" >> $script_file
echo "analysis_wmcsf1=$""TMPDIR/analysis_wmcsf1" >> $script_file
echo "confound_wmcsf2=$""TMPDIR/confound_wmcsf2" >> $script_file
echo "analysis_wmcsf2=$""TMPDIR/analysis_wmcsf2" >> $script_file
echo "confound_wmcsf3=$""TMPDIR/confound_wmcsf3" >> $script_file
echo "analysis_wmcsf3=$""TMPDIR/analysis_wmcsf3" >> $script_file

echo "confound_gsr1=$""TMPDIR/confound_gsr1" >> $script_file
echo "analysis_gsr1=$""TMPDIR/analysis_gsr1" >> $script_file
echo "confound_gsr2=$""TMPDIR/confound_gsr2" >> $script_file
echo "analysis_gsr2=$""TMPDIR/analysis_gsr2" >> $script_file
echo "confound_gsr3=$""TMPDIR/confound_gsr3" >> $script_file
echo "analysis_gsr3=$""TMPDIR/analysis_gsr3" >> $script_file

echo "confound_aCompCor1=$""TMPDIR/confound_aCompCor1" >> $script_file
echo "analysis_aCompCor1=$""TMPDIR/analysis_aCompCor1" >> $script_file
echo "confound_aCompCor2=$""TMPDIR/confound_aCompCor2" >> $script_file
echo "analysis_aCompCor2=$""TMPDIR/analysis_aCompCor2" >> $script_file
echo "confound_aCompCor3=$""TMPDIR/confound_aCompCor3" >> $script_file
echo "analysis_aCompCor3=$""TMPDIR/analysis_aCompCor3" >> $script_file

echo "fix_dir=$""TMPDIR/analysis_fix" >> $script_file



echo "mkdir -p $""preprocess" >> $script_file


echo "mkdir -p $""confound_wmcsf1" >> $script_file
echo "mkdir -p $""analysis_wmcsf1" >> $script_file 
echo "mkdir -p $""confound_wmcsf2" >> $script_file
echo "mkdir -p $""analysis_wmcsf2" >> $script_file 
echo "mkdir -p $""confound_wmcsf3" >> $script_file
echo "mkdir -p $""analysis_wmcsf3" >> $script_file 

echo "mkdir -p $""confound_gsr1" >> $script_file
echo "mkdir -p $""analysis_gsr1" >> $script_file 
echo "mkdir -p $""confound_gsr2" >> $script_file
echo "mkdir -p $""analysis_gsr2" >> $script_file 
echo "mkdir -p $""confound_gsr3" >> $script_file
echo "mkdir -p $""analysis_gsr3" >> $script_file 

echo "mkdir -p $""confound_aCompCor1" >> $script_file
echo "mkdir -p $""analysis_aCompCor1" >> $script_file 
echo "mkdir -p $""confound_aCompCor2" >> $script_file
echo "mkdir -p $""analysis_aCompCor2" >> $script_file 
echo "mkdir -p $""confound_aCompCor3" >> $script_file
echo "mkdir -p $""analysis_aCompCor3" >> $script_file 

echo "mkdir -p $""fix_dir" >> $script_file 

echo " " >> $script_file
echo "#### run RABIES preprocess ####" >> $script_file
echo " " >> $script_file


#run the preprocessing step of rabies
echo "apptainer run "${rabies}" --inclusion_ids "${func_file}" -p Linear preprocess "${bids}" $""{preprocess} "${bold_only}${prep_arg} >> $script_file 

if $train_fix; 
then
echo "#### run RABIES FIX training ####" >> $script_file
echo "source ~/code/FixRABIES/.venv/bin/activate" >> $script_file 
#run the FIX training step of rabies
echo "fixrabies prepare -fix_dir "$train_fix_dir" -preprocess_dir $""{preprocess} -species mouseDSURQE -FD_threshold 0.2 -tr "${trcor} >> $script_file 

echo "cp -r $""{preprocess}/preprocess_QC_report "$train_fix_dir >> $script_file 
echo "cp -r $""{preprocess}/motion_datasink "$train_fix_dir >> $script_file 
echo "cp -r $""{preprocess}/bold_datasink/tSNR_map_preprocess "$train_fix_dir >> $script_file 

#the real workflow
else


#copy the QC report, motion, and tSNR maps
echo "cp -r $""{preprocess}/preprocess_QC_report "$output_dir >> $script_file 
echo "cp -r $""{preprocess}/motion_datasink "$output_dir >> $script_file 
echo "cp -r $""{preprocess}/bold_datasink/tSNR_map_preprocess "$output_dir >> $script_file 

echo "source ~/code/FixRABIES/.venv/bin/activate" >> $script_file 

echo " " >> $script_file
echo "#### run RABIES confound/analysis for white matter / csf regression ####" >> $script_file
echo " " >> $script_file

#run the confound correction step of rabies
echo "apptainer run "${rabies}" --inclusion_ids "${func_file}" -p Linear confound_correction $""{preprocess} $""{confound_wmcsf1} "${conf_arg_wmcsf1}${conf_arg_gen}${trcor} >> $script_file 
#run the analysis step of fixrabies
echo "fixrabies import -fix_dir $""{analysis_wmcsf1} -confound_dir $""{confound_wmcsf1} -preprocess_dir $""{preprocess} -species mouseDSURQE" >> $script_file
echo "fixrabies analyze -fix_dir $""{analysis_wmcsf1} "${analysis_arg} >> $script_file 
#copy the analysis outputs and the data diagnosis to the output directory
echo "cp -r $""confound_wmcsf1/confound_correction_datasink/ "$output_dir"/wmcsf1" >> $script_file 
echo "cp -r $""analysis_wmcsf1/results "$output_dir"/"$conf_name_final >> $script_file 

#run the confound correction step of rabies
echo "apptainer run "${rabies}" --inclusion_ids "${func_file}" -p Linear confound_correction $""{preprocess} $""{confound_wmcsf2} "${conf_arg_wmcsf2}${conf_arg_gen}${trcor} >> $script_file 
#run the analysis step of fixrabies
echo "fixrabies import -fix_dir $""{analysis_wmcsf2} -confound_dir $""{confound_wmcsf2} -preprocess_dir $""{preprocess} -species mouseDSURQE" >> $script_file
echo "fixrabies analyze -fix_dir $""{analysis_wmcsf2} "${analysis_arg} >> $script_file 
#copy the analysis outputs and the data diagnosis to the output directory
echo "cp -r $""confound_wmcsf2/confound_correction_datasink/ "$output_dir"/wmcsf2" >> $script_file 
echo "cp -r $""analysis_wmcsf2/results "$output_dir"/"$conf_name_final >> $script_file 


#run the confound correction step of rabies
echo "apptainer run "${rabies}" --inclusion_ids "${func_file}" -p Linear confound_correction $""{preprocess} $""{confound_wmcsf3} "${conf_arg_wmcsf3}${conf_arg_gen}${trcor}  >> $script_file 
#run the analysis step of fixrabies
echo "fixrabies import -fix_dir $""{analysis_wmcsf3} -confound_dir $""{confound_wmcsf3} -preprocess_dir $""{preprocess} -species mouseDSURQE" >> $script_file
echo "fixrabies analyze -fix_dir $""{analysis_wmcsf3} "${analysis_arg} >> $script_file 
#copy the analysis outputs and the data diagnosis to the output directory
echo "cp -r $""confound_wmcsf3/confound_correction_datasink/ "$output_dir"/wmcsf3" >> $script_file 
echo "cp -r $""analysis_wmcsf3/results "$output_dir"/"$conf_name_final >> $script_file 



echo "#### run RABIES confound/analysis for global signal regression ####" >> $script_file
echo " " >> $script_file

#run the confound correction step of rabies
echo "apptainer run "${rabies}" --inclusion_ids "${func_file}" -p Linear confound_correction $""{preprocess} $""{confound_gsr1} "${conf_arg_gsr1}${conf_arg_gen}${trcor}  >> $script_file 
#run the analysis step of fixrabies
echo "fixrabies import -fix_dir $""{analysis_gsr1} -confound_dir $""{confound_gsr1} -preprocess_dir $""{preprocess} -species mouseDSURQE" >> $script_file
echo "fixrabies analyze -fix_dir $""{analysis_gsr1} "${analysis_arg} >> $script_file 
#copy the analysis outputs and the data diagnosis to the output directory
echo "cp -r $""confound_gsr1/confound_correction_datasink/ "$output_dir"/gsr1" >> $script_file 
echo "cp -r $""analysis_gsr1/results "$output_dir"/"$conf_name_final >> $script_file 


#run the confound correction step of rabies
echo "apptainer run "${rabies}" --inclusion_ids "${func_file}" -p Linear confound_correction $""{preprocess} $""{confound_gsr2} "${conf_arg_gsr2}${conf_arg_gen}${trcor}  >> $script_file 
#run the analysis step of fixrabies
echo "fixrabies import -fix_dir $""{analysis_gsr2} -confound_dir $""{confound_gsr2} -preprocess_dir $""{preprocess} -species mouseDSURQE" >> $script_file
echo "fixrabies analyze -fix_dir $""{analysis_gsr2} "${analysis_arg} >> $script_file 
#copy the analysis outputs and the data diagnosis to the output directory
echo "cp -r $""confound_gsr2/confound_correction_datasink/ "$output_dir"/gsr2" >> $script_file 
echo "cp -r $""analysis_gsr2/results "$output_dir"/"$conf_name_final >> $script_file 


#run the confound correction step of rabies
echo "apptainer run "${rabies}" --inclusion_ids "${func_file}" -p Linear confound_correction $""{preprocess} $""{confound_gsr3} "${conf_arg_gsr3}${conf_arg_gen}${trcor}  >> $script_file 
#run the analysis step of fixrabies
echo "fixrabies import -fix_dir $""{analysis_gsr3} -confound_dir $""{confound_gsr3} -preprocess_dir $""{preprocess} -species mouseDSURQE" >> $script_file
echo "fixrabies analyze -fix_dir $""{analysis_gsr3} "${analysis_arg} >> $script_file 
#copy the analysis outputs and the data diagnosis to the output directory
echo "cp -r $""confound_gsr3/confound_correction_datasink/ "$output_dir"/gsr3" >> $script_file 
echo "cp -r $""analysis_gsr3/results "$output_dir"/"$conf_name_final >> $script_file 


echo " " >> $script_file
echo "#### run RABIES confound/analysis for white matter / csf regression ####" >> $script_file
echo " " >> $script_file

#run the confound correction step of rabies
echo "apptainer run "${rabies}" --inclusion_ids "${func_file}" -p Linear confound_correction $""{preprocess} $""{confound_aCompCor1} "${conf_arg_aCompCor1}${conf_arg_gen}${trcor}  >> $script_file 
#run the analysis step of fixrabies
echo "fixrabies import -fix_dir $""{analysis_aCompCor1} -confound_dir $""{confound_aCompCor1} -preprocess_dir $""{preprocess} -species mouseDSURQE" >> $script_file
echo "fixrabies analyze -fix_dir $""{analysis_aCompCor1} "${analysis_arg} >> $script_file 
#copy the analysis outputs and the data diagnosis to the output directory
echo "cp -r $""confound_aCompCor1/confound_correction_datasink/ "$output_dir"/aCompCor1" >> $script_file 
echo "cp -r $""analysis_aCompCor1/results "$output_dir"/"$conf_name_final >> $script_file 

#run the confound correction step of rabies
echo "apptainer run "${rabies}" --inclusion_ids "${func_file}" -p Linear confound_correction $""{preprocess} $""{confound_aCompCor2} "${conf_arg_aCompCor2}${conf_arg_gen}${trcor}  >> $script_file 
#run the analysis step of fixrabies
echo "fixrabies import -fix_dir $""{analysis_aCompCor2} -confound_dir $""{confound_aCompCor2} -preprocess_dir $""{preprocess} -species mouseDSURQE" >> $script_file
echo "fixrabies analyze -fix_dir $""{analysis_aCompCor2} "${analysis_arg} >> $script_file 
#copy the analysis outputs and the data diagnosis to the output directory
echo "cp -r $""confound_aCompCor2/confound_correction_datasink/ "$output_dir"/aCompCor2" >> $script_file 
echo "cp -r $""analysis_aCompCor2/results "$output_dir"/"$conf_name_final >> $script_file 

#run the confound correction step of rabies
echo "apptainer run "${rabies}" --inclusion_ids "${func_file}" -p Linear confound_correction $""{preprocess} $""{confound_aCompCor3} "${conf_arg_aCompCor3}${conf_arg_gen}${trcor}  >> $script_file 
#run the analysis step of fixrabies
echo "fixrabies import -fix_dir $""{analysis_aCompCor3} -confound_dir $""{confound_aCompCor3} -preprocess_dir $""{preprocess} -species mouseDSURQE" >> $script_file
echo "fixrabies analyze -fix_dir $""{analysis_aCompCor3} "${analysis_arg} >> $script_file 
#copy the analysis outputs and the data diagnosis to the output directory
echo "cp -r $""confound_aCompCor3/confound_correction_datasink/ "$output_dir"/aCompCor3" >> $script_file 
echo "cp -r $""analysis_aCompCor3/results "$output_dir"/"$conf_name_final >> $script_file 

echo " " >> $script_file
echo "#### clean up####" >> $script_file
echo " " >> $script_file

#clean up scratch
echo "rm -rf $""TMPDIR/*" >> $script_file 

fi

#uncomment one of the following if you want to run the scripts automatically (do so if you are confident it will work)
##this is if you are using the new slurm system
#sbatch $script_file

##this is if you are using the old pbs system
#qsub -l 'nodes=1,mem=16gb,walltime=12:00:00' $script_file 

#end of the loop
done

