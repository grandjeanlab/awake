#!/bin/bash

#pip install git+https://github.com/brkraw/bruker 
#add to .bashrc --> PATH=$PATH:/groupshare/traaffneu/preclinimg/software/Bru2


#conda activate bruker
#cd /home/traaffneu/margal/awake_code/awake/assets/
# ./bids_convert.sh


# ---- Init varaibles ----
module load afni
root_dir='/project/4180000.36/to_convert'

cd $root_dir'/LabName/data/'                            #directory of data to be converted into bids format

output_dir=$root_dir'/LabName/export'                   #create the variable output_dir to be put in folder export  
mkdir -p $output_dir 

dataset_name='LabName'                                  #name dataset received 

ds_type=01                                              #database num, dataset num, scan num
ds_id=001                                              
id=0  

# --- Convert to BIDS ---- 

ls . | while read line                        
do                                            
cd $root_dir'/LabName/data/'$line             #go to directory of each line/folder

sub=$ds_type$ds_id'0'$id                      #subject ID
ses='1'                                       #session ID

output_sub_dir=$output_dir'/sub-'$sub'/ses-'$ses                                          
mkdir -p $output_sub_dir'/anat'                                                          
mkdir -p $output_sub_dir'/func'

anat_name='sub-'$sub'_ses-'$ses'_T2w.nii.gz'                                                   
func_name='sub-'$sub'_ses-'$ses'_run-1_bold.nii.gz'

#anat
Bru2 -a -z -o tmp anat/pdata/1/                                                                 #convert create tmp folder into anat/acqp directory -> convert
3dresample -inset tmp.nii.gz -prefix  $output_sub_dir'/anat/'$anat_name -orient LPI             #change scan orientation, save tp output directory 
rm tmp.nii.gz
# cp -u /anat/*nii.gz $output_sub_dir'/anat/tmp.nii.gz'                          #if no transformation needed


#func
Bru2 -z -a -o tmp func/pdata/1/                                                                 #convert data in folder 6 into BIDS format
3dresample -inset tmp.nii.gz -prefix  $output_sub_dir'/func/'$func_name -orient LPI             #change scan orientation, save into output directory 
rm tmp.nii.gz
# cp -u /func/*nii.gz $output_sub_dir'/func/tmp.nii.gz'                          #if no transformation needed


id=$((id + 1))                  #scan number increases by 1

cd ..                           #go to previous directory
done                            #end of the loop
cd ../export                    
tree                            #display organisation of data


