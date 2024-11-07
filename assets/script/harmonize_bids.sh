

cd /project/4180000.36/AwakeRodent

# correct nifti to nifti_gz

find BIDS_input-RABIES/ -name *.nii | while read line
do
fslchfiletype NIFTI_GZ $line
done


# reset qform and sform matrices. 
module load afni
find BIDS_input-RABIES/ -name *.nii.gz | while read line
do
echo $line
fslorient -deleteorient $line
3dresample -input $line -prefix tmp.nii.gz -orient ras
rm $line
mv tmp.nii.gz $line
done



## script to correct orientations when it is hard-coded in the file

module load afni
cd /project/4180000.36/AwakeRodent/BIDS_input-RABIES/to_correct/rat/
find . -name *.nii.gz | while read line
do
echo $line
fslorient -deleteorient $line
fslswapdim $line x -z y tmp.nii.gz
3dresample -input tmp.nii.gz -prefix tmp2.nii.gz -orient ras
rm $line
rm tmp.nii.gz
mv tmp2.nii.gz $line
done






