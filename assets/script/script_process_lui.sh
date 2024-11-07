
#unzip all the data files, only keep raw data
unzip alldatav1.1_part2.zip
rm -rf rat*/rfmri*
unzip alldatav1.1_part3.zip
rm -rf rat*/rfmri*
unzip alldatav1.1_part4.zip
rm -rf rat*/rfmri*
unzip alldatav1.1_part5.zip
rm -rf rat*/rfmri*
unzip alldatav1.1_part6.zip
rm -rf rat*/rfmri*

cd /project/4180000.36/AwakeRodent/to_convert/Online/04_Rats/YikangLiu/raw/data
output='/project/4180000.36/AwakeRodent/to_convert/Online/04_Rats/YikangLiu/raw/output'

ls | while read line
do
cd $line
echo $line
scan=$(find ./ -name *.nii -size +2M -size -10M)
cp $scan $output/$line'.nii'
cd ..
done


module load afni
cd /project/4180000.36/AwakeRodent/to_convert/Online/04_Rats/YikangLiu/raw/output/
ls | while read line
do
echo $line
fslchfiletype NIFTI_GZ $line
fslorient -deleteorient $line'.gz'
fslswapdim $line'.gz' x z -y tmp.nii.gz
3dresample -input tmp.nii.gz -prefix tmp2.nii.gz -orient ras
rm tmp.nii.gz
rm $line'.gz'
mv tmp2.nii.gz $line'.gz'
done


bids='/project/4180000.36/AwakeRodent/to_convert/Online/04_Rats/YikangLiu/raw/bids'
output='/project/4180000.36/AwakeRodent/to_convert/Online/04_Rats/YikangLiu/raw/output'
cd /project/4180000.36/AwakeRodent/to_convert/Online/04_Rats/YikangLiu/raw/
ls $output | while read line
do
echo $line
scan=$(remove_ext $line)
sub=$(grep $scan scan | cut -f1)
output_dir=$bids'/sub-0'$sub'/ses-1/anat/'
mkdir -p $output_dir
output_file=$output_dir'sub-0'$sub'_ses-1_T2w.nii.gz'
cp $output/$line $output_file

done
