module load afni
cd /project/4180000.36/awake/bids/to_correct/
ls datanat/*/*/anat/*.nii.gz | while read line
do

echo $line
pixdim1=$(fslhd $line | grep pixdim1 | cut -f3)
pixdim1=`echo $pixdim1 / 10 | bc -l`
pixdim2=$(fslhd $line | grep pixdim2 | cut -f3)
pixdim2=`echo $pixdim2 / 10 | bc -l`
pixdim3=$(fslhd $line | grep pixdim3 | cut -f3)
pixdim3=`echo $pixdim3 / 10 | bc -l`

fslchpixdim $line $pixdim1 $pixdim2 $pixdim3
fslswapdim $line x -y -z $line
fslorient -deleteorient $line

3dresample -input $line -prefix tmp.nii.gz -orient ras
rm $line
mv tmp.nii.gz $line

done
