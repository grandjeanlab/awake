

cd /project/4180000.36/awake/complete_output_mouse/
dr_template='/home/traaffneu/joagra/code/awake/assets/template/mouse/ica.nii.gz'
mask='/home/traaffneu/joagra/code/awake/assets/template/mouse/mask.nii.gz'


fd -a -t f -e .nii.gz cleaned_DR_maps ./*/commonspace_analysis_datasink/dual_regression_nii/ | while read line
do

output_dir=$(dirname $line)
dice_file="$outfile_dir/dice"


# If dice file does NOT exist → continue loop
if [[ ! -f "$dice_file" ]]; then
 echo "Dice file missing for $line → processing..."


sbatch --mem=1G --time=00:10:00 << EOF
#!/bin/bash
module load afni
cd \$TMPDIR

output_dir=$(dirname $line)

3dresample -input $dr_template -prefix rs_dr.nii.gz -master $line
3dresample -input $mask -prefix rs_mask.nii.gz -master $line

fslsplit rs_dr.nii.gz IC -t
fslsplit $line DR -t


ls IC* | while read IC
do
echo \$IC

DR="\${IC/IC/DR}"

echo \$DR

fslmaths \$DR -mul rs_mask.nii.gz -thrp 90 -bin bin_stage2.nii.gz
fslmaths \$IC -mul rs_mask.nii.gz -thrp 90 -bin bin_IC.nii.gz
fslmaths bin_IC.nii.gz -mul bin_stage2.nii.gz intersection.nii.gz

# Count voxels (sum of values; assumes binary masks 0/1)
fslstats bin_stage2.nii.gz -V > statA
fslstats bin_IC.nii.gz -V > statB
fslstats intersection.nii.gz -V > statI

read volA _ < statA
read volB _ < statB
read volI _ < statI

# Compute Dice coefficient
echo "scale=6; (2 * \$volI) / (\$volA + \$volB)" | bc >> dice

cp dice \$output_dir 
rm bin_stage2.nii.gz
rm bin_IC.nii.gz
rm intersection.nii.gz
rm statA
rm statB
rm statI
done

rm -rf \$TMPDIR

EOF

sleep 2

else
 echo "Dice file exists for $line → skipping"
 continue
fi


done






