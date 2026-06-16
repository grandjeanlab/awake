module load afni


check_pixdim_and_continue() {
    local file="$1"

    # Run the command and capture output
    local result
    result=$(fslinfo "$file" | grep pixdim1 | grep 0.4)

    # If empty -> continue processing (return 0)
    # If not empty -> signal to skip (return 1)
    if [[ -z "$result" ]]; then
        return 0   # OK to continue
    else
        return 1   # Skip item
    fi
}


cd /project/4180000.36/awake/output_mouse/
fd -t f -e .nii.gz cleaned ./*/confound_correction_datasink/cleaned_timeseries/ | while read line
do

check_pixdim_and_continue "$line" || continue

echo $line

sbatch --mem=8G --time=00:5:00 << EOF
#!/bin/bash
module load afni

echo \$TMPDIR
3dresample -input $line -prefix \$TMPDIR/tmp.nii.gz -dxyz 0.4 0.4 0.4

rm $line
mv \$TMPDIR/tmp.nii.gz $line

EOF

done








