    #!/bin/bash
    
    #qsub -l 'procs=1,mem=24gb,walltime=12:00:00' -I 

    #cd /home/traaffneu/margal/awake_code/awake/assets/
    # ./fslorient.sh

    base_dir='/project/4180000.36/to_convert/Online/03_Mice/Gutierrez-Barragan_Gozzi/BIDS_fslorient/'
    # Loop through subdirectories
    for sub_dir in "$base_dir"sub-*/ses-*/func/; do
        cd "$sub_dir" || exit

        # Get the T2w file in the current subdirectory
        for t2w_file in *nii; do 
        if [ -f "$t2w_file" ]; then
            output_dir="$sub_dir"
            mkdir -p "$output_dir"

            echo "Processing: $sub_dir$t2w_file"

            # Copy the T2w file into the output directory and rename it accordingly
            cp "$sub_dir$t2w_file" "$output_dir/tmp1_$t2w_file"                                 
            rm "$sub_dir$t2w_file"

            # Apply orientation correction
            #fslorient -deleteorient imagename "$output_dir/tmp1_$t2w_file"   
            fslorient -setqformcode 1 "$output_dir/tmp1_$t2w_file"                            # Reset the qform code  
            fslorient -setsform -1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1 "$output_dir/tmp1_$t2w_file"   # Correct the labelled of the orientation (Q form), the -1 changes the Posterior-Anterior to Anterior-Posterior
            fslorient -setqform -1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1 "$output_dir/tmp1_$t2w_file"
            fslreorient2std "$output_dir/tmp1_$t2w_file" "$output_dir/$t2w_file"
            #fslswapdim "$output_dir/tpm2_$t2w_file" -x y z "$output_dir/$t2w_file"        

        
            # Optionally remove temporary files if needed
            rm "$output_dir/tmp1_$t2w_file"
            #rm "$output_dir/tmp2_$t2w_file"
            
        echo "Done!"
        else
            echo "No T2w file found in: $sub_dir"
        fi
    done
done



    # base_dir='/project/4180000.36/to_convert/Mice/LTN/export_1005_fslorient/'

    # # Loop through subdirectories
    # for sub_dir in "$base_dir"sub-*/ses-*/func/; do
    #     cd "$sub_dir" || exit

    #     # Get the T2w file in the current subdirectory
    #     t2w_file=$(ls -1 *.nii 2>/dev/null | head -n 1)
        
    #     # Check if T2w file exists in the current subdirectory
    #     if [ -n "$t2w_file" ]; then
    #         output_dir="$sub_dir"
    #         mkdir -p "$output_dir"

    #         echo "Processing: $sub_dir$t2w_file"

    #         # Copy the T2w file into the output directory and rename it accordingly
    #         cp "$sub_dir$t2w_file" "$output_dir/tmp1_$t2w_file"                                 
    #         #rm "$sub_dir$t2w_file"

    #         # Apply orientation correction
    #         fslorient -setsform 1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 -1 "$output_dir/tmp1_$t2w_file"   # correct the labelled of the orientation (Q form), the -1 changes the Posterior-Anterior to Anterior-Posterior
    #         fslorient -setqform 1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 -1 "$output_dir/tmp1_$t2w_file"
    #         fslreorient2std "$output_dir/tmp1_$t2w_file" "$output_dir/$t2w_file"
    #         #fslswapdim "$output_dir/tpm2_$t2w_file" -x y z "$output_dir/tpm3_$t2w_file"        

        
    #         # Optionally remove temporary files if needed
    #         rm "$output_dir/tmp1_$t2w_file"
    #         #rm "$output_dir/tmp2_$t2w_file"
            
    #         echo "Done!"
    #     else
    #         echo "No T2w file found in: $sub_dir"
    #     fi
    # done

