#!/bin/bash

# Loop through all .h5 files in the current directory and subdirectories
find "eon/runs/1e8m/gprd/final_gprd_wparam" -name "*.h5" -print0 | while IFS= read -r -d $'\0' h5file; do
    # Extract the directory
    h5dir=$(dirname "$h5file")

    # Determine if it's a singlet or doublet based on directory name
    if [[ "$h5dir" == *"doublets"* ]]; then
        prefix="d"
    elif [[ "$h5dir" == *"singlets"* ]]; then
        prefix="s"
    else
        echo "Skipping file $h5file: Could not determine singlet/doublet type."
        continue
    fi

    # Extract the number reliably regardless of directory depth
    number=$(echo "$h5dir" | grep -oE '/[0-9]+[^/]*$' | grep -oE '[0-9]+')

    if [[ -z "$number" ]]; then
        echo "Skipping file $h5file: Could not extract number from directory."
        continue
    fi

    output_name="${prefix}_${number}.traj"

    # Run the Python script
    python plot_gprdruns.py hdf5-to-traj "$h5file" --output_name "$output_name" --output_dir .

    echo "Processed: $h5file -> $output_name"

done

echo "Finished processing all .h5 files."
