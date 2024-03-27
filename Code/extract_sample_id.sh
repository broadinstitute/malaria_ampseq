#!/bin/bash

# Function to extract prefix before the second underscore
extract_prefix() {
    # Extract the part before the second underscore
    prefix=$(echo "$1" | cut -d'_' -f1-2)
    echo "$prefix"
}

# Check if a directory argument is provided
if [ $# -eq 1 ]; then
    directory="$1"

    # Ensure the directory exists
    if [ -d "$directory" ]; then
        # Process each file in the directory
        for file_path in "$directory"/*_L001_R1_001.fastq.gz; do
            # Extract the file name from the path
            file_name=$(basename "$file_path")
            
            # Call the function to extract the prefix
            result=$(extract_prefix "$file_name")
            
            # Print the result
            echo "$result"
        done
    else
        echo "Directory does not exist: $directory"
    fi
else
    echo "Usage: $0 <directory>"
fi
