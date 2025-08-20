#!/bin/bash

# Max file size in MB
MAX_SIZE=10

# Create .gitignore if it doesn't exist
touch .gitignore

# Find all files larger than MAX_SIZE MB
find . -type f -size +"${MAX_SIZE}M" | while read -r file; do
    # Remove leading ./ from find output
    file=${file#./}
    # Check if file is already in .gitignore
    if ! grep -Fxq "$file" .gitignore; then
        echo "$file" >> .gitignore
        echo "Added $file to .gitignore"
    fi
done
