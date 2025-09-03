#!/bin/bash

# Maximum allowed file size in MB
MAX_SIZE=150

# --- Step 1: Add large files to .gitignore ---
touch .gitignore
find . -type f -size +"${MAX_SIZE}M" | while read -r file; do
    file=${file#./}  # Remove leading ./
    if ! grep -Fxq "$file" .gitignore; then
        echo "$file" >> .gitignore
        echo "Added $file to .gitignore"
    fi
done

# --- Step 2: Create pre-commit hook to block large files ---
HOOK_FILE=".git/hooks/pre-commit"
cat > "$HOOK_FILE" << 'EOF'
#!/bin/bash
MAX_SIZE=10
STAGED_FILES=$(git diff --cached --name-only)

for file in $STAGED_FILES; do
    if [ -f "$file" ]; then
        FILE_SIZE_MB=$(du -m "$file" | cut -f1)
        if [ "$FILE_SIZE_MB" -gt "$MAX_SIZE" ]; then
            echo "Error: $file is $FILE_SIZE_MB MB (larger than $MAX_SIZE MB). Commit aborted."
            exit 1
        fi
    fi
done

exit 0
EOF

chmod +x "$HOOK_FILE"
echo "Pre-commit hook installed. Commits with files larger than $MAX_SIZE MB will be blocked."
