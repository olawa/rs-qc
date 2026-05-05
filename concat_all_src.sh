#!/bin/bash
# concat_src.sh - Concatenate rsnap source code into a single text file

set -euo pipefail

OUTPUT_FILE="rs-qc_full_source.txt"
PROJECT_ROOT="$(pwd)"

echo "Generating $OUTPUT_FILE from $PROJECT_ROOT ..."
echo "================================================================================" > "$OUTPUT_FILE"
echo "rsnap Full Source Concatenation" >> "$OUTPUT_FILE"
echo "Generated on: $(date)" >> "$OUTPUT_FILE"
echo "Project Root: $PROJECT_ROOT" >> "$OUTPUT_FILE"
echo "================================================================================" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

# Find and concatenate files
find . \
    -type f \
    \( -name "*.rs" -o -name "Cargo.toml" -o -name "Cargo.lock" \
       -o -name "*.md" -o -name "*.sh" -o -name "*.toml" \
       -o -name "*.rs" \) \
    -not -path "./target/*" \
    -not -path "./.git/*" \
    -not -path "./data/*" \
    -not -path "./region_plot/target/*" \
    -not -path "*/__pycache__/*" \
    -not -name "*.png" \
    -not -name "*.sh" \
    -not -name "*.txt" \
    -not -name "*.gz" \
    -not -name "*.py" \
    -not -name "*.md" \
    -not -name "rsnap" \
    -not -name "*.tar.gz" \
    -not -name "concat_src.sh" \
    -not -name "make_tar.sh" \
    -not -name "rsnap_full_source.txt" \
    | sort | while read -r file; do
    echo "================================================================================" >> "$OUTPUT_FILE"
    echo "FILE: ${file#./}" >> "$OUTPUT_FILE"
    echo "================================================================================" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
    cat "$file" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
done

echo "================================================================================" >> "$OUTPUT_FILE"
echo "End of rsnap source concatenation" >> "$OUTPUT_FILE"
echo "Total files processed: $(grep -c "^FILE:" "$OUTPUT_FILE")" >> "$OUTPUT_FILE"

echo "✅ Done! Created $OUTPUT_FILE ($(du -h "$OUTPUT_FILE" | cut -f1))"
