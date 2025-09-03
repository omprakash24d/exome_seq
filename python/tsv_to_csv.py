import os
import pandas as pd

# Get script folder
base_dir = os.path.dirname(__file__)

# Paths
sam_file = os.path.join(base_dir, "..", "mapped", "father.sam")
csv_file = os.path.join(base_dir, "..", "python", "father_parsed.csv")

# Mandatory SAM columns
mandatory_cols = [
    "QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR",
    "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"
]

rows = []

with open(sam_file, "r") as f:
    for line in f:
        # skip header lines
        if line.startswith("@"):
            continue
        
        parts = line.strip().split("\t")
        
        # first 11 fields
        record = {col: parts[i] for i, col in enumerate(mandatory_cols)}
        
        # optional fields
        for opt in parts[11:]:
            try:
                tag, typ, val = opt.split(":", 2)
                record[tag] = val
            except ValueError:
                continue
        
        rows.append(record)

# Convert to DataFrame (pandas will align columns automatically)
df = pd.DataFrame(rows)

# Save to CSV
df.to_csv(csv_file, index=False)

print(f"✅ CSV saved to {csv_file}")
print(f"   Rows: {len(df)}, Columns: {df.shape[1]}")
