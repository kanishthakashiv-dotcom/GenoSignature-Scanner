from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# --- STEP 1 & 2: Load the Genome ---
def load_genome(file_path):
    try:
        record = SeqIO.read(file_path, "fasta")
        return str(record.seq).upper()
    except FileNotFoundError:
        print(f"Error: Could not find {file_path}")
        return None

# !!! CHANGE THIS to your actual filename !!!
FILENAME = "GCF_000005845.2_ASM584v2_genomic.fna"
genome_dna = load_genome(FILENAME)

if genome_dna:
    # --- STEP 3: Sliding Window Calculation ---
    window_size = 1000
    step_size = 500
    gc_values = [] # <--- THIS IS WHAT WAS MISSING

    for i in range(0, len(genome_dna) - window_size, step_size):
        subseq = genome_dna[i : i + window_size]
        gc_content = (subseq.count("G") + subseq.count("C")) / window_size
        gc_values.append(gc_content)

    # --- STEP 4: Heatmap Visualization ---
    # Now gc_values EXISTS, so Step 4 will work!
    data_array = np.array(gc_values)
    rows = 80
    cols = len(data_array) // rows
    trimmed_data = data_array[:rows * cols].reshape((rows, cols))

    plt.figure(figsize=(12, 6))
    sns.heatmap(trimmed_data, cmap="viridis")
    plt.title("Final GC Content Heatmap")
    plt.show()
    
# --- STEP 5: Biological Stats ---
total_g = genome_dna.count("G")
total_c = genome_dna.count("C")
total_len = len(genome_dna)
global_gc = (total_g + total_c) / total_len * 100

print("\n--- BIOLOGICAL SUMMARY ---")
print(f"Species: {FILENAME}")
print(f"Global GC Content: {global_gc:.2f}%")

# Check against E. coli baseline
if 50 <= global_gc <= 52:
    print("Validation: This matches the expected range for E. coli K-12.")
else:
    print("Validation: GC content deviates from E. coli baseline. Check species or data.")
    
# --- STEP 6: CpG Island Detection ---
window_size = 200  # CpG islands are smaller than general GC shifts
step_size = 100
cpg_islands = []

for i in range(0, len(genome_dna) - window_size, step_size):
    subseq = genome_dna[i : i + window_size]
    
    # 1. Calculate GC %
    g_count = subseq.count("G")
    c_count = subseq.count("C")
    gc_content = (g_count + c_count) / window_size
    
    # 2. Calculate Observed/Expected Ratio
    # We count "CG" pairs specifically
    cg_pairs = subseq.count("CG")
    
    # Avoid division by zero if G or C is missing
    if g_count > 0 and c_count > 0:
        oe_ratio = (cg_pairs * window_size) / (g_count * c_count)
    else:
        oe_ratio = 0

    # 3. Apply the Biological Thresholds
    if gc_content > 0.50 and oe_ratio > 0.6:
        # Save the index and the ratio for plotting
        cpg_islands.append(oe_ratio)
    else:
        # If it doesn't meet criteria, we mark it as 0
        cpg_islands.append(0)

print(f"Found {len([x for x in cpg_islands if x > 0])} potential CpG windows!")

plt.figure(figsize=(15, 5))
plt.fill_between(range(len(cpg_islands)), cpg_islands, color="red", alpha=0.7)
plt.title("CpG Island Locations (O/E Ratio > 0.6)")
plt.ylabel("Observed/Expected Ratio")
plt.xlabel("Genomic Window")
plt.show()
