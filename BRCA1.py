from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# --- CONFIGURATION ---
FILENAME = "brca1_region.fna"  # The file you downloaded from NCBI
WINDOW_SIZE = 300              # Human CpG islands are sharp and narrow
STEP_SIZE = 50                 # Small steps give us high-resolution "skyscrapers"

# --- STEP 1 & 2: Load the Human Genomic Data ---
def load_human_data(file_path):
    try:
        record = SeqIO.read(file_path, "fasta")
        sequence = str(record.seq).upper()
        print(f" Successfully Loaded: {record.description}")
        print(f" Sequence Length: {len(sequence):,} base pairs")
        return sequence
    except FileNotFoundError:
        print(f" Error: Could not find '{file_path}'. Please check the filename.")
        return None

dna_seq = load_human_data(FILENAME)

if dna_seq:
    # --- STEP 3 & 6: Calculate GC Content and CpG O/E Ratio ---
    gc_values = []
    oe_ratios = []
    
    print(f" Scanning BRCA1 region for CpG Skyscrapers...")

    for i in range(0, len(dna_seq) - WINDOW_SIZE, STEP_SIZE):
        subseq = dna_seq[i : i + WINDOW_SIZE]
        
        # Calculate GC Content
        g = subseq.count("G")
        c = subseq.count("C")
        gc_perc = (g + c) / WINDOW_SIZE
        gc_values.append(gc_perc)
        
        # Calculate Observed/Expected CpG Ratio
        cg_pairs = subseq.count("CG")
        if g > 0 and c > 0:
            # Formula: (Count of CG * WindowSize) / (Count of C * Count of G)
            oe = (cg_pairs * WINDOW_SIZE) / (g * c)
        else:
            oe = 0
            
        # Biological Filter: Only record the "Skyscraper" if it meets the criteria
        # (O/E Ratio > 0.6 AND GC Content > 50%)
        if oe > 0.6 and gc_perc > 0.5:
            oe_ratios.append(oe)
        else:
            oe_ratios.append(0)

    # --- STEP 4 & 7: Dual Visualization ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), sharex=True)
    
    # Plot 1: GC Content Heatmap (The "Signature")
    # We reshape the GC data for a barcode look
    rows = 1
    data_array = np.array(gc_values).reshape(rows, -1)
    sns.heatmap(data_array, cmap="viridis", ax=ax1, cbar_kws={'label': 'GC %'})
    ax1.set_title("BRCA1 Region: GC Content 'Barcode'")
    ax1.set_yticks([]) # Hide Y-axis for the barcode

    # Plot 2: CpG O/E Ratio (The "Skyscrapers")
    ax2.fill_between(range(len(oe_ratios)), oe_ratios, color="crimson", alpha=0.8)
    ax2.set_title("BRCA1 Region: CpG Island 'Skyscrapers' (Regulatory Switches)")
    ax2.set_ylabel("Observed/Expected Ratio")
    ax2.set_xlabel(f"Genomic Window (Step Size: {STEP_SIZE}bp)")
    ax2.axhline(0.6, color='black', linestyle='--', alpha=0.5, label='O/E Threshold (0.6)')
    ax2.legend()

    plt.tight_layout()
    print(" Displaying Genomic Analysis...")
    plt.show()

    # --- STEP 5: Final Biological Stats ---
    print("\n--- FINAL GENOMIC SUMMARY ---")
    print(f"Global GC Content of BRCA1 region: {(dna_seq.count('G') + dna_seq.count('C')) / len(dna_seq) * 100:.2f}%")
    print(f"Total Windows Scanned: {len(oe_ratios)}")
    print(f"Potential Regulatory Elements Found: {len([x for x in oe_ratios if x > 0])}")
