import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

class GenoSignatureScanner:
    """
    A bioinformatics tool for scanning genomic sequences to identify 
    structural and regulatory signatures.
    """
    def __init__(self, window_size=300, step_size=100):
        self.window_size = window_size
        self.step_size = step_size

    def analyze(self, file_path):
        """Extracts GC content and CpG O/E ratios from a FASTA file."""
        try:
            record = SeqIO.read(file_path, "fasta")
            sequence = str(record.seq).upper()
            
            gc_values = []
            oe_ratios = []

            for i in range(0, len(sequence) - self.window_size, self.step_size):
                sub = sequence[i : i + self.window_size]
                g, c, cg = sub.count("G"), sub.count("C"), sub.count("CG")
                
                # GC Content calculation
                gc_values.append((g + c) / self.window_size)
                
                # CpG O/E Ratio calculation
                if g > 0 and c > 0:
                    oe = (cg * self.window_size) / (g * c)
                    # Filter for biologically significant CpG Islands
                    oe_ratios.append(oe if (oe > 0.6 and (g+c)/self.window_size > 0.5) else 0)
                else:
                    oe_ratios.append(0)
            
            return np.array(gc_values), np.array(oe_ratios), record.description
        except Exception as e:
            print(f"Error processing file: {e}")
            return None, None, None

    def plot_signature(self, gc_vals, oe_vals, title):
        """Visualizes the genomic landscape."""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8))
        
        # Plot 1: GC Barcode (Heatmap)
        sns.heatmap(gc_vals.reshape(1, -1), cmap="viridis", ax=ax1, cbar_kws={'label': 'GC%'})
        ax1.set_title(f"GenoSignature-Scanner: GC Content ({title})")
        ax1.set_yticks([])

        # Plot 2: CpG Skyscrapers (Bar Chart)
        ax2.fill_between(range(len(oe_vals)), oe_vals, color="crimson", alpha=0.8)
        ax2.set_title("GenoSignature-Scanner: CpG O/E Ratio (Regulatory Skyscrapers)")
        ax2.set_ylabel("O/E Ratio")
        ax2.axhline(0.6, color='black', linestyle='--', alpha=0.3, label='Threshold (0.6)')
        ax2.legend()
        
        plt.tight_layout()
        plt.show()

# --- Example Usage ---
# --- FINAL INTEGRATED EXAMPLE USAGE ---
if __name__ == "__main__":
    # 1. Initialize the scanner
    scanner = GenoSignatureScanner(window_size=300, step_size=50)

    # 2. RUN E. COLI ANALYSIS (Prokaryotic Example)
    # This shows the "Rolling Hills" / dense genome structure
    print("Starting E. coli Analysis...")
    g_eco, o_eco, n_eco = scanner.analyze("ecoli_genome.fna")
    if g_eco is not None:
        scanner.plot_signature(g_eco, o_eco, "E. coli (Bacterial Architecture)")

    # 3. RUN BRCA1 ANALYSIS (Eukaryotic Example)
    # This shows the "Skyscraper Effect" / regulatory switches
    print("Starting Human BRCA1 Analysis...")
    g_brca, o_brca, n_brca = scanner.analyze("brca1_region.fna")
    if g_brca is not None:
        scanner.plot_signature(g_brca, o_brca, "Human BRCA1 Region (Regulatory Skyscrapers)")
