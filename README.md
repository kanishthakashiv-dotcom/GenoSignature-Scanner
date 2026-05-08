# GenoSignature-Scanner
The GenoSignature-Scanner is a Python-based bioinformatics tool for extracting localized genomic features. Using Biopython and customizable sliding windows, it calculates GC-density and CpG O/E ratios. The repository includes validated examples identifying bacterial genomic structures and human BRCA1 regulatory "skyscrapers"—a vital tool for gene "on-switch" detection.
It is basically a modular Python-based bioinformatics tool designed for the automated detection of genomic signatures and regulatory elements. By utilizing a high-resolution sliding window methodology, the tool extracts localized chemical and statistical properties of DNA sequences.

🚀 Core Functionality
GenoSignature-Scanner processes genomic data in FASTA/FNA format to calculate two primary biological metrics:
1. GC-Density Mapping: Visualizes the percentage of G-C bases to identify regions of high thermal stability and potential gene density.
2. CpG Island Detection: Calculates the Observed/Expected (O/E) Ratio of CG dinucleotides. It identifies functional "on-switches" (promoters) based on established biological thresholds (O/E > 0.6 and GC > 50%).

🔬 Research Examples Included
The repository provides two pre-configured datasets to demonstrate the tool's versatility:
1. Prokaryotic Case Study (E. coli): Illustrates the dense, uniform coding architecture characteristic of bacterial genomes.
2. Eukaryotic Case Study (Human BRCA1): Demonstrates the "Skyscraper Effect," pinpointing the BRCA1 promoter region within a highly suppressed genomic background.

🛠️ Requirements
- Python 3.8+
- Biopython
- Matplotlib / Seaborn
- NumPy

🔬 The Results
1. The E. coli Result :
Because E. coli is a bacterium, its DNA is packed tightly with genes. Your graph will look like a series of rolling hills. There isn't much "empty space," so the signal is scattered across the whole genome.

2. The BRCA1 Result :
Because this is human DNA, most of the sequence is "quiet" (suppressed). You will see a lot of zeros, followed by a massive vertical spike (the "Skyscraper"). This spike is the Promoter, the molecular "On-Button" for the BRCA1 gene.
