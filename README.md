# Velvet Worm Diet Analysis from Midsection RNA-seq

**Repository:** [gyagu98/velvet_worm_diet_2](https://github.com/gyagu98/velvet_worm_diet_2)  
**Research area:** Computational genomics · diet-associated sequence screening · Onychophora  
**Maintainer:** Genrietta Y. Rozenberg (née Yagudayeva)

## Project overview

This repository preserves Python and R scripts developed during an exploratory analysis of **diet-associated sequence signals in wild-collected Andean velvet worms** (Onychophora: Peripatidae). The study uses RNA-seq-derived midsection assemblies, collected with gut contents retained, to search for potential non-host mitochondrial **cytochrome c oxidase subunit I (COX1)** sequences and other taxonomic matches.

The objective is to distinguish plausible non-host signals from **velvet worm host sequences, short or ambiguous matches, environmental DNA/RNA, and possible laboratory contamination**. Sequence detection does **not**, on its own, demonstrate that a particular taxon was eaten. The analyses and biological interpretations remain subject to validation.

## Samples

| Sample | Source | General analysis |
| --- | --- | --- |
| **F734** | Wild velvet worm midsection, gut contents retained | COX1 candidate screening and filtering |
| **F739** | Wild velvet worm midsection, gut contents retained | Candidate screening and taxonomic review |
| **F744** | Wild velvet worm midsection, gut contents retained | COX1 screening, contamination review, and selected read-mapping checks |

Specimens were dissected promptly, with separate thoroughly cleaned instruments. These precautions do not exclude all sources of contamination.

## Repository organization

The project evolved through multiple trials. **Scripts are retained under their original filenames and trial names**, rather than presented as one fully automated production pipeline.

```text
velvet_worm_diet_2/
├── README.md
├── .gitignore
└── scripts/
    ├── F734/                      # Initial F734 candidate-screening scripts
    ├── Trial_2/                   # Early COX1 reference retrieval and screening
    ├── Trial_3/                   # Per-sample screens, BLAST input preparation, plots
    │   ├── F734/
    │   ├── F739/
    │   ├── F744/
    │   └── Final_results/
    ├── Trial_4/                   # Candidate trimming, validation, visualizations
    ├── Trial_5/                   # Additional BLAST checks and figure revisions
    ├── Manuscript_v1/            # Manuscript-associated plotting script
    ├── velvet_worm_diet/         # Historical nested copies of earlier scripts
    └── original_layout/          # Backup of the earlier directory structure
```

Some paths contain **duplicate historical copies**. A higher script-file count does not mean that every file is a unique analysis. Not every folder shown above is guaranteed to be present in each clone until all local changes have been committed and pushed.

## Analysis strategy

The archived scripts support portions of the following exploratory workflow (specific actions varied across trials):

1. **Prepare assembled sequence inputs** from velvet worm midsection RNA-seq data.
2. **Collect and search reference sequences**, particularly animal COX1 markers, using reference retrieval and alignment tools.
3. **Screen candidate sequences** for potential taxonomic assignments to invertebrates or other organisms.
4. **Filter likely host sequences** by comparison against onychophoran references and available host mitochondrial sequences.
5. **Review contamination signals**, especially human or vertebrate matches; snake-like hits were treated as contamination candidates rather than dietary evidence.
6. **Validate shortlisted candidates** by examining sequence identity, alignment length, coverage, E-values, alternative hits, and taxonomic specificity.
7. **Check read support where applicable**, including Bowtie2 mapping and SAMtools-based inspection for selected F744 candidates.
8. **Prepare screening summaries and visualizations** for manuscript development.

### Selected scripts

| Script or group | Role |
| --- | --- |
| `scripts/F734/prey_cox1_screen.py` and `prey_cox1_screen_v2.py` | F734 COX1 candidate screening versions |
| `scripts/Trial_2/download_cox1.py` | COX1 reference retrieval |
| `scripts/Trial_2/velvet_prey_cox1.py` | Early prey-associated COX1 screening |
| `scripts/Trial_2/F744/velvet_prey_cox1_v2.py` | Revised F744 screening |
| `scripts/Trial_3/F734/`, `F739/`, `F744/` | Trial 3 sample-specific screens |
| `scripts/Trial_3/Final_results/make_blast_fasta.py` | FASTA preparation for subsequent BLAST assessment |
| `scripts/Trial_4/trim_for_blast.py`, `trial4_validate.py` | Candidate trimming and validation |
| `scripts/Trial_5/trial5_blast_verify.py` | Later BLAST verification |
| `scripts/Trial_3/Final_results/*figures*.R`, `scripts/Trial_4/*figures*.R`, `scripts/Trial_5/*figures*.R` | Figure generation and revisions |

Script roles are descriptive; inspect each script's arguments, dependencies, and hard-coded paths before running.

## Software and reproducibility

Tools used across the broader analyses include **Python, R, BLAST+, MMseqs2, Bowtie2, and SAMtools**. Work was performed on macOS and Linux systems. Different trials may have used different dependencies and configurations.

This repository is an **analysis archive, not yet a turn-key reproducible workflow**. In particular:

- Scripts may reference local input/output paths that must be changed for another machine.
- Input assemblies, raw reads, and large reference databases are not included in this code backup.
- Trial folders represent iterations rather than a definitive execution order.
- Exact software versions, commands, reference accessions, and filtering cutoffs must be verified against scripts and original logs before citing them as final manuscript methods.

For reproduction, users need the appropriate input files, software environments, and trial-specific command history.

## Data and size limits

This GitHub repository is intended primarily for **source code and lightweight documentation**, not bulk sequencing data. Very large FASTA assemblies and MMseqs2 outputs are kept outside normal Git history. In particular, the F734 and F744 assembly/MMseqs2 FASTA files that exceed GitHub's 100 MB per-file limit should not be committed to this repository.

Small reviewed summary tables can be added separately, with relevant sample identifiers and explanatory metadata. Raw reads and full assemblies should ultimately be shared using an appropriate scientific data repository when publication and permissions allow.

## Scientific interpretation

An apparent non-host COX1 match in a midsection transcriptome is **not definitive prey identification**. Alternative explanations include conserved or short alignments, host mitochondrial sequence, incomplete reference databases, contamination, or nucleic acids present in the gut without confirmed consumption. Confidence requires independent validation, robust alignments, appropriate controls, and reproducible read support.

## Citation and reuse

The manuscript and data availability details will be updated if and when the study is published. No open-source license has been selected; absent an explicit license, standard copyright applies. Contact the maintainer regarding reuse.
