# ViralJEPA pretraining corpus

`build_corpus.py` builds an **unlabeled** corpus of 300 bp DNA windows for JEPA pretraining.
The output uses the ViraMiner CSV format (`id,sequence,label`), with the label always `-1`.
The numbers below come from **corpus_v3** (`corpus3_cluster/corpus_v3/manifest.json`).

## Usage

```bash
python build_corpus.py download --raw-dir RAW
python build_corpus.py build --raw-dir RAW --viraminer-dir dataset/ --out-dir OUT \
    --mix virus_eukaryotic=0.18,virus_prokaryotic=0.325,human=0.25,bacteria=0.20,viraminer=0.045 \
    --max-per-species 10 --max-genome-weight-bp 50000 --holdout-fraction 0
```

Requirements:
- `numpy`
- `zstandard` (or the `zstd` binary)
- the NCBI `datasets` CLI, for the `ncbi_virus` source

On the cluster, use `build_corpus.slurm`.

## Sources

| Source | Downloaded from | Kept | Category | Records / Mbp kept | Windows |
|---|---|---|---|---|---|
| `refseq_viral` | `ftp.ncbi.nlm.nih.gov/refseq/release/viral/`, all `viral.*.genomic.fna.gz` | all, minus 175 anellovirus and 139 genomes also in `ncbi_virus` | header contains "phage" → `virus_prokaryotic`, otherwise `virus_eukaryotic` | 19,275 / 580 | shared, see below |
| `ncbi_virus` | `datasets download virus genome taxon <T> --host human --complete-only`; 18 taxa, listed below | ≤ 10 genomes per ICTV species (53,487 capped); Anelloviridae dropped | `virus_eukaryotic` | 1,266 / 33 | shared, see below |
| `uhgv` | Zenodo 17402089: `votus_full.fna.zst` + `votus_metadata_extended.tsv.zst` | vOTU representatives with CheckV quality `Complete` or `High-quality`; 927 Anelloviridae dropped | `virus_prokaryotic` | 57,734 / 2,598 | shared, see below |
| `uhgg` | EBI MGnify human-gut v2.0.2: `genomes-all_metadata.tsv`, then `species_catalogue/.../<id>.fna` | species representatives with completeness ≥ 90% and contamination ≤ 5%; 1,000 genomes, picked round-robin across phyla (seed 0) | `bacteria` | 92,449 contigs / 2,612 | 800,000 |
| `human` | T2T-CHM13v2.0 analysis set (`chm13v2.0.fa.gz`) | all 25 sequences | `human` | 25 / 3,117 | 1,000,000 |
| ViraMiner | `dataset/fullset_train.csv`, labels dropped | random rows | `viraminer` | — | 180,000 |

How the shared rows split:
- **`virus_eukaryotic` (720,000 windows):** the non-phage part of `refseq_viral`, plus all of `ncbi_virus`.
- **`virus_prokaryotic` (1,300,000 windows):** the phage part of `refseq_viral`, plus all of `uhgv`.

**NCBI Virus taxa (host = human):**
- Families: Papillomaviridae, Polyomaviridae, Anelloviridae (downloaded but not used), Orthoherpesviridae, Adenoviridae, Parvoviridae, Circoviridae, Genomoviridae, Hepadnaviridae, Poxviridae, Picornaviridae, Caliciviridae, Astroviridae, Flaviviridae.
- Retroviruses: HIV-1, HIV-2, HTLV-1 and HTLV-2, downloaded one by one because the whole Retroviridae package fails.

The species is the lineage entry below the genus in each genome's `data_report.jsonl`. RefSeq has no cap of its own.

## Build steps

1. **Index.** Drop records shorter than 300 bp and apply the filters in the table above.
2. **Allocate.** Each category gets `num_windows × mix` windows (total 4,000,000). Inside a category, windows are spread over genomes in proportion to genome length. That length is capped at 50 kb for every category except human, so jumbo phages and poxviruses don't dominate. The draw takes 1.3× the target, to cover later rejections.
3. **Extract.** Take random 300 bp windows. Windows with characters other than ACGT are redrawn up to 4 times. Lower-case repeats are upper-cased and kept. Each window is reverse-complemented with probability 0.5.
4. **Screen.** Remove leakage using canonical k-mers. A window is dropped if at least 10% of its k-mers appear in the screen set:

   | Screen | Applied to | Screen set | k | Removed in v3 |
   |---|---|---|---|---|
   | `evaluation_standard` | external genomes | `fullset_validation` + `fullset_test` + `test_data_anell_300_binary` | 31 | human 93,397; virus_euk 970; virus_prok 534; bacteria 439 |
   | `anell_test_overlap` | ViraMiner windows | `test_data_anell_300_binary`, which re-splits the same pool as `fullset_train` | 31 | 25,829 |
   | `anello_unseen_strict` | all windows | anellovirus positives in `test_data_anell_300_binary` | 15 | 447 |

   Anelloviruses are removed by taxonomy, by the header regex (`anellovir|torque teno|TTV|…`) and by this screen. That keeps the anell test an unseen-family benchmark.
5. **Write.** Trim each category to its target, shuffle, and write `pretrain_corpus_000…003.csv` (1M rows each) plus `manifest.json`. v3 has no holdout file.

**Window ID format:** `category|source|record|start|strand`, for example `human|human|chr1|217010|+`.

## Known limitations

- `refseq_viral` includes plant and insect viruses. `uhgv` is labelled phage as a whole, but contains some eukaryotic gut viruses.
- The species cap is taxonomic, not sequence clustering, so near-identical strains of one species can still repeat.
- All bacteria come from the gut (UHGG). There are no skin or cervical genomes, although ViraMiner samples include those sites.
- Sources not used:
  - IMG/VR: needs a JGI account and is mostly environmental.
  - OpenGenome2: eukaryotic viruses were removed from it.
