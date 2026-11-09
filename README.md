# ViralJEPA pretraining corpus

`build_corpus.py` builds an **unlabeled** corpus of 300 bp DNA windows for JEPA pretraining.
The output uses the ViraMiner CSV format (`id,sequence,label`), with the label always `-1`.
The numbers below come from **corpus_v3** (`corpus3_cluster/corpus_v3/manifest.json`).
**corpus_v4** is the same recipe plus the BLAST screen (step 4b). Its numbers will be in
`corpus_v4/manifest.json` once the cluster build finishes.

## Usage

```bash
python build_corpus.py download --raw-dir RAW
python build_corpus.py build --raw-dir RAW --viraminer-dir dataset/ --out-dir OUT \
    --mix virus_eukaryotic=0.18,virus_prokaryotic=0.325,human=0.25,bacteria=0.20,viraminer=0.045 \
    --max-per-species 10 --max-genome-weight-bp 50000 --holdout-fraction 0 \
    --blast-screen --blast-jobs 32        # v4 only; v3 was built without these two flags
```

Requirements:
- `numpy`
- `zstandard` (or the `zstd` binary)
- the NCBI `datasets` CLI, for the `ncbi_virus` source
- BLAST+ (`blastn`, `makeblastdb`), for `--blast-screen`. Pass `--blast-bin DIR` if they are not on `PATH`.

On the cluster, use `build_corpus.slurm`. It reuses the existing `/save/corpus_raw` without downloading again, and builds v4 with the BLAST screen on 32 cores in about 2–2.5 h.

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

4b. **BLAST screen (`--blast-screen`, v4 only).** This step removes close homologs of the test reads. It runs after the k-mer screens, on the windows they keep, in every category including ViraMiner.
   - **Search:** the 26,405 `fullset_test` reads form a blastn database, and every corpus window is a query. Parameters: `-task blastn -word_size 7 -evalue 1e-3 -strand both -dust yes -max_target_seqs 500`.
   - **Rule:** a window is dropped if, for at least one test read, its HSPs with identity > 70% and E ≤ 1e-3 cover at least 80% of that read (≥ 240 of 300 bp).
     - Coverage is the union of the HSP intervals on the test read, counted per strand. Two partial alignments that together span the read therefore count.
     - Options to change it: `--blast-identity`, `--blast-coverage`, `--blast-evalue`, `--blast-word-size`, and `--blast-categories` to choose which categories are screened.
   - **E-value:** E depends on database size, so it is not comparable to the opposite search direction (test reads against a corpus database). With this small database the same alignment gets a ~100× smaller E. That doesn't matter here: any alignment covering ≥ 240 bp at > 70% identity is far below 1e-3, so the coverage and identity thresholds decide.
   - **Output:** every dropped window is written to `blast_screen_hits.csv` together with the test read that removed it, the coverage, the best identity and the E-value. The manifest records the rule and the number removed per category, under `windows_screened_out`, `<category>:blast_fullset_test`.
   - **Expected effect,** from a local pilot on 20,000 v3 windows:

     | Category | Dropped |
     |---|---|
     | human | 12.8% |
     | viraminer | 2.2% |
     | virus_eukaryotic | 0.03% |
     | virus_prokaryotic | 0.02% |
     | bacteria | 0% |

     The human windows are not copies of test reads. Each one matches a median of 5 different non-viral test reads at about 83% identity, which points to shared repeat families. Every category is drawn at 1.3× its target, so all of them still reach their targets after this step.
   - **Cost:** 31.8 windows/s per core, about 45 core-hours for the ~5.1M windows left after the k-mer screens.

5. **Write.** Trim each category to its target, shuffle, and write `pretrain_corpus_000…003.csv` (1M rows each) plus `manifest.json`, and `blast_screen_hits.csv` when the BLAST screen ran. v3 has no holdout file.

**Window ID format:** `category|source|record|start|strand`, for example `human|human|chr1|217010|+`.

## Known limitations

- `refseq_viral` includes plant and insect viruses. `uhgv` is labelled phage as a whole, but contains some eukaryotic gut viruses.
- The species cap is taxonomic, not sequence clustering, so near-identical strains of one species can still repeat.
- In v4, the BLAST screen removes repeat-rich human windows, so the human category has fewer repeats than the genome itself. To leave human windows out of the BLAST screen, pass `--blast-categories virus_eukaryotic,virus_prokaryotic,bacteria,viraminer`.
- All bacteria come from the gut (UHGG). There are no skin or cervical genomes, although ViraMiner samples include those sites.
- Sources not used:
  - IMG/VR: needs a JGI account and is mostly environmental.
  - OpenGenome2: eukaryotic viruses were removed from it.
