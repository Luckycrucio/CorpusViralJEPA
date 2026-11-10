# ViralJEPA pretraining corpus

`build_corpus.py` builds an **unlabeled** corpus of 300 bp DNA windows for JEPA pretraining.
The output uses the ViraMiner CSV format (`id,sequence,label`), with the label always `-1`.
The numbers below come from **corpus_v3** (`corpus3_cluster/corpus_v3/manifest.json`).
**corpus_v4** is the same recipe plus the BLAST screen (step 4b). The section
[corpus_v4: what BLASTn changed](#corpus_v4-what-blastn-changed) reports how it differs from v3.

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
   - **Effect:** see [corpus_v4: what BLASTn changed](#corpus_v4-what-blastn-changed). Every category is drawn at 1.3× its target, so all of them still reach their targets after this step.
   - **Cost:** in v4, 113 min on 32 cores for the 5.05M windows left after the k-mer screens (about 23 windows/s per blastn process).

5. **Write.** Trim each category to its target, shuffle, and write `pretrain_corpus_000…003.csv` (1M rows each) plus `manifest.json`, and `blast_screen_hits.csv` when the BLAST screen ran. v3 has no holdout file.

**Window ID format:** `category|source|record|start|strand`, for example `human|human|chr1|217010|+`.

## corpus_v4: what BLASTn changed

corpus_v4 was built on cluster taking 2 h 14 min (mostly for BLASTn), it's log and results are in `corpus4_cluster/corpus_v4/`. Up to the BLAST screen, the build is identical to v3: same raw files, record counts, targets and k-mer screen counts. Both corpora have 4,000,000 windows with the same total per category.

| Category | Removed by BLAST in v4 | v3 windows that fail the rule | Removals caused by viral test reads | Windows shared with v3 |
|---|---|---|---|---|
| human | 161,897 | 134,065 (13.4%) | 0 | 83% |
| viraminer | 4,231 | 4,108 (2.3%) | 424 | 97% |
| virus_eukaryotic | 414 | 327 (0.05%) | 383 | 77% |
| bacteria | 393 | 303 (0.04%) | 0 | 77% |
| virus_prokaryotic | 145 | 122 (0.01%) | 73 | 77% |

v3 and v4 draw from the same window pool, so each v3 window was checked against v4's `blast_screen_hits.csv`. That gives the "fail the rule" column. No v4 window fails the rule.

- **Human:** these windows are repeat homologs, not copies of test reads. Median identity is 78% over 94% of the read. The removals come from 1,020 non-viral test reads, and a few reads account for most of them (the top one alone removes 14,087 windows). Only 571 hits are ≥ 95% identical. Mean GC barely changes (0.405 → 0.403), but v4's human category has fewer repeats than the genome itself.
- **Eukaryotic viruses:** few in number, but these are the removals that matter for the viral task. 383 of the 414 come from viral test reads, at 70–83% identity, none ≥ 95%.
- **ViraMiner:** 381 hits are ≥ 95% identical, i.e. pieces of the same contigs as the test reads.
- **Overall:** 2,127 test reads caused at least one removal, and 160 of them are viral.
- **v4 is not v3 minus the flagged windows.** After screening, each category is randomly trimmed to its target, and that draw differs from v3's. So only 77–83% of the genome-category windows are the same as in v3, and 97% of ViraMiner windows.

## Known limitations

- `refseq_viral` includes plant and insect viruses. `uhgv` is labelled phage as a whole, but contains some eukaryotic gut viruses.
- The species cap is taxonomic, not sequence clustering, so near-identical strains of one species can still repeat.
- In v4, the BLAST screen removes repeat-rich human windows, so the human category has fewer repeats than the genome itself. To leave human windows out of the BLAST screen, pass `--blast-categories virus_eukaryotic,virus_prokaryotic,bacteria,viraminer`.
- All bacteria come from the gut (UHGG). There are no skin or cervical genomes, although ViraMiner samples include those sites.
- Sources not used:
  - IMG/VR: needs a JGI account and is mostly environmental.
  - OpenGenome2: eukaryotic viruses were removed from it.
