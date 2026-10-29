# ViralJEPA pretraining corpus

`build_corpus.py` builds a large **unlabeled** corpus of 300 bp DNA windows for JEPA
pretraining. The windows have the same format as the ViraMiner files (`fullset_*.csv`).
They come from public genome collections chosen to match what ViraMiner samples
contain: human DNA, bacteria, viruses that infect humans, and phages.

The corpus exists to fix the problem found in the exp01 diagnostics
(`diagnostics/diag_exp01.py`). Pretraining only on `fullset_train` (63 Mbp, 98% non-viral)
produced an encoder no better than random initialisation on a frozen linear probe. The two
levers here are **more data** and **far more viral diversity**, at no labelling cost.

---

## 1. Quick start

```bash
# 1) Download raw collections (about 6–8 GB; network-bound)
python corpus/build_corpus.py download --raw-dir /path/to/

# 2) Build the corpus (CPU and RAM; run as a cluster job, not on the login node)
python corpus/build_corpus.py build \
    --raw-dir     /path/to \
    --viraminer-dir ~/path/to/viraminer/dataset \
    --out-dir     /path/to \
    --num-windows 4000000 \
    --screen standard
```

Dependencies are `numpy` plus the Python standard library. Reading `.zst` files needs either
the `zstandard` Python module or the `zstd` binary. The NCBI Virus source also needs the
[NCBI `datasets` CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/).
Without it, that source is skipped with a message.

## 2. Sources

| Source key | What | Where it comes from | Size (download) | Corpus category |
|---|---|---|---|---|
| `refseq_viral` | All NCBI RefSeq viral reference genomes | `https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/` (every `viral.*.genomic.fna.gz`) | ~180 MB gz; **19,414 records, 583 Mbp** (Oct 2026) | `virus_eukaryotic` (`virus_prokaryotic` if the header says "phage") |
| `ncbi_virus` | Complete genomes of selected virus families with **host = human** (many strains per species) | NCBI `datasets download virus genome taxon <family> --host human --complete-only` | varies by family | `virus_eukaryotic` |
| `uhgv` | Unified Human Gut Virome catalogue v1.0, **vOTU representatives** (one genome per species-level cluster) | Zenodo record [17402089](https://zenodo.org/records/17402089): `votus_full.fna.zst` + `votus_metadata_extended.tsv.zst` | ~965 MB + 14 MB; 168,536 vOTUs, of which ~58.7k are Complete/High-quality | `virus_prokaryotic` |
| `uhgg` | Unified Human Gastrointestinal Genome catalogue v2.0.2, **species representatives** | EBI MGnify: `.../human-gut/v2.0.2/species_catalogue/<id>/genome/<id>.fna` | metadata 115 MB, then ~3 MB per genome (default 1000 genomes) | `bacteria` |
| `human` | T2T-CHM13v2.0 complete human genome (analysis set) | `s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz` | 981 MB | `human` |
| (built in) | ViraMiner `fullset_train.csv`, **labels ignored** | your `--viraminer-dir` | — | `viraminer` |

Default NCBI Virus taxa (`--ncbi-taxa`): Papillomaviridae, Polyomaviridae,
Anelloviridae, Orthoherpesviridae, Adenoviridae, Parvoviridae, Circoviridae, Genomoviridae,
Hepadnaviridae, Poxviridae, Picornaviridae, Caliciviridae, Astroviridae, Flaviviridae, and the
human retroviruses HIV-1, HIV-2, HTLV-1 and HTLV-2. These are the families seen in human
clinical metagenomes (serum, skin, cervix, condylomata), which are the ViraMiner sample types.
Retroviridae is requested species by species because NCBI returns the whole-family package
(tens of thousands of HIV-1 genomes) as an invalid zip. Each taxon is saved as
`<taxon>.fna` plus its `<taxon>.jsonl` data report, which holds the taxonomy used by the
species cap; a taxon whose report is missing is downloaded again. Failed downloads are
retried 3 times.

**Sources deliberately not used**

- **IMG/VR v4 / MetaVR.** These are huge and mostly environmental phages, and downloading
  needs a JGI account. UHGV already covers human-associated phages. You can add a
  human-associated IMG/VR slice later as an extra source.
- **OpenGenome2 (Evo 2).** Its dataset card says eukaryotic viruses were removed on purpose,
  and those are exactly the viruses this task needs.

---

## 3. What `download` does

`python build_corpus.py download --raw-dir RAW [--sources ...] [--uhgg-genomes 1000] [--workers 8]`

- Every file goes to `RAW/<source>/` through a `.part` temporary file that is renamed when
  complete. Files that already exist are skipped, so re-running resumes an interrupted
  download. Failed downloads are retried 3 times with back-off.
- **`refseq_viral`:** parses the FTP index page and downloads every `viral.N.M.genomic.fna.gz`.
  RefSeq occasionally splits into more files; the index parsing handles that.
- **`ncbi_virus`:** runs the `datasets` CLI once per family, extracts `genomic.fna` from each
  zip into `RAW/ncbi_virus/<Family>.fna`, then deletes the zip. A family that fails is
  logged and skipped.
- **`uhgv`:** downloads the vOTU representative FASTA and the extended vOTU metadata, which
  includes CheckV quality and ICTV taxonomy.
- **`uhgg`:**
  1. Downloads `genomes-all_metadata.tsv`.
  2. Keeps species representatives (`Genome == Species_rep`) with completeness ≥ 90% and
     contamination ≤ 5%.
  3. Picks `--uhgg-genomes` of them **round-robin across phyla**, so Firmicutes and
     Bacteroidota don't swamp the rarer phyla. Seeded and reproducible; the chosen IDs are
     written to `RAW/uhgg/selected_genomes.txt`.
  4. Downloads them in parallel with `--workers` threads.
- **`human`:** downloads the CHM13 FASTA.

---

## 4. What `build` does, step by step

### 4.1 Indexing (pass 1)
Each source's FASTA files are streamed once, and only the record ID, length and category are
kept. That keeps memory low even for multi-GB catalogues. The filters applied here:

| Filter | Applies to | Rule |
|---|---|---|
| Minimum length | all | record ≥ 300 bp |
| UHGV quality | `uhgv` | `checkv_quality` in `--uhgv-quality` (default `Complete,High-quality`) |
| Per-species cap | `ncbi_virus` | at most `--max-per-species` records (default 10) per species, drawn at random (seeded). The species is the ICTV species taxid from the data report's lineage (the entry after the genus), so old and new names ("Monkeypox virus" / "Orthopoxvirus monkeypox", "WU Polyomavirus" / "Betapolyomavirus quartihominis") count as one species. Genomes without a report fall back to the FASTA description cut at "isolate / strain / clone / …". This stops HIV, HBV or HPV16, with thousands of near-identical genomes, from dominating. |
| RefSeq duplicates | `refseq_viral` | RefSeq genomes that also appear in the NCBI Virus reports are dropped from `refseq_viral`, so they count once, under the species cap. RefSeq is otherwise uncapped (it has about one genome per species, and segmented viruses need all their segments). |
| Anelloviridae exclusion | all viral sources | see §4.5 |

### 4.2 Mixture and window budget
`--num-windows` (default 4,000,000 = 1.2 Gbp) is split across categories by `--mix`:

```
virus_eukaryotic=0.30, virus_prokaryotic=0.15, human=0.25, bacteria=0.20, viraminer=0.10
```

Viruses make up 45% of the corpus, against 2% in ViraMiner. **This oversampling is the main
class-imbalance fix on the pretraining side.** The encoder spends capacity on viral sequence
diversity without needing any labels. Human DNA and bacteria stay well represented because
the classifier must reject them, and they are the hard negatives.

If a category has no data (e.g. a source wasn't downloaded), it is dropped and the remaining
weights are renormalised. The manifest records this.

### 4.3 Allocating windows to genomes
Within each category, windows are allocated with a multinomial draw, weighted by genome
length. For every category except `human`, the weight is **capped at
`--max-genome-weight-bp`** (default 50 kb). Without the cap, a few giant genomes (jumbo
phages, poxviruses, 6 Mb bacteria) would dominate, and small viruses (2–8 kb papillomaviruses,
polyomaviruses, circoviruses) would hardly appear. With it, a 5 kb virus is sampled at its
length and anything larger than 50 kb is treated as 50 kb. Human chromosomes are sampled
in proportion to their true length.

To leave room for later rejections, the draw requests `--oversample` × the target
(default 1.3).

### 4.4 Window extraction (pass 2)
- Each file is streamed again, and only the selected records are decoded.
- Start positions are uniform. Windows containing anything other than A/C/G/T (e.g. `N`
  gaps, IUPAC codes) are rejected and redrawn up to 4 times, because ViraMiner contains no N.
  Soft-masked (lower-case) human repeats are upper-cased and **kept**: Alu and LINE fragments
  are a real part of the negative class.
- Each window is flipped to its reverse complement with probability 0.5.
- **Optional:** `--substitution-rate r` mutates a fraction `--noise-fraction` of the windows
  at rate `r`, to mimic sequencing and assembly errors (compare `dataset/SIM/noise_*`).
  It is **off by default**: it's better to add noise as an on-the-fly augmentation during
  training, so every epoch sees fresh noise.

### 4.5 Excluding Anelloviridae (on by default)
`test_data_anell_300_binary.csv` is ViraMiner's **unseen-family** benchmark: anelloviruses
were removed from training and appear only in this test set. For that benchmark to stay
meaningful after pretraining, the corpus must not contain anelloviruses either. They are
removed in three layers:

1. **Taxonomy:** UHGV vOTUs whose `ictv_taxonomy` or `family_name` contains `Anelloviridae`,
   and the `Anelloviridae` file from NCBI Virus, are dropped. The test run dropped 927 UHGV
   vOTUs.
2. **Names:** any viral FASTA header matching
   `anellovir|torque teno|torquetenovirus|TTV|TTMV|TTMDV|gyrovirus|chicken anemia`
   is dropped. This removed 175 RefSeq records.
3. **Sequence:** any window sharing ≥10% of its 15-mers with an anellovirus *positive* in
   `anell_test` is dropped. This is the "strict" screen in §4.6, which catches unnamed
   relatives at roughly ≥90% identity.

`--keep-anelloviridae` turns all three off. Use it only for a corpus that won't be evaluated
on `anell_test`.

### 4.6 Screening against evaluation data (decontamination)
Unlabeled sequences from the test set (or ones nearly identical to them) shouldn't appear in
pretraining. Screening uses **canonical k-mers**: 2-bit packed, and strand-independent
because each k-mer is replaced by the smaller of itself and its reverse complement. A sorted
`uint64` array of all k-mers in the evaluation sequences is built once. Every candidate
window is then checked with a vectorised binary search.

A window is rejected if the **fraction of its k-mers found in the evaluation set** is at
least the threshold:

| `--screen` | k | threshold | What it removes |
|---|---|---|---|
| `none` | — | — | nothing (only the Anelloviridae screen runs) |
| `standard` (default) | 31 | ≥ 10% | near-exact overlaps: any window sharing about 60+ identical bases with an evaluation sequence |
| `strict` | 15 | ≥ 10% | anything at roughly ≥ 90% identity. At 10% divergence about 20% of 15-mers survive (0.9¹⁵), giving about 57 shared k-mers, while random 15-mer matches stay around 1–2% per window. |

**Which sequences are screened against**
- External genomes are screened against `fullset_validation` + `fullset_test` +
  `test_data_anell_300_binary`. The `SIM/noise_*` sets are mutated copies of `fullset_test`,
  so they're covered too.
- The `viraminer` category (windows from `fullset_train`) is **not** screened against the
  fullset val/test splits. Fine-tuning trains on all of `fullset_train` with labels anyway,
  so using it unlabeled in pretraining adds no leakage. It *is* screened against
  `anell_test` (standard preset), because `anell_test` re-splits the same pool: about 12% of
  `fullset_train` sequences appear in `anell_test`.

**Which mode to use.** Report results for both:
- `standard` is the realistic setting. Public reference databases contain the known viruses,
  and ViraMiner's labels were themselves made with BLAST against them. The model can benefit
  from having seen related genomes, but never the test sequences themselves.
- `strict` shows the gain isn't just memorisation of close relatives. It is much more
  aggressive. In the test run it removed **~29% of eukaryotic-virus windows, ~17% of bacterial
  windows, and 0.1% of human windows**, because the test positives are known viruses. Raise
  `--oversample` (e.g. 1.8) with `strict` so categories still reach their targets.

### 4.7 Trimming, shuffling, writing
- Each category is trimmed at random to its target. A category below target is kept as is,
  with a warning.
- Everything is shuffled together. A fraction `--holdout-fraction` (default 1%) is written to
  `pretrain_holdout.csv`; with `0` no holdout file is written. The holdout is drawn per window,
  so its windows come from genomes (and strains) also in training: in corpus_v2, 78% of
  eukaryotic-virus holdout windows overlapped a training window of the same genome. Its loss
  therefore does not measure generalisation. Monitor pretraining with a frozen probe on
  `fullset_validation` instead.
- The rest goes into shards `pretrain_corpus_000.csv`, `_001.csv`, … of `--shard-size` rows
  (default 1M).

---

## 5. Output

```
OUT/
├── pretrain_corpus_000.csv   # training shards
├── pretrain_corpus_001.csv
├── ...
├── pretrain_holdout.csv      # only if --holdout-fraction > 0 (see 4.7)
└── manifest.json             # everything needed to reproduce and report the corpus
```

**CSV format:** no header and three columns, exactly like `fullset_*.csv`. The label is
always `-1`, since there are no labels:

```
virus_eukaryotic|refseq_viral|NC_075702.1|100246|+,TACGGTATTTCT...,-1
human|human|chr1|217010|+,ATCGCTAGATAA...,-1
bacteria|uhgg|MGYG000002155_4|86042|-,TCCCGCTGCGAC...,-1
```

The ID is `category|source|record|start|strand`. `start` is the 0-based position in the
source record, and `strand` is `-` if the window was reverse-complemented. These files load
unchanged with `viral_jepa.DNADataset(pd.read_csv(path, header=None), augment=True)`. This
was checked: 300×4 one-hot tensors, label `-1`.

**`manifest.json`** records:
- all arguments;
- the per-category targets and achieved counts;
- records and Mbp per source;
- Anelloviridae exclusions per source, per-species cap counts, RefSeq genomes dropped as
  duplicates of NCBI Virus, and genomes whose species came from the header name;
- windows rejected for N or short records;
- windows removed by each screen per category;
- screen definitions (k, threshold, number of k-mers);
- source URLs and run time.

Cite the numbers from it when reporting the corpus.

---

## 6. Resource estimates

| Item | Approx. |
|---|---|
| Download | ~6–8 GB raw (CHM13 1 GB, UHGV 1 GB, RefSeq 0.2 GB, UHGG ~3–4 GB for 1000 genomes, NCBI Virus variable) |
| Output, 4M windows | ~1.25 GB CSV |
| RAM | ~6–10 GB peak (window pool ≈ 1.3 × 4M × 300 B, evaluation k-mer arrays ≈ 0.2 GB, one human chromosome decoded at a time). Request 16 GB. |
| Time | Local test with RefSeq viral, 4 UHGG genomes and small stand-ins, 100k windows: **~40 s**, about half of it building the evaluation k-mer set. A full build is dominated by streaming UHGV and CHM13 twice: expect tens of minutes. The manifest records the actual time. |

---

## 8. Known limitations
- **RefSeq viral contains more than human viruses.** It also includes plant, insect and other
  viruses. They still teach "viral-ness", but you can raise the human-relevant share by
  giving NCBI Virus more weight.
- **UHGV is labelled `virus_prokaryotic` as a whole.** It is overwhelmingly phage, but includes
  some eukaryotic gut viruses.
- **The species cap is taxonomic.** Species with many strains are capped, but strains of one
  species are still near-identical, and unclassified records ("Polyomavirus sp.") each form
  their own species. Clustering genomes with mmseqs2 at 95% identity would be cleaner, but
  adds a dependency.
- **Bacterial coverage is gut-centred (UHGG).** ViraMiner also has skin, cervical and other
  samples. Skin and vaginal reference genomes (e.g. from HMP) can be added as another
  `bacteria` source.
- **ViraMiner's train and test splits are themselves very similar at ~90% identity.** When
  `fullset_train` was screened in strict mode against val/test/anell_test (in an earlier
  test run), about 65% of its windows matched, of which about 12% are exact overlaps with
  `anell_test`. The same organisms recur
  across samples. That's a property of the benchmark, not of this corpus, but keep it in mind
  when interpreting any ViraMiner result.

---
