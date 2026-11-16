import argparse
import csv
import gzip
import io
import json
import re
import shutil
import subprocess
import sys
import time
import urllib.request
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np


WINDOW = 300

REFSEQ_VIRAL_INDEX = "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"
CHM13_URL = (
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/"
    "assemblies/analysis_set/chm13v2.0.fa.gz"
)
UHGV_RECORD = "https://zenodo.org/api/records/17402089/files/{name}/content"
UHGV_FILES = ("votus_full.fna.zst", "votus_metadata_extended.tsv.zst")
UHGG_ROOT = "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/"

# Taxa of eukaryotic viruses commonly found in human clinical metagenomes. Mostly families;
# Retroviridae is split into its human species because the whole-family package (tens of
# thousands of HIV-1 genomes) comes back from NCBI as an invalid zip.
NCBI_VIRUS_TAXA = (
    "Papillomaviridae",
    "Polyomaviridae",
    "Anelloviridae",
    "Orthoherpesviridae",
    "Adenoviridae",
    "Parvoviridae",
    "Circoviridae",
    "Genomoviridae",
    "Hepadnaviridae",
    "Poxviridae",
    "Picornaviridae",
    "Caliciviridae",
    "Astroviridae",
    "Flaviviridae",
    "Human immunodeficiency virus 1",
    "Human immunodeficiency virus 2",
    "Primate T-lymphotropic virus 1",
    "Primate T-lymphotropic virus 2",
)

CATEGORIES = ("virus_eukaryotic", "virus_prokaryotic", "human", "bacteria", "viraminer")
DEFAULT_MIX = "virus_eukaryotic=0.30,virus_prokaryotic=0.15,human=0.25,bacteria=0.20,viraminer=0.10"

ANELLO_PATTERN = re.compile(
    r"anellovir|torque[ _-]?teno|torquetenovirus|\bttv\b|\bttmv|\bttmdv|gyrovirus|chicken anemia",
    re.IGNORECASE,
)
PHAGE_PATTERN = re.compile(r"phage", re.IGNORECASE)
SPECIES_CUT = re.compile(
    r"\s+(isolate|strain|clone|segment|genomic|complete|partial|from|genotype|variant)\b|,",
    re.IGNORECASE,
)

SCREEN_PRESETS = {
    "none": None,
    # Near-exact overlap: >=10% of the window's 31-mers appear in an evaluation sequence.
    "standard": (31, 0.10),
    # ~90% identity: 15-mers survive 10% divergence often enough to exceed this fraction,
    # while random 15-mer matches stay around 1-2% of a window.
    "strict": (15, 0.10),
}

CODE = np.full(256, 4, dtype=np.uint8)
for index, base in enumerate(b"ACGT"):
    CODE[base] = index
    CODE[base + 32] = index
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


# --------------------------------------------------------------------------- I/O


def log(message):
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    if path.suffix == ".zst":
        # UHGV files use long-distance matching (512 MB window), so the window limit is raised.
        try:
            import zstandard

            decompressor = zstandard.ZstdDecompressor(max_window_size=2**31)
            return io.TextIOWrapper(decompressor.stream_reader(open(path, "rb")))
        except ImportError:
            if shutil.which("zstd") is None:
                raise RuntimeError(f"{path} needs the 'zstandard' module or the 'zstd' binary")
            command = ["zstd", "-dc", "--long=31", str(path)]
            return subprocess.Popen(command, stdout=subprocess.PIPE, text=True).stdout
    return open(path)


def read_fasta(path):
    header, chunks = None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def download(url, destination, retries=3):
    destination = Path(destination)
    if destination.exists() and destination.stat().st_size > 0:
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "viraljepa-corpus/1.0"})
            with urllib.request.urlopen(request, timeout=120) as response, open(partial, "wb") as out:
                shutil.copyfileobj(response, out, length=1 << 20)
            partial.rename(destination)
            return destination
        except Exception as error:  # noqa: BLE001
            log(f"  download failed ({attempt}/{retries}) {url}: {error}")
            time.sleep(5 * attempt)
    raise RuntimeError(f"could not download {url}")


# --------------------------------------------------------------------- download


def download_refseq_viral(raw_dir, args):
    index = urllib.request.urlopen(REFSEQ_VIRAL_INDEX, timeout=60).read().decode()
    names = sorted(set(re.findall(r"viral\.\d+\.\d+\.genomic\.fna\.gz", index)))
    for name in names:
        log(f"  RefSeq viral: {name}")
        download(REFSEQ_VIRAL_INDEX + name, raw_dir / "refseq_viral" / name)


def download_ncbi_virus(raw_dir, args):
    if shutil.which("datasets") is None:
        log(
            "  NCBI Virus skipped: the NCBI 'datasets' CLI is not on PATH "
            "(https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)"
        )
        return
    target = raw_dir / "ncbi_virus"
    target.mkdir(parents=True, exist_ok=True)
    for taxon in args.ncbi_taxa.split(","):
        name = taxon.replace(" ", "_")
        fasta, report = target / f"{name}.fna", target / f"{name}.jsonl"
        # The data report carries the species taxonomy used by the per-species cap.
        if fasta.exists() and report.exists():
            continue
        archive = target / f"{name}.zip"
        log(f"  NCBI Virus: {taxon} (host human, complete genomes)")
        command = [
            "datasets", "download", "virus", "genome", "taxon", taxon,
            "--host", "human", "--complete-only", "--no-progressbar", "--filename", str(archive),
        ]
        for attempt in range(1, 4):
            if subprocess.run(command).returncode == 0 and zipfile.is_zipfile(archive):
                break
            log(f"  NCBI Virus: {taxon} failed ({attempt}/3)")
            time.sleep(30 * attempt)
        else:
            archive.unlink(missing_ok=True)
            continue
        with zipfile.ZipFile(archive) as bundle:
            members = {
                destination: next((n for n in bundle.namelist() if n.endswith(suffix)), None)
                for destination, suffix in ((fasta, "genomic.fna"), (report, "data_report.jsonl"))
            }
            if None in members.values():
                log(f"  NCBI Virus: {taxon} archive lacks genomic.fna or data_report.jsonl")
                continue
            for destination, member in members.items():
                with bundle.open(member) as source, open(destination, "wb") as out:
                    shutil.copyfileobj(source, out)
        archive.unlink()


def download_uhgv(raw_dir, args):
    for name in UHGV_FILES:
        log(f"  UHGV: {name}")
        download(UHGV_RECORD.format(name=name), raw_dir / "uhgv" / name)


def download_uhgg(raw_dir, args):
    target = raw_dir / "uhgg"
    metadata = download(UHGG_ROOT + "genomes-all_metadata.tsv", target / "genomes-all_metadata.tsv")
    representatives = {}
    with open(metadata) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["Genome"] != row["Species_rep"]:
                continue
            if float(row["Completeness"]) < 90 or float(row["Contamination"]) > 5:
                continue
            phylum = next((p for p in row["Lineage"].split(";") if p.startswith("p__")), "p__unknown")
            representatives[row["Genome"]] = phylum

    # Round-robin over phyla so that the subset is not dominated by Firmicutes/Bacteroidota.
    rng = np.random.default_rng(args.seed)
    by_phylum = defaultdict(list)
    for genome, phylum in representatives.items():
        by_phylum[phylum].append(genome)
    for genomes in by_phylum.values():
        rng.shuffle(genomes)
    chosen = []
    while len(chosen) < min(args.uhgg_genomes, len(representatives)):
        for genomes in by_phylum.values():
            if genomes and len(chosen) < args.uhgg_genomes:
                chosen.append(genomes.pop())

    (target / "selected_genomes.txt").write_text("\n".join(chosen) + "\n")
    log(f"  UHGG: {len(chosen)} species representatives from {len(by_phylum)} phyla")

    def fetch(genome):
        url = f"{UHGG_ROOT}species_catalogue/{genome[:11]}/{genome}/genome/{genome}.fna"
        download(url, target / "genomes" / f"{genome}.fna")

    with ThreadPoolExecutor(args.workers) as pool:
        list(pool.map(fetch, chosen))


def download_human(raw_dir, args):
    log("  Human: T2T-CHM13v2.0 analysis set")
    download(CHM13_URL, raw_dir / "human" / "chm13v2.0.fa.gz")


DOWNLOADERS = {
    "refseq_viral": download_refseq_viral,
    "ncbi_virus": download_ncbi_virus,
    "uhgv": download_uhgv,
    "uhgg": download_uhgg,
    "human": download_human,
}


def run_download(args):
    raw_dir = Path(args.raw_dir)
    for source in args.sources.split(","):
        log(f"Downloading {source}")
        DOWNLOADERS[source](raw_dir, args)
    log("Downloads complete")


# ------------------------------------------------------------- record indexing


class Record:
    __slots__ = ("source", "category", "path", "record_id", "length")

    def __init__(self, source, category, path, record_id, length):
        self.source = source
        self.category = category
        self.path = path
        self.record_id = record_id
        self.length = length


def record_id_of(header):
    return header.split()[0].replace(",", ";")


def species_key(header):
    description = header.split(maxsplit=1)[1] if " " in header else header
    return SPECIES_CUT.split(description, maxsplit=1)[0].strip().lower()


def lineage_species(virus):
    """ICTV species of an NCBI virus data report: the lineage entry right after the genus.

    Old common names ("Monkeypox virus", "Human herpesvirus 5") sit below the binomial species
    ("Orthopoxvirus monkeypox", "Cytomegalovirus humanbeta5"), so keying on it merges them.
    """
    lineage = virus.get("lineage", [])
    # The last one-word "...virus" name is the genus (or a subgenus below it).
    for index in reversed(range(len(lineage) - 1)):
        name = lineage[index]["name"]
        if " " not in name and name.endswith("virus"):
            return f"taxid:{lineage[index + 1]['taxId']}"
    return f"taxid:{virus.get('taxId')}"


_ncbi_species = None


def ncbi_virus_species(raw_dir):
    """Accession -> species key for every genome in the NCBI Virus data reports (cached)."""
    global _ncbi_species
    if _ncbi_species is None:
        _ncbi_species = {}
        for path in sorted((raw_dir / "ncbi_virus").glob("*.jsonl")):
            with open(path) as handle:
                for line in handle:
                    report = json.loads(line)
                    _ncbi_species[report["accession"]] = lineage_species(report.get("virus", {}))
    return _ncbi_species


def index_refseq_viral(raw_dir, args, stats):
    # RefSeq genomes also present in the NCBI Virus download are counted there, under the
    # shared species cap, instead of a second time here.
    in_ncbi = ncbi_virus_species(raw_dir) if "ncbi_virus" in args.sources.split(",") else {}
    for path in sorted((raw_dir / "refseq_viral").glob("viral.*.genomic.fna.gz")):
        for header, sequence in read_fasta(path):
            if args.exclude_anelloviridae and ANELLO_PATTERN.search(header):
                stats["excluded_anelloviridae"]["refseq_viral"] += 1
                continue
            if record_id_of(header) in in_ncbi:
                stats["refseq_also_in_ncbi_virus"]["refseq_viral"] += 1
                continue
            category = "virus_prokaryotic" if PHAGE_PATTERN.search(header) else "virus_eukaryotic"
            yield Record("refseq_viral", category, path, record_id_of(header), len(sequence))


def index_ncbi_virus(raw_dir, args, stats):
    per_species = Counter()
    species = ncbi_virus_species(raw_dir)
    # Strains kept per species are a seeded random draw, not the first ones in file order.
    by_species = defaultdict(list)
    for accession, key in sorted(species.items()):
        by_species[key].append(accession)
    rng = np.random.default_rng(args.seed)
    selected = {
        accession
        for accessions in by_species.values()
        for accession in rng.permutation(accessions)[: args.max_per_species]
    }
    for path in sorted((raw_dir / "ncbi_virus").glob("*.fna")):
        family = path.stem
        if args.exclude_anelloviridae and family == "Anelloviridae":
            stats["excluded_anelloviridae"]["ncbi_virus_family_files"] += 1
            continue
        for header, sequence in read_fasta(path):
            if args.exclude_anelloviridae and ANELLO_PATTERN.search(header):
                stats["excluded_anelloviridae"]["ncbi_virus"] += 1
                continue
            record_id = record_id_of(header)
            if record_id in species:
                capped = record_id not in selected
            else:
                stats["species_from_header_name"]["ncbi_virus"] += 1
                key = species_key(header)
                per_species[key] += 1
                capped = per_species[key] > args.max_per_species
            if capped:
                stats["capped_per_species"]["ncbi_virus"] += 1
                continue
            yield Record("ncbi_virus", "virus_eukaryotic", path, record_id_of(header), len(sequence))


def index_uhgv(raw_dir, args, stats):
    directory = raw_dir / "uhgv"
    allowed_quality = set(args.uhgv_quality.split(","))
    keep = set()
    with open_text(directory / "votus_metadata_extended.tsv.zst") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["checkv_quality"] not in allowed_quality:
                continue
            taxonomy = f"{row['ictv_taxonomy']};{row['family_name']}"
            if args.exclude_anelloviridae and "Anelloviridae" in taxonomy:
                stats["excluded_anelloviridae"]["uhgv"] += 1
                continue
            keep.add(row["uhgv_genome"])
    log(f"  UHGV: {len(keep)} vOTU representatives pass quality {sorted(allowed_quality)}")
    path = directory / "votus_full.fna.zst"
    for header, sequence in read_fasta(path):
        record_id = record_id_of(header)
        if record_id in keep:
            yield Record("uhgv", "virus_prokaryotic", path, record_id, len(sequence))


def index_uhgg(raw_dir, args, stats):
    for path in sorted((raw_dir / "uhgg" / "genomes").glob("*.fna")):
        for header, sequence in read_fasta(path):
            yield Record("uhgg", "bacteria", path, record_id_of(header), len(sequence))


def index_human(raw_dir, args, stats):
    for path in sorted((raw_dir / "human").glob("*.fa.gz")):
        for header, sequence in read_fasta(path):
            yield Record("human", "human", path, record_id_of(header), len(sequence))


INDEXERS = {
    "refseq_viral": index_refseq_viral,
    "ncbi_virus": index_ncbi_virus,
    "uhgv": index_uhgv,
    "uhgg": index_uhgg,
    "human": index_human,
}

SOURCE_DIRS = {
    "refseq_viral": "refseq_viral",
    "ncbi_virus": "ncbi_virus",
    "uhgv": "uhgv",
    "uhgg": "uhgg/genomes",
    "human": "human",
}


# ------------------------------------------------------------ k-mer screening


def kmer_codes(codes, k):
    """Canonical 2-bit k-mers for an (n, length) array of base codes 0..3."""
    count = codes.shape[1] - k + 1
    forward = np.zeros((codes.shape[0], count), dtype=np.uint64)
    reverse = np.zeros_like(forward)
    complement = (3 - codes[:, ::-1]).astype(np.uint64)
    codes = codes.astype(np.uint64)
    for offset in range(k):
        shift = np.uint64(2 * (k - 1 - offset))
        forward |= codes[:, offset : offset + count] << shift
        reverse |= complement[:, offset : offset + count] << shift
    return np.minimum(forward, reverse[:, ::-1])


def sequences_to_codes(sequences):
    sequences = [s for s in sequences if len(s) == WINDOW]
    raw = np.frombuffer("".join(sequences).encode(), dtype=np.uint8).reshape(len(sequences), WINDOW)
    codes = CODE[raw]
    return codes[(codes < 4).all(1)]


class KmerScreen:
    def __init__(self, name, sequences, k, min_fraction):
        self.name = name
        self.k = k
        self.min_fraction = min_fraction
        codes = sequences_to_codes(sequences)
        kmers = [np.unique(kmer_codes(codes[i : i + 20000], k)) for i in range(0, len(codes), 20000)]
        self.kmers = np.unique(np.concatenate(kmers))
        log(f"  screen '{name}': k={k}, min shared fraction={min_fraction}, {len(self.kmers):,} k-mers")

    def hits(self, windows, ids):
        hit = np.zeros(len(windows), dtype=bool)
        for start in range(0, len(windows), 50000):
            queries = kmer_codes(windows[start : start + 50000], self.k)
            position = np.searchsorted(self.kmers, queries).clip(max=len(self.kmers) - 1)
            hit[start : start + 50000] = (self.kmers[position] == queries).mean(1) >= self.min_fraction
        return hit

    def describe(self):
        return {"name": self.name, "type": "kmer", "k": self.k, "min_fraction": self.min_fraction,
                "kmers": int(len(self.kmers))}


# ------------------------------------------------------------ BLAST screening


BLAST_FIELDS = "qseqid sseqid pident evalue sstart send"


def find_blast_tool(name, folder):
    path = Path(folder) / name if folder else shutil.which(name)
    if path is None or not Path(path).exists():
        raise RuntimeError(f"{name} not found: install BLAST+ or pass --blast-bin")
    return str(path)


def sequence_number(token):
    """'w12', 't12', 'lcl|t12' or 'gb|t12|' -> 12. makeblastdb reads ids such as t12345 as GenBank accessions."""
    return int(re.search(r"(\d+)\D*$", token).group(1))


def merge_intervals(intervals):
    """Union of closed (start, end) intervals, as a sorted list of disjoint closed intervals."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged, (start, end) = [], intervals[0]
    for a, b in intervals[1:]:
        if a > end + 1:
            merged.append((start, end))
            start, end = a, b
        else:
            end = max(end, b)
    merged.append((start, end))
    return merged


def covered_length(intervals):
    """Number of positions covered by a list of closed (start, end) intervals."""
    return sum(b - a + 1 for a, b in merge_intervals(intervals))


MASK_BACKGROUNDS = {"human": ("human", "*.fa.gz"), "uhgg": ("uhgg/genomes", "*.fna")}


def background_counts(kmers, k, raw_dir, backgrounds, fraction, seed, block=1000, batch=20000):
    """Occurrences of each k-mer of the sorted array `kmers` in a random sample of the background.

    Each background sequence is cut into `block`-bp blocks and a block is scanned with probability
    `fraction`; every k-mer that starts in a scanned block is counted, on both strands. Counts
    are therefore about `fraction` x the true copy number. Contig tails shorter than a block and
    k-mers with non-ACGT bases are not counted.
    """
    # Own generator: the mask must not shift the random stream that samples the corpus windows.
    rng = np.random.default_rng([seed, k])
    counts = np.zeros(len(kmers), dtype=np.int64)
    scanned = Counter()

    def scan(rows):
        valid = (rows < 4).astype(np.int32)
        valid = np.concatenate([np.zeros((len(rows), 1), np.int32), np.cumsum(valid, 1, dtype=np.int32)], 1)
        valid = (valid[:, k:] - valid[:, :-k]) == k
        queries = kmer_codes(np.where(rows < 4, rows, 0), k)[valid]
        position = np.searchsorted(kmers, queries).clip(max=len(kmers) - 1)
        position = position[kmers[position] == queries]
        return np.bincount(position, minlength=len(kmers))

    for name in backgrounds:
        folder, pattern = MASK_BACKGROUNDS[name]
        paths = sorted((Path(raw_dir) / folder).glob(pattern))
        if not paths:
            raise RuntimeError(f"mask background '{name}' not found in {Path(raw_dir) / folder}")
        pending, pending_rows = [], 0

        def flush():
            nonlocal pending, pending_rows, counts
            if pending:
                counts += scan(np.concatenate(pending))
            pending, pending_rows = [], 0

        for path in paths:
            for _, sequence in read_fasta(path):
                raw = CODE[np.frombuffer(sequence.encode(), dtype=np.uint8)]
                blocks = (len(raw) - (k - 1)) // block
                if blocks <= 0:
                    continue
                starts = np.flatnonzero(rng.random(blocks) < fraction) * block
                scanned[name] += len(starts) * block
                view = np.lib.stride_tricks.sliding_window_view(raw, block + k - 1)
                for first in range(0, len(starts), batch):
                    pending.append(view[starts[first : first + batch]])
                    pending_rows += len(pending[-1])
                    if pending_rows >= batch:
                        flush()
        flush()
        log(f"    mask background '{name}': {scanned[name] / 1e6:,.1f} Mbp scanned")
    return counts, dict(scanned)


class InformativeMask:
    """Stage 0 of the specificity-weighted BLAST rule: which bases of each test read are informative.

    A base is informative if at least one canonical k-mer covering it occurs fewer than max_copies
    times in the sampled background; otherwise it is generic (repeat). `cumulative[r, i]` is the
    number of informative bases among the first i bases of read r.
    """

    def __init__(self, rows, args):
        self.k = args.mask_k
        self.max_copies = args.mask_max_copies
        self.fraction = args.mask_sample_fraction
        self.backgrounds = args.mask_background.split(",")
        unknown = set(self.backgrounds) - set(MASK_BACKGROUNDS)
        if unknown:
            raise ValueError(f"unknown --mask-background {sorted(unknown)}, expected some of {sorted(MASK_BACKGROUNDS)}")
        k = self.k
        length = max(len(row[1]) for row in rows)
        codes = np.full((len(rows), length), 4, dtype=np.uint8)
        for index, row in enumerate(rows):
            codes[index, : len(row[1])] = CODE[np.frombuffer(row[1].encode(), dtype=np.uint8)]
        valid = np.lib.stride_tricks.sliding_window_view(codes < 4, k, axis=1).all(2)
        read_kmers = kmer_codes(np.where(codes < 4, codes, 0), k)
        table = np.unique(read_kmers[valid])
        log(f"  informativeness mask: {len(table):,} distinct {k}-mers in {len(rows):,} reads, background "
            f"{'+'.join(self.backgrounds)} sampled at {self.fraction:g}x")
        counts, self.scanned = background_counts(table, k, args.raw_dir, self.backgrounds, self.fraction, args.seed)
        # A k-mer with a non-ACGT base cannot be counted, so it never makes a base informative.
        specific = valid & (counts[np.searchsorted(table, read_kmers).clip(max=len(table) - 1)] < self.max_copies)
        change = np.zeros((len(rows), length + 1), dtype=np.int32)
        change[:, : specific.shape[1]] += specific
        change[:, k:] -= specific
        self.mask = np.cumsum(change, 1)[:, :length] > 0
        self.cumulative = np.concatenate([np.zeros((len(rows), 1), np.int32), np.cumsum(self.mask, 1, dtype=np.int32)], 1)
        self.informative = self.cumulative[:, -1]
        self.repeat_kmers = int((counts >= self.max_copies).sum())
        self.kmers = int(len(table))

    def overlap(self, read, intervals):
        """Informative bases of `read` inside disjoint closed 1-based intervals."""
        row = self.cumulative[read]
        return int(sum(row[b] - row[a - 1] for a, b in intervals))

    def intervals(self, read):
        edges = np.flatnonzero(np.diff(np.concatenate([[0], self.mask[read].astype(np.int8), [0]])))
        return ";".join(f"{a + 1}-{b}" for a, b in zip(edges[::2], edges[1::2]))

    def describe(self):
        return {"k": self.k, "max_copies": self.max_copies, "background": self.backgrounds,
                "sample_fraction": self.fraction, "megabases_scanned": {n: round(v / 1e6, 1) for n, v in self.scanned.items()},
                "distinct_read_kmers": self.kmers, "read_kmers_at_or_above_max_copies": self.repeat_kmers}


class BlastScreen:
    """Drop windows whose blastn alignments to an evaluation read amount to leakage.

    The evaluation reads are the database and the corpus windows are the queries. For each
    (window, read, strand), C is the union on the read of the HSPs with identity > min_identity
    and E-value <= evalue. Two rules (--blast-rule):

    coverage (corpus_v4): drop if C covers at least min_coverage of the read.

    informative (corpus_v5, specificity-weighted coverage): drop if C covers at least
    informative_bp informative bases of the read (see InformativeMask), or if the HSPs with
    identity >= dup_identity alone cover at least min_coverage of the read (near-duplicate, whatever
    the repeat content). Reads with fewer than min_informative informative bases are left out of
    the database.

    Every dropped window is listed in blast_screen_hits.csv with the read that removed it.
    """

    def __init__(self, name, rows, args):
        self.name = name
        self.categories = set(args.blast_categories.split(","))
        unknown = self.categories - set(CATEGORIES)
        if unknown:
            raise ValueError(f"unknown --blast-categories {sorted(unknown)}, expected some of {CATEGORIES}")
        self.min_identity = args.blast_identity
        self.min_coverage = args.blast_coverage
        self.evalue = args.blast_evalue
        self.word_size = args.blast_word_size
        self.jobs = args.blast_jobs
        self.chunk = args.blast_chunk
        self.blastn = find_blast_tool("blastn", args.blast_bin)
        version = subprocess.run([self.blastn, "-version"], capture_output=True, text=True, check=True)
        self.version = version.stdout.splitlines()[0]
        self.read_ids = [row[0] for row in rows]
        self.read_labels = [row[2] for row in rows]
        self.read_lengths = np.array([len(row[1]) for row in rows])
        self.rule = args.blast_rule
        self.informative_bp = args.blast_informative_bp
        self.dup_identity = args.blast_dup_identity
        self.min_informative = args.mask_min_informative
        # Reads the informative-bp rule applies to; the near-duplicate rule always applies to all reads.
        self.homology_labels = sorted(set(args.blast_homology_labels.split(",")) - {""})
        self.homology = np.array([not self.homology_labels or label in self.homology_labels for label in self.read_labels])
        self.mask = None
        in_database = np.ones(len(rows), dtype=bool)
        if self.rule == "informative":
            self.mask = InformativeMask(rows, args)
            in_database = self.mask.informative >= self.min_informative
            with open(Path(args.out_dir) / "informative_mask.csv", "w") as out:
                out.write("test_id,test_label,length,informative_bp,in_database,informative_intervals\n")
                for index, row in enumerate(rows):
                    out.write(f"{row[0]},{row[2]},{len(row[1])},{self.mask.informative[index]},"
                              f"{int(in_database[index])},{self.mask.intervals(index)}\n")
            for label in sorted(set(self.read_labels)):
                chosen = np.array([l == label for l in self.read_labels])
                log(f"    label {label}: {int(chosen.sum()):,} reads, median {int(np.median(self.mask.informative[chosen]))} "
                    f"informative bp, {int((chosen & ~in_database).sum()):,} left out of the database "
                    f"(< {self.min_informative} informative bp)")
        self.database_reads = int(in_database.sum())

        self.work = Path(args.out_dir) / "blast_screen_work"
        self.work.mkdir(parents=True, exist_ok=True)
        fasta = self.work / "evaluation.fa"
        with open(fasta, "w") as out:
            for index, row in enumerate(rows):
                if in_database[index]:
                    out.write(f">lcl|t{index}\n{row[1]}\n")
        self.db = self.work / "evaluation"
        command = [find_blast_tool("makeblastdb", args.blast_bin), "-in", str(fasta), "-dbtype", "nucl",
                   "-parse_seqids", "-out", str(self.db)]
        subprocess.run(command, check=True, capture_output=True)
        self.report = open(Path(args.out_dir) / "blast_screen_hits.csv", "w")
        if self.rule == "informative":
            self.report.write("window_id,test_id,test_label,rule,informative_bp,read_coverage,max_identity,min_evalue\n")
        else:
            self.report.write("window_id,test_id,test_label,read_coverage,max_identity,min_evalue\n")
        log(f"  screen '{name}': {self.version}, {self.database_reads:,} of {len(rows):,} reads in the database, "
            f"{self.rule_text()}, word size {self.word_size}, {self.jobs} blastn processes")

    def rule_text(self):
        hsps = f"HSPs with identity > {self.min_identity:g}% and E-value <= {self.evalue:g}"
        if self.rule == "informative":
            reads = f"test read with label {'/'.join(self.homology_labels)}" if self.homology_labels else "test read"
            return (f"drop a window if its {hsps} cover >= {self.informative_bp} informative bp of one {reads}, "
                    f"or if its HSPs with identity >= {self.dup_identity:g}% cover >= {self.min_coverage:g} of any test read")
        return f"drop a window if its {hsps} cover >= {self.min_coverage:g} of one test read"

    def search(self, windows, offset):
        """Run blastn on one chunk; return {window index: (score, report fields)} for the dropped windows.

        The reported read is the one with the highest score: read coverage under the coverage rule,
        informative bp covered under the informative rule (near-duplicates rank above all others).
        """
        fasta = self.work / f"windows_{offset}.fa"
        table = self.work / f"windows_{offset}.tsv"
        with open(fasta, "w") as out:
            for index, row in enumerate(windows):
                out.write(f">w{offset + index}\n{BASES[row].tobytes().decode()}\n")
        # -perc_identity only shrinks the output; the strict "> min_identity" test is applied below.
        command = [self.blastn, "-task", "blastn", "-word_size", str(self.word_size), "-query", str(fasta),
                   "-db", str(self.db), "-evalue", str(self.evalue), "-perc_identity", str(self.min_identity),
                   "-strand", "both", "-dust", "yes", "-max_target_seqs", "500", "-num_threads", "1",
                   "-outfmt", f"6 {BLAST_FIELDS}", "-out", str(table)]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"blastn failed on {fasta}: {result.stderr.strip()}")
        spans = defaultdict(list)
        with open(table) as handle:
            for line in handle:
                query, read, identity, evalue, start, end = line.split("\t")
                identity, evalue, start, end = float(identity), float(evalue), int(start), int(end)
                if identity > self.min_identity and evalue <= self.evalue:
                    key = (sequence_number(query), sequence_number(read), start <= end)
                    spans[key].append((min(start, end), max(start, end), identity, evalue))
        best = {}
        for (window, read, _), hsps in spans.items():
            union = merge_intervals([(a, b) for a, b, _, _ in hsps])
            coverage = sum(b - a + 1 for a, b in union) / self.read_lengths[read]
            identity, evalue = max(h[2] for h in hsps), min(h[3] for h in hsps)
            if self.rule == "coverage":
                if coverage < self.min_coverage:
                    continue
                score, fields = coverage, f"{coverage:.3f},{identity:.2f},{evalue:.3g}"
            else:
                informative = self.mask.overlap(read, union) if self.homology[read] else 0
                close = covered_length([(a, b) for a, b, i, _ in hsps if i >= self.dup_identity] or [(1, 0)])
                duplicate = close / self.read_lengths[read] >= self.min_coverage
                if informative < self.informative_bp and not duplicate:
                    continue
                reason = "informative" if informative >= self.informative_bp else "near_duplicate"
                score = (duplicate, informative, coverage)
                fields = f"{reason},{informative},{coverage:.3f},{identity:.2f},{evalue:.3g}"
            if window not in best or score > best[window][0]:
                best[window] = (score, f"{self.read_ids[read]},{self.read_labels[read]},{fields}")
        fasta.unlink()
        table.unlink()
        return best

    def hits(self, windows, ids):
        hit = np.zeros(len(windows), dtype=bool)
        if not len(windows):
            return hit
        started = time.time()
        with ThreadPoolExecutor(self.jobs) as pool:
            futures = [pool.submit(self.search, windows[s : s + self.chunk], s)
                       for s in range(0, len(windows), self.chunk)]
            for done, future in enumerate(as_completed(futures), 1):
                for window, (_, fields) in future.result().items():
                    hit[window] = True
                    self.report.write(f"{ids[window]},{fields}\n")
                if done % max(1, len(futures) // 10) == 0 or done == len(futures):
                    elapsed = time.time() - started
                    log(f"    blastn {done}/{len(futures)} chunks, {int(hit.sum()):,} windows hit, "
                        f"{elapsed / 60:.1f} min, ETA {elapsed / done * (len(futures) - done) / 60:.1f} min")
        self.report.flush()
        return hit

    def close(self):
        self.report.close()
        shutil.rmtree(self.work, ignore_errors=True)

    def describe(self):
        description = {"name": self.name, "type": "blastn", "blastn": self.version, "database": "fullset_test.csv",
                       "reads": len(self.read_ids), "reads_in_database": self.database_reads,
                       "categories": sorted(self.categories), "rule": self.rule_text(),
                       "word_size": self.word_size, "strand": "both", "dust": "yes", "max_target_seqs": 500}
        if self.mask is not None:
            description["informative_mask"] = dict(
                self.mask.describe(), min_informative_bp=self.min_informative,
                reads_left_out=len(self.read_ids) - self.database_reads,
                median_informative_bp=int(np.median(self.mask.informative)))
        return description


def read_viraminer(path):
    with open(path) as handle:
        return [(row[0], row[1], row[2]) for row in csv.reader(handle) if len(row) == 3]


def build_screens(args):
    """Return (screens for external genomes, screens for the ViraMiner train windows, BLAST screen or None)."""
    viraminer = Path(args.viraminer_dir)
    anell_test = viraminer / "test_data_anell_300_binary.csv"
    external, own = [], []
    preset = SCREEN_PRESETS[args.screen]
    if preset is not None:
        sequences = []
        for name in ("fullset_validation.csv", "fullset_test.csv", anell_test.name):
            if (viraminer / name).exists():
                sequences.extend(row[1] for row in read_viraminer(viraminer / name))
            else:
                log(f"  warning: {name} not found, not screened against")
        external.append(KmerScreen(f"evaluation_{args.screen}", sequences, *preset))
    if anell_test.exists():
        # fullset_train is fine-tuned on anyway, so it is only screened against anell_test,
        # which re-splits the same pool and would otherwise leak into that evaluation.
        rows = read_viraminer(anell_test)
        own.append(KmerScreen("anell_test_overlap", [row[1] for row in rows], *SCREEN_PRESETS["standard"]))
        if args.exclude_anelloviridae:
            positives = KmerScreen("anello_unseen_strict", [r[1] for r in rows if r[2] == "1"], *SCREEN_PRESETS["strict"])
            external.append(positives)
            own.append(positives)
    blast = None
    if args.blast_screen:
        # Runs after the k-mer screens, on the windows they keep.
        blast = BlastScreen("blast_fullset_test", read_viraminer(viraminer / "fullset_test.csv"), args)
    return external, own, blast


# ------------------------------------------------------------ window sampling


def parse_mix(text):
    mix = {}
    for item in text.split(","):
        name, value = item.split("=")
        if name not in CATEGORIES:
            raise ValueError(f"unknown category '{name}', expected one of {CATEGORIES}")
        mix[name] = float(value)
    return mix


def sample_from_sequence(sequence, count, rng, attempts=4):
    raw = CODE[np.frombuffer(sequence.encode(), dtype=np.uint8)]
    if len(raw) < WINDOW or count == 0:
        return np.empty((0, WINDOW), dtype=np.uint8), np.empty(0, dtype=np.int64)
    view = np.lib.stride_tricks.sliding_window_view(raw, WINDOW)
    windows, starts = [], []
    missing = count
    for _ in range(attempts):
        candidates = rng.integers(0, len(view), size=missing)
        sampled = view[candidates]
        valid = (sampled < 4).all(1)
        windows.append(sampled[valid])
        starts.append(candidates[valid])
        missing -= int(valid.sum())
        if missing == 0:
            break
    return np.concatenate(windows), np.concatenate(starts)


class WindowPool:
    def __init__(self):
        self.codes = defaultdict(list)
        self.ids = defaultdict(list)

    def add(self, category, codes, ids):
        if len(codes):
            self.codes[category].append(codes)
            self.ids[category].extend(ids)

    def collect(self, category):
        if not self.codes[category]:
            return np.empty((0, WINDOW), dtype=np.uint8), []
        return np.concatenate(self.codes[category]), self.ids[category]


def finalize_windows(codes, ids, category, args, rng, screens, stats):
    strand = rng.random(len(codes)) < 0.5
    codes[strand] = 3 - codes[strand, ::-1]
    ids = [f"{i}|{'-' if s else '+'}" for i, s in zip(ids, strand)]
    if args.substitution_rate > 0:
        noisy = rng.random(len(codes)) < args.noise_fraction
        mutate = (rng.random(codes.shape) < args.substitution_rate) & noisy[:, None]
        codes[mutate] = (codes[mutate] + rng.integers(1, 4, size=int(mutate.sum()))) % 4
    # Each screen only sees the windows that the previous screens kept.
    keep = np.arange(len(codes))
    for screen in screens:
        log(f"  {category}: screen '{screen.name}' on {len(keep):,} windows")
        hit = screen.hits(codes[keep], [ids[i] for i in keep])
        stats["screened_out"][f"{category}:{screen.name}"] = int(hit.sum())
        keep = keep[~hit]
    return codes[keep], [ids[i] for i in keep]


def run_build(args):
    raw_dir = Path(args.raw_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    stats = defaultdict(Counter)
    started = time.time()

    log("Indexing sources")
    records = defaultdict(list)
    for source in args.sources.split(","):
        if not (raw_dir / SOURCE_DIRS[source]).exists():
            log(f"  {source}: not downloaded, skipped")
            continue
        for record in INDEXERS[source](raw_dir, args, stats):
            if record.length >= WINDOW:
                records[record.category].append(record)
                stats["records"][source] += 1
                stats["bases"][source] += record.length
        log(f"  {source}: {stats['records'][source]:,} records, {stats['bases'][source] / 1e6:,.1f} Mbp")

    viraminer_train = Path(args.viraminer_dir) / "fullset_train.csv"
    viraminer_rows = read_viraminer(viraminer_train) if viraminer_train.exists() else []

    mix = parse_mix(args.mix)
    available = {c: w for c, w in mix.items() if w > 0 and (records[c] or (c == "viraminer" and viraminer_rows))}
    missing = sorted(set(c for c, w in mix.items() if w > 0) - set(available))
    if missing:
        log(f"  warning: no data for {missing}; mixture renormalised over {sorted(available)}")
    total = sum(available.values())
    targets = {c: int(round(args.num_windows * w / total)) for c, w in available.items()}
    log(f"Targets: {targets}")

    external_screens, own_screens, blast_screen = build_screens(args)
    pool = WindowPool()

    if "viraminer" in targets:
        # anell_test re-splits the same pool, so ~1.5% of fullset_train is screened out; oversample to compensate.
        wanted = min(len(viraminer_rows), int(targets["viraminer"] * args.oversample))
        chosen = rng.choice(len(viraminer_rows), wanted, replace=False)
        codes = sequences_to_codes([viraminer_rows[i][1] for i in chosen])
        pool.add("viraminer", codes, [f"viraminer|fullset_train|{viraminer_rows[i][0]}|0" for i in chosen])

    # Number of windows per record: multinomial over records, weighted by length; for the
    # genome-collection categories the weight is capped so huge genomes do not dominate.
    requests = defaultdict(dict)
    for category, target in targets.items():
        if category == "viraminer":
            continue
        lengths = np.array([r.length for r in records[category]], dtype=np.float64)
        if category != "human":
            lengths = np.minimum(lengths, args.max_genome_weight_bp)
        counts = rng.multinomial(int(target * args.oversample), lengths / lengths.sum())
        for record, count in zip(records[category], counts):
            if count:
                requests[record.path][record.record_id] = (category, record.source, int(count))

    log(f"Sampling windows from {sum(len(r) for r in requests.values()):,} records")
    for path, wanted in requests.items():
        for header, sequence in read_fasta(path):
            record_id = record_id_of(header)
            if record_id not in wanted:
                continue
            category, source, count = wanted[record_id]
            codes, starts = sample_from_sequence(sequence, count, rng)
            stats["windows_with_n_or_short"][category] += count - len(codes)
            pool.add(category, codes, [f"{category}|{source}|{record_id}|{s}" for s in starts])

    log("Strand flipping, optional noise, screening, trimming to targets")
    final_codes, final_ids = [], []
    for category, target in targets.items():
        codes, ids = pool.collect(category)
        screens = own_screens if category == "viraminer" else external_screens
        if blast_screen is not None and category in blast_screen.categories:
            screens = screens + [blast_screen]
        codes, ids = finalize_windows(codes, ids, category, args, rng, screens, stats)
        if len(codes) > target:
            keep = np.sort(rng.choice(len(codes), target, replace=False))
            codes, ids = codes[keep], [ids[i] for i in keep]
        stats["windows"][category] = len(codes)
        if len(codes) < target:
            log(f"  warning: {category} has {len(codes):,} windows, below target {target:,}")
        final_codes.append(codes)
        final_ids.extend(ids)
    if blast_screen is not None:
        blast_screen.close()

    codes = np.concatenate(final_codes)
    order = rng.permutation(len(codes))
    holdout = int(len(order) * args.holdout_fraction)
    if holdout:
        write_csv(out_dir / "pretrain_holdout.csv", codes, final_ids, order[:holdout])
    shards = []
    for shard, start in enumerate(range(holdout, len(order), args.shard_size)):
        path = out_dir / f"pretrain_corpus_{shard:03d}.csv"
        write_csv(path, codes, final_ids, order[start : start + args.shard_size])
        shards.append(path.name)

    manifest = {
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "arguments": {k: str(v) for k, v in vars(args).items() if k != "func"},
        "window_length": WINDOW,
        "targets": targets,
        "windows_per_category": dict(stats["windows"]),
        "holdout_windows": holdout,
        "shards": shards,
        "records_per_source": dict(stats["records"]),
        "megabases_per_source": {k: round(v / 1e6, 2) for k, v in stats["bases"].items()},
        "excluded_anelloviridae": dict(stats["excluded_anelloviridae"]),
        "capped_per_species": dict(stats["capped_per_species"]),
        "refseq_also_in_ncbi_virus": dict(stats["refseq_also_in_ncbi_virus"]),
        "species_from_header_name": dict(stats["species_from_header_name"]),
        "windows_rejected_n_or_short": dict(stats["windows_with_n_or_short"]),
        "windows_screened_out": dict(stats["screened_out"]),
        "screens": {
            "external_genomes": [s.describe() for s in external_screens],
            "viraminer_train": [s.describe() for s in own_screens],
            "blast": blast_screen.describe() if blast_screen is not None else None,
        },
        "sources": {
            "refseq_viral": REFSEQ_VIRAL_INDEX,
            "ncbi_virus": "NCBI datasets CLI, host=human, complete genomes",
            "uhgv": UHGV_RECORD.format(name="votus_full.fna.zst"),
            "uhgg": UHGG_ROOT,
            "human": CHM13_URL,
        },
        "minutes": round((time.time() - started) / 60, 1),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    log(f"Wrote {len(codes) - holdout:,} training + {holdout:,} holdout windows to {out_dir}")


def write_csv(path, codes, ids, order):
    with open(path, "w") as handle:
        for index in order:
            handle.write(f"{ids[index]},{BASES[codes[index]].tobytes().decode()},-1\n")


# ------------------------------------------------------------------------ CLI


def parse_args():
    parser = argparse.ArgumentParser(description="Build an unlabeled 300 bp pretraining corpus for ViralJEPA")
    commands = parser.add_subparsers(dest="command", required=True)

    fetch = commands.add_parser("download", help="download raw genome collections")
    fetch.add_argument("--raw-dir", required=True)
    fetch.add_argument("--sources", default="refseq_viral,ncbi_virus,uhgv,uhgg,human")
    fetch.add_argument("--ncbi-taxa", default=",".join(NCBI_VIRUS_TAXA))
    fetch.add_argument("--uhgg-genomes", type=int, default=1000)
    fetch.add_argument("--workers", type=int, default=8)
    fetch.add_argument("--seed", type=int, default=0)
    fetch.set_defaults(func=run_download)

    build = commands.add_parser("build", help="sample, screen and write 300 bp windows")
    build.add_argument("--raw-dir", required=True)
    build.add_argument("--viraminer-dir", required=True, help="folder with fullset_*.csv and *_anell_* files")
    build.add_argument("--out-dir", required=True)
    build.add_argument("--sources", default="refseq_viral,ncbi_virus,uhgv,uhgg,human")
    build.add_argument("--num-windows", type=int, default=4_000_000)
    build.add_argument("--mix", default=DEFAULT_MIX)
    build.add_argument("--screen", choices=sorted(SCREEN_PRESETS), default="standard")
    build.add_argument(
        "--keep-anelloviridae",
        dest="exclude_anelloviridae",
        action="store_false",
        help="keep anelloviruses (breaks the anell_test 'unseen family' evaluation)",
    )
    build.add_argument(
        "--blast-screen",
        action="store_true",
        help="also drop windows that blastn-align to fullset_test (needs BLAST+); see --blast-rule",
    )
    build.add_argument(
        "--blast-rule",
        choices=("coverage", "informative"),
        default="coverage",
        help="coverage (v4): HSPs cover --blast-coverage of a test read. informative (v5): HSPs cover "
        "--blast-informative-bp non-repeat bases of a test read, or the read is a near-duplicate",
    )
    build.add_argument("--blast-informative-bp", type=int, default=150,
                       help="informative rule: informative bases of the read the HSPs must cover")
    build.add_argument("--blast-dup-identity", type=float, default=95.0,
                       help="informative rule: HSPs with at least this identity (%%) covering --blast-coverage "
                       "of the read drop the window whatever the repeat content")
    build.add_argument("--blast-homology-labels", default="",
                       help="informative rule: comma list of test labels the informative-bp rule applies to "
                       "(1 = viral reads only); empty = all reads. The near-duplicate rule always uses all reads")
    build.add_argument("--mask-background", default="human",
                       help="informative rule: comma list of backgrounds for the k-mer counts (human, uhgg)")
    build.add_argument("--mask-k", type=int, default=31)
    build.add_argument("--mask-max-copies", type=int, default=50,
                       help="a k-mer seen at least this often in the sampled background is generic")
    build.add_argument("--mask-sample-fraction", type=float, default=0.1,
                       help="fraction of the background scanned; --mask-max-copies applies to the sampled counts")
    build.add_argument("--mask-min-informative", type=int, default=60,
                       help="test reads with fewer informative bases are left out of the BLAST database")
    build.add_argument("--blast-identity", type=float, default=70.0, help="HSP identity must be above this (%%)")
    build.add_argument("--blast-coverage", type=float, default=0.8, help="fraction of the test read the HSPs must cover")
    build.add_argument("--blast-evalue", type=float, default=1e-3, help="HSP E-value must be at most this")
    build.add_argument("--blast-word-size", type=int, default=7)
    build.add_argument("--blast-categories", default=",".join(CATEGORIES), help="categories the BLAST screen applies to")
    build.add_argument("--blast-jobs", type=int, default=4, help="blastn processes run in parallel")
    build.add_argument("--blast-chunk", type=int, default=2000, help="windows per blastn call")
    build.add_argument("--blast-bin", default=None, help="folder with blastn and makeblastdb (default: PATH)")
    build.add_argument("--uhgv-quality", default="Complete,High-quality")
    build.add_argument("--max-per-species", type=int, default=10)
    build.add_argument("--max-genome-weight-bp", type=int, default=50_000)
    build.add_argument("--oversample", type=float, default=1.3)
    build.add_argument("--substitution-rate", type=float, default=0.0)
    build.add_argument("--noise-fraction", type=float, default=0.5)
    build.add_argument("--holdout-fraction", type=float, default=0.01)
    build.add_argument("--shard-size", type=int, default=1_000_000)
    build.add_argument("--seed", type=int, default=0)
    build.set_defaults(func=run_build)
    return parser.parse_args()


def main():
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())
