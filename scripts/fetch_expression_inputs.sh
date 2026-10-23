#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
raw_dir="$project_dir/data/expression"
python_bin="${RETROEVENT_PYTHON:-python3}"
mkdir -p "$raw_dir"/{encode,fantom5,gse137490,reference,riboseq}
mkdir -p "$project_dir/data/event/annotation"

download() {
    local url="$1"
    local output="$2"
    if [[ -s "$output" ]]; then
        return
    fi
    curl --fail --location --retry 4 --retry-delay 2 "$url" --output "$output"
}

download_encode() {
    local accession="$1"
    local extension="$2"
    download "https://www.encodeproject.org/files/${accession}/@@download/${accession}.${extension}" "$raw_dir/encode/${accession}.${extension}"
}

verify_sha256() {
    local expected="$1"
    local file="$2"
    local observed
    observed="$(sha256sum "$file" | awk '{print $1}')"
    if [[ "$observed" != "$expected" ]]; then
        printf 'SHA-256 mismatch for %s\nexpected: %s\nobserved: %s\n' "$file" "$expected" "$observed" >&2
        exit 1
    fi
}

download "https://rest.ensembl.org/lookup/id/ENSMUSG00000086727?expand=1;content-type=application/json" "$project_dir/data/event/annotation/ENSMUSG00000086727.lookup.json"
download "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000067.6&seq_start=31090000&seq_stop=31230000&rettype=fasta&strand=1" "$raw_dir/reference/mm10_NC_000067.6_31090000-31230000.fasta"

for accession in ENCFF183XSN ENCFF448JGR ENCFF957IRE ENCFF802BII ENCFF926TOA ENCFF392NZM; do
    download_encode "$accession" tsv
done
verify_sha256 3ae27665a86ee87d1e94199137b8f11fdcbfa14590b7ed8889d70a4b091ecad0 "$raw_dir/encode/ENCFF183XSN.tsv"
verify_sha256 f46c3cdee671ac088a0025e5ec2c94989cd51cde35a18d1bf5f0ed77739b396d "$raw_dir/encode/ENCFF448JGR.tsv"
verify_sha256 ca35910ccc7555c319dfd4b9a56703ff01d447ed90a9b51a5a742be07ef9df7c "$raw_dir/encode/ENCFF957IRE.tsv"
verify_sha256 a6ecae5243c49cd50702b322ba005948d722bbe57eea22ffab6484f3cb1a5a33 "$raw_dir/encode/ENCFF802BII.tsv"
verify_sha256 b83e30f915e3016df607d263e78b4039b2477d3d5229ad16f850917b41ed567e "$raw_dir/encode/ENCFF926TOA.tsv"
verify_sha256 055b587bda6ca4a78aacf175cbd6b895c1e52cde89db810180b8081ff9de9894 "$raw_dir/encode/ENCFF392NZM.tsv"

fantom_full="$raw_dir/fantom5/mm9.cage_peak_phase1and2combined_tpm.osc.txt.gz"
download "https://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/mm9.cage_peak_phase1and2combined_tpm.osc.txt.gz" "$fantom_full"
verify_sha256 db9321da9befef691bde1e8f831ca2945a485d1501ec591a390a5257ec6cc2b9 "$fantom_full"
gzip -dc "$fantom_full" | awk -F '\t' 'NR == 1 || $1 == "chr1:31279449..31279465,-" || $1 == "chr1:31279472..31279488,-" || $1 == "chr1:31279489..31279522,-" || $1 == "chr1:31279627..31279637,+" || $1 == "chr1:31279653..31279666,+"' > "$raw_dir/fantom5/mm9_Alyreffm1_Dnaaf6rt_CAGE_peaks.tsv"
verify_sha256 9e21472d94af75ef363c365d2f8d173930f95ef91ac8100c062f6705a2a1a2df "$raw_dir/fantom5/mm9_Alyreffm1_Dnaaf6rt_CAGE_peaks.tsv"
rm "$fantom_full"

geo_base="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137490/suppl"
download "$geo_base/GSE137490_TestisPacBio_FLNC.fasta.gz" "$raw_dir/gse137490/GSE137490_TestisPacBio_FLNC.fasta.gz"
download "$geo_base/GSE137490_TranscriptomeAssembly.MouseTestis.mm10.gtf.gz" "$raw_dir/gse137490/GSE137490_TranscriptomeAssembly.MouseTestis.mm10.gtf.gz"
download "$geo_base/GSE137490_Summary.MousePacBio.ONT.xlsx" "$raw_dir/gse137490/GSE137490_Summary.MousePacBio.ONT.xlsx"
verify_sha256 0b8fb3adb0fab561ab2906cd410743a3369e16cfa24e6176904439f41d0258fa "$raw_dir/gse137490/GSE137490_TestisPacBio_FLNC.fasta.gz"
verify_sha256 afa7ac45fc5d6354a0da0a75cc2b5704ad9171459d7fff4c7f201da70c4fdda0 "$raw_dir/gse137490/GSE137490_TranscriptomeAssembly.MouseTestis.mm10.gtf.gz"
verify_sha256 e09f9c6a05e37cffa2d0d27e29d1a27eab8f866d8752cd7e3aaa64b20a946aae "$raw_dir/gse137490/GSE137490_Summary.MousePacBio.ONT.xlsx"

ribo_base="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4160nnn"
download "$ribo_base/GSM4160756/suppl/GSM4160756_Riboseq.mouse.testis.wt.adult.rep1.gz.x_rRNA.x_hairpin.mm10v1.unique.%2Bjxn.bed13.txt.gz" "$raw_dir/riboseq/GSM4160756.bed13.txt.gz"
download "$ribo_base/GSM4160757/suppl/GSM4160757_Riboseq.mouse.testis.wt.adult.rep2.gz.x_rRNA.x_hairpin.mm10v1.unique.%2Bjxn.bed13.txt.gz" "$raw_dir/riboseq/GSM4160757.bed13.txt.gz"
verify_sha256 c7e4d231995ec5831716a779afd4da333303eeaa419bc0895f8b7d366548311c "$raw_dir/riboseq/GSM4160756.bed13.txt.gz"
verify_sha256 72349708c9a84129c9e57d64e9b3b8f741da33c0efc0834f59ea0b623897338e "$raw_dir/riboseq/GSM4160757.bed13.txt.gz"

for accession in ENCFF001JPX ENCFF001JPZ; do
    output="$project_dir/results/expression/${accession}.junction_kmer_hits.tsv"
    temporary_output="${output}.tmp"
    curl --fail --location --retry 4 "https://www.encodeproject.org/files/${accession}/@@download/${accession}.fastq.gz" \
        | gzip -dc \
        | "$python_bin" "$project_dir/scripts/scan_junction_kmers.py" --fastq - --source "$accession" --output "$temporary_output"
    mv "$temporary_output" "$output"
done

echo "Public inputs are ready. Run: python scripts/03_analyse_expression_full.py"

