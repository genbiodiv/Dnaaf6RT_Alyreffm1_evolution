# Data provenance

Every retained biological input is a snapshot of an accessioned public
resource. `MANIFEST.sha256` records the exact bytes distributed here.

## Event reconstruction

| Local material | Public source | Assembly/accession | Use |
|---|---|---|---|
| chromosome-1 recipient FASTA (shared with the synteny module) | NCBI EFetch/RefSeq | GRCm39, `NC_000067.7:31,155,362–31,303,806` | insertion, TSD and recipient sequence |
| chromosome-X source FASTA | NCBI EFetch/RefSeq | GRCm39, `NC_000086.8:138,955,000–139,015,000` | source-homologous 5′ block |
| parent RNA/CDS/protein | NCBI Datasets/RefSeq | `NM_177921.3`, `NP_808589.1` | parent–retrogene comparisons |
| retrogene CDS/protein | NCBI Datasets/RefSeq | `NM_029062.2`, `NP_083338.1` | processed-copy comparisons |
| expanded gene lookup | Ensembl REST | `ENSMUSG00000086727` | current `Alyreffm1` models |
| BED and unmapped records | UCSC liftOver | `mm39ToRn7.over.chain.gz` | mouse–rat synteny audit |

The UCSC chain description is available at
<https://genome.ucsc.edu/goldenpath/help/chain.html>; the chain used was
downloaded from the UCSC mm39 liftOver directory.

## Multispecies synteny

`data/synteny/species.tsv` is the authoritative table of the nine species,
taxa, `Lgsn` Gene IDs, assemblies, sequence accessions, strand and interval
coordinates. The included FASTA files were retrieved through NCBI EFetch.
`data/synteny/metadata/locus_provenance.tsv` records the retrieval URL, time
and SHA-256 hash for every locus. The set comprises:

- *Mus musculus*, *M. caroli* and *M. pahari* (Murini);
- *Mastomys coucha* (Praomyini);
- *Apodemus sylvaticus* (Apodemini);
- *Arvicanthis niloticus* and *Grammomys surdaster* (Arvicanthini); and
- *Rattus rattus* and *R. norvegicus* (Rattini).

## Expression

### Included compact inputs

| Input | Public source | Reference |
|---|---|---|
| shared expanded Ensembl gene JSON (`data/event/annotation/`) | Ensembl REST | `ENSMUSG00000086727` (contains the transcript models) |
| five-peak CAGE subset | FANTOM5 phase 1+2 mouse CAGE matrix | mm9 focal interval |
| local reference FASTA | NCBI EFetch | mm10/GRCm38 `NC_000067.6:31,090,000–31,230,000` |
| transcript assembly GTF | GEO | GSE137490 |

### Large inputs retrieved by the optional workflow

| Evidence | Accessions/files | Observed SHA-256 where applicable |
|---|---|---|
| RSEM tables | ENCFF183XSN, ENCFF448JGR, ENCFF957IRE, ENCFF802BII, ENCFF926TOA, ENCFF392NZM | hashes are checked against the values below |
| exact junction FASTQ scans | ENCFF001JPX, ENCFF001JPZ | positive-read tables are retained in `results/expression/` |
| PacBio FLNC reads | `GSE137490_TestisPacBio_FLNC.fasta.gz` | `0b8fb3adb0fab561ab2906cd410743a3369e16cfa24e6176904439f41d0258fa` |
| PacBio transcript assembly | `GSE137490_TranscriptomeAssembly.MouseTestis.mm10.gtf.gz` | `afa7ac45fc5d6354a0da0a75cc2b5704ad9171459d7fff4c7f201da70c4fdda0` |
| Ribo-seq replicate 1 | GSM4160756 / SRX7141430 | `c7e4d231995ec5831716a779afd4da333303eeaa419bc0895f8b7d366548311c` |
| Ribo-seq replicate 2 | GSM4160757 / SRX7141431 | `72349708c9a84129c9e57d64e9b3b8f741da33c0efc0834f59ea0b623897338e` |

Observed RSEM-table SHA-256 values, in the order listed above:

```text
ENCFF183XSN  3ae27665a86ee87d1e94199137b8f11fdcbfa14590b7ed8889d70a4b091ecad0
ENCFF448JGR  f46c3cdee671ac088a0025e5ec2c94989cd51cde35a18d1bf5f0ed77739b396d
ENCFF957IRE  ca35910ccc7555c319dfd4b9a56703ff01d447ed90a9b51a5a742be07ef9df7c
ENCFF802BII  a6ecae5243c49cd50702b322ba005948d722bbe57eea22ffab6484f3cb1a5a33
ENCFF926TOA  b83e30f915e3016df607d263e78b4039b2477d3d5229ad16f850917b41ed567e
ENCFF392NZM  055b587bda6ca4a78aacf175cbd6b895c1e52cde89db810180b8081ff9de9894
```

The complete FANTOM5 matrix observed during the analysis had SHA-256
`db9321da9befef691bde1e8f831ca2945a485d1501ec591a390a5257ec6cc2b9`.
The retained five-peak subset had SHA-256
`9e21472d94af75ef363c365d2f8d173930f95ef91ac8100c062f6705a2a1a2df`.

## Retrieval behavior

Download helpers stop on HTTP errors. They do not synthesize, impute or mock a
missing source file. Because public resources can update metadata or retire
URLs, the distributed snapshots and their hashes define the exact inputs used
for this release.
