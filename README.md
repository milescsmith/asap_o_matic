# asap-o-matic

[![codecov](https://codecov.io/gh/milescsmith/asap_o_matic/graph/badge.svg?token=hu7JytF189)](https://codecov.io/gh/milescsmith/asap_o_matic)
[![PyPI version](https://badge.fury.io/py/asap-o-matic.svg)](https://badge.fury.io/py/asap-o-matic)
[![](https://img.shields.io/pypi/l/asap-o-matic.svg)](https://github.com/milescsmith/asap_o_matic/blob/main/LICENSE)
[![image](https://img.shields.io/pypi/pyversions/asap-o-matic.svg)](https://pypi.python.org/pypi/asap-o-matic)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

![Alt](https://repobeats.axiom.co/api/embed/f1b4bb2f7e8274fa08b39844f2ef37bfe02c0daf.svg "Repobeats analytics image")

asap-o-matic provides the ability to process [ASAP-seq](https://www.nature.com/articles/s41587-021-00927-2) FASTQs for
downstream processing and counting of the antibody-dependent reads using [Salmon Alevin](https://salmon.readthedocs.io/en/latest/alevin.html).

A heavily modified version of [asap_to_kite](https://github.com/caleblareau/asap_to_kite).

## About

ASAP-seq uses a few tricks to bridge the oligo sequences attached to CITE-seq/Total-seq antibodies with the oligo tails
on the beads of 10x Genomics scATAC-seq kits; however, the reads produced don't match anything that Cellranger understands
how to count. asap-o-matic reformats those reads so that they appear like those coming from the feature library of a
10x Genomics scRNA-seq library.

## Installation

The easiest way is to run via [uv](https://github.com/astral-sh/uv):

```console
uv tool install asap-o-matic
```

Alteratively, it can be installed using `pip` 

```console
pip install asap-o-matic
```

or `uv`:
```console
uv pip install asap_o_matic
```

## Requirements

* Python >= 3.11
  - currently, asap-o-matic is tested against 3.11-3.14
* [Rust](https://rust-lang.org/)
* R1/R2/I1/I2 files output by `bcl-convert/bcl2fastq` or the R1/R2/R3/I3 produced by `cellranger 
mkfastq`

## Usage:

```console
asap-o-matic [OPTIONS] COMMAND [ARGS]...
```

### Options :
* `-f, --fastqs DIRECTORY`: Path of folder created by mkfastq or bcl2fastq; can be comma separated that will be collapsed into one output  [required]
* `-s, --sample TEXT`: Prefix of the filenames of FASTQs to select; can be comma separated that will be collapsed into one output  [required]
* `-o, --id TEXT`: A unique run id, used to name output.  [required]
* `-a, --fastq_source [cellranger|bcl-convert]`: Name of the program used to convert bcls to FASTQs. Cellranger mkfastq creates R1, R2, R3, and I3 files while bcl-convert creates R1, I1, R2, I2 files.  [default: cellranger]
* `-d, --outdir DIRECTORY`: Directory to save files to.  If none is give, save in the directory from which the script was called.
* `-c, --cores INTEGER`: Number of cores to use for parallel processing.  [default: 18]
* `-r, --rc-R2 / -R, --no-rc-R2`: Should the reverse complement of R2 be used? Pass &#x27;--rc-R2&#x27; if the reads were generated on a NextSeq or v1.0 chemistry NovaSeq.  [default: no-rc-R2]
* `-j, --conjugation [TotalSeqA|TotalSeqB]`: String specifying antibody conjugation; either TotalSeqA or TotalSeqB  [default: TotalSeqA]
* `--debug`: Print extra information for debugging.
* `--save_log`: Save the log to a file
* `--version`: Print version number.
* `--help`: Show this message and exit.


### Example usage:

Assuming we have FASTQs from bcl-convert in the folder `/path/to/fastq/folder/sample_1` that are named:
* sample_1_prot_S11_L004_R1_001.fastq.gz
* sample_1_prot_S11_L004_R2_001.fastq.gz
* sample_1_prot_S11_L004_I1_001.fastq.gz
* sample_1_prot_S11_L004_I2_001.fastq.gz

```console
asap-o-matic \
    --fastqs /path/to/fastq/folder \
    --sample sample_1_prot \
    --id sample_1_reformatted \
    --conjugation TotalSeqB \
    --outdir /path/to/output/sample_1 \
    --cores 24 \
    --no-rc-R2
```

The resulting reformatted reads will be output as:
* /path/to/output/sample_1/sample_1_reformatted_R1.fastq.gz
* /path/to/output/sample_1/sample_1_reformatted_R2.fastq.gz

