import gzip
from collections.abc import Iterator
from enum import Enum
from functools import partial
from itertools import chain
from pathlib import Path
from typing import Annotated

import fastq as fq
import typer
from loguru import logger
from multiprocess import Pool, cpu_count
from revseq import revseq
from rich.traceback import install

from asap_o_matic import app, verbosity_level, version_callback
from asap_o_matic.logger import init_logger

install()


class Conjugation(str, Enum):
    TotalSeqA = "TotalSeqA"
    TotlaSeqB = "TotalSeqB"


DEFAULT_NUMBER_OF_THREADS = cpu_count()
DEFAULT_MAX_READS_PER_ITERATION = 10000000


# Function to verify R1/R2/R3 are present for nominated samples
def verify_sample_from_R1(list_of_R1s: list[Path]) -> list[Path]:
    verified_read1s = []
    for read1_file in list_of_R1s:
        read2_file = read1_file.parent.joinpath(read1_file.name.replace("R1", "R2"))
        read3_file = read1_file.parent.joinpath(read1_file.name.replace("R1", "R3"))
        if read2_file.exists() and read3_file.exists():
            verified_read1s.append(read1_file)
        else:
            logger.warning(f"matching R2 and R3 not found for {read1_file}")
    return verified_read1s


# identify all sequencing data that should be parsed for conversion
def parse_directories(folder_list: list[Path], sample_list: list[str]) -> list[Path]:
    all_read1s: list[Path] = []

    # Look into all supplied folders for specific files:
    for folder in folder_list:
        # Look at all of the possible sample names
        for sample in sample_list:
            matching_read1s = [
                f for f in folder.glob("*R1_001.fastq.gz") if sample in f.name
            ]
            all_read1s = list(chain(all_read1s, matching_read1s))
    return verify_sample_from_R1(all_read1s)


# Process through iterator
def batch_iterator(iterator: Iterator[fq.fastq_object], batch_size: int) -> Iterator[list[fq.fastq_object]]:
    entry = True  # Make sure we loop once
    while entry:
        batch: list[fq.fastq_object] = []
        while len(batch) < batch_size:
            try:
                entry = iterator.__next__()
                # logger.info(f"Adding {entry}")
            except StopIteration:
                entry = False
            if not entry:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


# Reformat read for export
def formatRead(title: str, sequence: str, quality: str) -> str:
    return f"@{title}\n{sequence}\n+\n{quality}\n"


def asap_to_kite_v1(trio: list[fq.fastq_object], rc_R2: bool, conjugation: str) -> list[list[str]]: #noqa FBT001
    read1, read2, read3 = trio

    # Parse aspects of existing read
    title1 = read1.head
    sequence1 = read1.body
    quality1 = read1.qstr

    title2 = read2.head
    sequence2 = read2.body
    quality2 = read2.qstr

    title3 = read3.head
    sequence3 = read3.body
    quality3 = read3.qstr

    # process R2
    if rc_R2:
        # Update sequence with reverse complement
        sequence2 = revseq(sequence2)
        # update the quality
        quality2 = quality2[::-1]

    # Recombine attributes based on conjugation logic
    if conjugation == "TotalSeqA":
        new_sequence1 = sequence2 + sequence1[:10]
        new_sequence2 = sequence3

        new_quality1 = quality2 + quality1[:10]
        new_quality2 = quality3

    elif conjugation == "TotalSeqB":
        new_sequence1 = sequence2 + sequence3[:10] + sequence3[25:34]
        new_sequence2 = sequence3[10:25]

        new_quality1 = quality2 + quality3[:10] + quality3[25:34]
        new_quality2 = quality3[10:25]

    # Prepare reads for exporting
    out_fq1 = formatRead(title1, new_sequence1, new_quality1)
    out_fq2 = formatRead(title2, new_sequence2, new_quality2)

    return [[out_fq1, out_fq2]]


@app.callback(invoke_without_command=True)
@app.command(name="asap_o_matic")
def main(
    folder_of_fastqs: Annotated[
        list[Path],
        typer.Option(
            "-f",
            "--fastqs",
            help="Path of folder created by mkfastq or bcl2fastq; can be comma separated that will be collapsed into one output",
            file_okay=False,
            resolve_path=True,
            dir_okay=True,
            readable=True,
        ),
    ],
    sample_name: Annotated[
        list[str],
        typer.Option(
            "--sample",
            "-s",
            help="Prefix of the filenames of FASTQs to select; can be comma separated that will be collapsed into one output",
        ),
    ],
    out: Annotated[
        str,
        typer.Option(
            "--id",
            "-o",
            help="A unique run id, used to name output.",
        ),
    ],
    n_cpu: Annotated[
        int,
        typer.Option(
            "--cores", "-c", help="Number of cores to use for parallel processing."
        ),
    ] = DEFAULT_NUMBER_OF_THREADS,
    n_reads: Annotated[
        int,
        typer.Option(
            "-n",
            "--nreads",
            help="Maximum number of reads to process in one iteration. Decrease this if in a low memory environment.",
        ),
    ] = DEFAULT_MAX_READS_PER_ITERATION,
    rc_R2: Annotated[ #noqa FBT002
        bool,
        typer.Option(
            "--rc-R2/--no-rc-R2",
            "-r/-R",
            help="Should the reverse complement of R2 be used? Pass '--rc-R2' if the reads were generated on a NextSeq or v1.0 chemistry NovaSeq.",
        ),
    ] = False,
    conjugation: Annotated[
        Conjugation,
        typer.Option(
            "-j",
            "--conjugation",
            help="String specifying antibody conjugation; either TotalSeqA or TotalSeqB",
        ),
    ] = Conjugation.TotalSeqA,
    debug: Annotated[ #noqa FBT002
        bool, typer.Option("--debug", help="Print extra information for debugging.")
    ] = False,
    version: Annotated[  # noqa ARG001
        bool,
        typer.Option(
            "--version",
            callback=version_callback,
            help="Print version number.",
        ),
    ] = False,
) -> None:
    """
    IT SLICES
    IT DICES
    IT REFORMATES RAW SEQUENCING DATA FROM CELLRANGER-ATAC INTO SOMETHING USABLE BY OTHER TOOLS
    """

    if not isinstance(conjugation, str):
        conjugation = conjugation.value

    if debug:
        init_logger(verbose=verbosity_level)
    else:
        init_logger(verbose=1)

    read1s_for_analysis = parse_directories(folder_of_fastqs, sample_name)

    # Main loop -- process input reads and write out the processed fastq files
    logger.info("Processing these fastq samples: ")
    for r in read1s_for_analysis:
        logger.info(r.stem)

    outfq1file = f"{out}_R1.fastq.gz"
    outfq2file = f"{out}_R2.fastq.gz"
    with gzip.open(outfq1file, "wt") as out_f1:
        with gzip.open(outfq2file, "wt") as out_f2:
            for read1_file in read1s_for_analysis:
                read2_file = read1_file.parent.joinpath(read1_file.name.replace("R1", "R2"))
                read3_file = read1_file.parent.joinpath(read1_file.name.replace("R1", "R3"))

                # Read in fastq in chunks the size of the maximum user tolerated number
                logger.warning(f"Creating batches for reads in {read1_file}")
                it1 = batch_iterator(fq.read(read1_file), n_reads)
                logger.warning(f"Creating batches for reads in {read2_file}")
                it2 = batch_iterator(fq.read(read2_file), n_reads)
                logger.warning(f"Creating batches for reads in {read3_file}")
                it3 = batch_iterator(fq.read(read3_file), n_reads)

                for batch_read1 in it1:
                    batch_read2 = it2.__next__()
                    batch_read3 = it3.__next__()

                    pool = Pool(processes=n_cpu)
                    convert_partial = partial(
                        asap_to_kite_v1, rc_R2=rc_R2, conjugation=conjugation
                    )
                    pm = pool.map(convert_partial, zip(batch_read1, batch_read2, batch_read3, strict=True))
                    pool.close()

                    # process and write out
                    fq_data = list(map("".join, zip(*[item.pop(0) for item in pm], strict=True)))
                    out_f1.writelines(fq_data[0])
                    out_f2.writelines(fq_data[1])

    logger.info("Done!")
