import gzip
import sys
import tempfile
import pysam
from enum import Enum
from itertools import chain
from multiprocessing import cpu_count
from pathlib import Path
from typing import Annotated, Literal, Optional

import fastq as fq
import typer
from joblib import Parallel, delayed
from loguru import logger
from revseq import revseq

# from rich.traceback import install
from tqdm.auto import tqdm

from asap_o_matic import app, verbosity_level, version_callback
from asap_o_matic.asap_o_matic import format_read, rearrange_reads
from asap_o_matic.logger import init_logger

# install()


class Conjugation(str, Enum):
    TotalSeqA = "TotalSeqA"
    TotlaSeqB = "TotalSeqB"


class FastqSource(str, Enum):
    cellranger = "cellranger"
    bclconvert = "bcl-convert"


DEFAULT_NUMBER_OF_THREADS = cpu_count()
DEFAULT_MAX_READS_PER_ITERATION = 1000000


def verify_sample_from_R1(
    list_of_R1s: list[Path], fastq_source: Literal["cellranger", "bcl-convert"] = "bcl-convert"
) -> list[Path]:
    """Verify R1/R2/R3 are present for nominated samples

    Parameters
    ----------
    list_of_R1s : list[Path]
        Location(s) of the FASTQs corresponding to read1

    Returns
    -------
    list[Path]
        Location(s) of read1 FASTQs that have matching read2 and read3
    """
    verified_read1s = []
    for read1_file in list_of_R1s:
        logger.debug(f"looking for matches for {read1_file}...")
        match fastq_source:
            case "bcl-convert":
                read2_file = read1_file.parent.joinpath(read1_file.name.replace("R1", "R2"))
                index1_file = read1_file.parent.joinpath(read1_file.name.replace("R1", "I1"))
                index2_file = read1_file.parent.joinpath(read1_file.name.replace("R1", "I2"))
                if read2_file.exists() and index1_file.exists() and index2_file.exists():
                    logger.debug(
                        f"found read2 at {read2_file!s}, index1 at {index1_file!s}, and index2 at {index2_file!s}"
                    )
                    verified_read1s.append(read1_file)
                else:
                    logger.warning(f"matching R2, I1, and/or I2 not found for {read1_file}")
            case "cellranger":
                read2_file = read1_file.parent.joinpath(read1_file.name.replace("R1", "R2"))
                read3_file = read1_file.parent.joinpath(read1_file.name.replace("R1", "R3"))
                if read2_file.exists() and read3_file.exists():
                    logger.debug(f"found read2 at {read2_file!s}, read3 at {read3_file!s}")
                    verified_read1s.append(read1_file)
                else:
                    logger.warning(f"matching R2 and/or R3 not found for {read1_file!s}")
            case _:
                msg = f"{fastq_source} is not a recognized source for bcl->fastq conversion"
                raise ValueError(msg)
    return verified_read1s


def parse_directories(
    folder_list: list[Path], sample_list: list[str], fastq_source: Literal["cellranger", "bcl-convert"] = "bcl-convert"
) -> list[Path]:
    """Identify all sequencing data that should be parsed for conversion

    Parameters
    ----------
    folder_list : list[Path]
        _description_
    sample_list : list[str]
        _description_

    Returns
    -------
    list[Path]
        _description_
    """
    all_read1s: list[Path] = []

    # Look into all supplied folders for specific files:
    for folder in folder_list:
        # Look at all of the possible sample names
        for sample in sample_list:
            matching_read1s = [f for f in folder.glob("*R1_001.fastq.gz") if sample in f.name]
            all_read1s = list(chain(all_read1s, matching_read1s))
    logger.debug(f"Found {', '.join([str(_) for _ in all_read1s])}")
    return verify_sample_from_R1(all_read1s, fastq_source)


# def formatRead(title: str, sequence: str, quality: str) -> str:
#     # Reformat read for export
#     return f"@{title}\n{sequence}\n+\n{quality}\n"


def asap_to_kite(
    fastq_list: list[fq.fastq_object],
    rc_R2: bool,
    conjugation: str,
    new_read1_handle: gzip.GzipFile,
    new_read2_handle: gzip.GzipFile,
) -> None:  # list[list[str]]:  # nFBT001
    """Rearrange the disparate portions of CITE-seq reads that are split among the R1, R2, and R3 of ASAP-seq data
    into something that Kallisto/Bustools or CITE-seq-Count can process

    Parameters
    ----------
    trio : list[fq.fastq_object]
        _description_
    rc_R2 : bool
        _description_
    conjugation : str
        The type of CITE-seq antibodies used, either TotalSeqA (if using the 10x Genomics scATAC-seq kit) or TotalSeqB
        (if using the 10x Genomics Multiome kit)

    Returns
    -------
    list[list[str]]
        A list of the reformatted read1 and read2 pairs
    """
    read1, read2, read3 = fastq_list

    # Parse aspects of existing read
    title1 = read1.head
    sequence1 = read1.body
    quality1 = read1.qstr

    title2 = read2.head
    sequence2 = read2.body
    quality2 = read2.qstr

    # title3 = read3.head
    sequence3 = read3.body
    quality3 = read3.qstr

    # process R2
    if rc_R2:
        # Update sequence with reverse complement
        sequence2 = revseq(sequence2)
        # update the quality
        quality2 = quality2[::-1]

    # Recombine attributes based on conjugation logic

    new_sequence1, new_sequence2, new_quality1, new_quality2 = rearrange_reads(
        sequence1, sequence2, sequence3, quality1, quality2, quality3, conjugation
    )

    with open(new_read1_handle, "ab") as f1, open(new_read2_handle, "ab") as f2:
        f1.write(format_read(title1, new_sequence1, new_quality1).encode())
        f2.write(format_read(title2, new_sequence2, new_quality2).encode())
    # out_fq1 = formatRead(title1, new_sequence1, new_quality1)
    # out_fq2 = formatRead(title2, new_sequence2, new_quality2)

    # return [[out_fq1, out_fq2]]


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
            exists=True,
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
    fastq_source: Annotated[
        FastqSource,
        typer.Option(
            "--fastq_source",
            "-a",
            help="Name of the program used to convert bcls to FASTQs. Cellranger mkfastq creates R1, R2, R3, and I3 files while bcl-convert creates R1, I1, R2, I2 files.",
        ),
    ] = FastqSource.cellranger,
    outdir: Annotated[
        Optional[Path],
        typer.Option(
            "--outdir",
            "-d",
            help="Directory to save files to.  If none is give, save in the directory from which the script was called.",
            file_okay=False,
            resolve_path=True,
            dir_okay=True,
            readable=True,
        ),
    ] = None,
    n_cpu: Annotated[
        int,
        typer.Option("--cores", "-c", help="Number of cores to use for parallel processing."),
    ] = DEFAULT_NUMBER_OF_THREADS,
    rc_R2: Annotated[  # nFBT002
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
    debug: Annotated[bool, typer.Option("--debug", help="Print extra information for debugging.")] = False,  # nFBT002
    version: Annotated[  # ARG001
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
    IT REFORMATES RAW SEQUENCING DATA FROM CELLRANGER-ATAC INTO SOMETHING USABLE BY OTHER TOOLS\n
    USE THE WHOLE READ WITHOUT ANY FISH WASTE
    WITHOUT ANY SCALING, CUTTING, OR GUTTING
    """

    if not isinstance(conjugation, str):
        conjugation = conjugation.value

    logger.remove()
    if debug:
        logger.add(
            sys.stderr,
            format="* <red>{elapsed}</red> - <cyan>{module}:{file}:{function}</cyan>:<green>{line}</green> - <yellow>{message}</yellow>",
            colorize=True,
        )
        init_logger(verbose=verbosity_level)
    else:
        logger.add(sys.stderr, format="* <yellow>{message}</yellow>", colorize=True)
        init_logger(verbose=1, msg_format="<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>")

    read1s_for_analysis = parse_directories(folder_of_fastqs, sample_name, fastq_source)

    # Main loop -- process input reads and write out the processed fastq files
    logger.info("Processing these fastq samples: ")
    for r in read1s_for_analysis:
        logger.info(r.name)

    if outdir is None:
        outdir = Path().cwd()
    outfq1file = outdir.joinpath(f"{out}_R1.fastq.gz")
    outfq2file = outdir.joinpath(f"{out}_R2.fastq.gz")
    tempfq1file = tempfile.NamedTemporaryFile()
    tempfq2file = tempfile.NamedTemporaryFile()
    
    for read1_file in read1s_for_analysis:
        logger.info(f"Processing reads associated with {read1_file.name}")
        match fastq_source:
            case "bcl-convert":
                read2_file = read1_file.parent.joinpath(read1_file.name.replace("R1", "I2"))
                read3_file = read1_file.parent.joinpath(read1_file.name.replace("R1", "R2"))
            case "cellranger":
                read2_file = read1_file.parent.joinpath(read1_file.name.replace("R1", "R2"))
                read3_file = read1_file.parent.joinpath(read1_file.name.replace("R1", "R3"))

        # No need to chunk when we use an iterator
        read1 = fq.read(read1_file)
        read2 = fq.read(read2_file)
        read3 = fq.read(read3_file)

        if n_cpu > 1:
            parallel = Parallel(n_jobs=n_cpu, return_as="list")
            _ = parallel(
                delayed(asap_to_kite)(
                    (a, b, c),
                    rc_R2=rc_R2,
                    conjugation=conjugation,
                    new_read1_handle=tempfq1file.name,
                    new_read2_handle=tempfq2file.name,
                )
                for a, b, c in tqdm(zip(read1, read2, read3, strict=True), unit="Reads")
            )
        else:
            _ = [
                asap_to_kite(
                    (a, b, c),
                    rc_R2=rc_R2,
                    conjugation=conjugation,
                    new_read1_handle=tempfq1file.name,
                    new_read2_handle=tempfq2file.name,
                )
                for a, b, c in tqdm(zip(read1, read2, read3, strict=True), unit="Reads")
            ]

        logger.info("Finished rearranging reads. Now compressing...")
        pysam.tabix_compress(tempfq1file.file.name, outfq1file, force=True)
        pysam.tabix_compress(tempfq2file.file.name, outfq2file, force=True)
        logger.info("Finished compressing.")

    logger.info("Done!")
