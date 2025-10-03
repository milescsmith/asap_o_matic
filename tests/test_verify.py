import importlib.resources as ir
from pathlib import Path

import pytest
from asap_o_matic import verify_sample_from_R1
from asap_o_matic.__main__ import FastqSource


@pytest.fixture
def read1_from_cellranger():
    return ir.files("tests").joinpath("data", "cellranger", "test_S_S3_R1_001.fastq.gz") # type: ignore


@pytest.fixture
def read1_from_bcl_convert():
    return ir.files("tests").joinpath("data", "bcl-convert", "test_S_S3_R1_001.fastq.gz") # type: ignore


@pytest.fixture
def read1_incomplete():
    return ir.files("tests").joinpath("data", "incomplete", "test_S_S3_R1_001.fastq.gz") # type: ignore


def test_verify_sample_from_R1_bcl_convert(
    read1_from_bcl_convert: Path,
):
    assert verify_sample_from_R1([read1_from_bcl_convert], fastq_source=FastqSource.bclconvert) == [read1_from_bcl_convert]


def test_verify_sample_from_R1_cellranger(read1_from_cellranger: Path):
    assert verify_sample_from_R1([read1_from_cellranger], fastq_source=FastqSource.cellranger) == [read1_from_cellranger]


def test_verify_sample_incomplete(read1_incomplete: Path):
    assert verify_sample_from_R1([read1_incomplete]) == []


def test_verify_bad_source(read1_incomplete: Path):
    with pytest.raises(ValueError):
        assert verify_sample_from_R1([read1_incomplete], fastq_source="narnia")
