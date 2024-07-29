import importlib.resources as ir

import pytest
from asap_o_matic import verify_sample_from_R1


@pytest.fixture
def read1_from_cellranger():
    return ir.files("tests").joinpath("data", "cellranger", "test_S_S3_R1_001.fastq.gz")


@pytest.fixture
def read1_from_bcl_convert():
    return ir.files("tests").joinpath("data", "bcl-convert", "test_S_S3_R1_001.fastq.gz")


@pytest.fixture
def read1_incomplete():
    return ir.files("tests").joinpath("data", "incomplete", "test_S_S3_R1_001.fastq.gz")


def test_verify_sample_from_R1_bcl_convert(
    read1_from_bcl_convert,
):
    assert verify_sample_from_R1([read1_from_bcl_convert], fastq_source="bcl-convert") == [read1_from_bcl_convert]


def test_verify_sample_from_R1_cellranger(read1_from_cellranger):
    assert verify_sample_from_R1([read1_from_cellranger], fastq_source="cellranger") == [read1_from_cellranger]


def test_verify_sample_incomplete(read1_incomplete):
    assert verify_sample_from_R1([read1_incomplete]) == []


def test_verify_bad_source(read1_incomplete):
    with pytest.raises(ValueError):
        assert verify_sample_from_R1([read1_incomplete], fastq_source="narnia")
