import importlib.resources as ir

import pytest
from asap_o_matic import parse_directories


@pytest.fixture
def bcl_convert_read_path():
    return ir.files("tests").joinpath("data", "bcl-convert")


@pytest.fixture
def cellranger_read_path():
    return ir.files("tests").joinpath("data", "cellranger")


def test_parse_directories_bcl_convert(bcl_convert_read_path):
    file_list = parse_directories(
        folder_list=[bcl_convert_read_path], sample_list=["test_S_S3_"], fastq_source="bcl-convert"
    )
    assert len(file_list) == 1
    assert "bcl-convert/test_S_S3_R1_001.fastq.gz" in str(file_list[0])


def test_parse_directories_cellranger(cellranger_read_path):
    file_list = parse_directories(
        folder_list=[cellranger_read_path], sample_list=["test_S_S3_"], fastq_source="cellranger"
    )
    assert len(file_list) == 1
    assert "cellranger/test_S_S3_R1_001.fastq.gz" in str(file_list[0])


def test_parse_directories_wrong_sample(bcl_convert_read_path):
    file_list = parse_directories(
        folder_list=[bcl_convert_read_path.joinpath("bcl-fastq")],
        sample_list=["experiment_S_S3_"],
        fastq_source="bcl-convert",
    )
    assert len(file_list) == 0
