import datetime
from pathlib import Path
from sys import modules, stdout
from typing import TextIO

from loguru import logger

parent_module = modules[".".join(__name__.split(".")[:-1]) or "__main__"]


def init_logger(verbose: int = 0, msg_format: str | None = None, save: bool = True) -> int:
    if msg_format is None:
        msg_format = "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>"

    timezone = datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo

    match verbose:
        case 3:
            log_level = "DEBUG"
            backtrace = True
            diagnose = True
        case 2:
            log_level = "INFO"
            backtrace = False
            diagnose = False
        case 1:
            log_level = "WARNING"
            backtrace = False
            diagnose = False
        case _:
            log_level = "ERROR"
            backtrace = False
            diagnose = False

    output_sink: Path | TextIO = (
        Path(f"asap_o_matic_{datetime.datetime.now(tz=timezone).strftime('%d-%m-%Y--%H-%M-%S')}.log")
        if save
        else stdout
    )

    return logger.add(
        sink=output_sink,
        format=msg_format,
        level=log_level,
        backtrace=backtrace,
        diagnose=diagnose,
        filter="asap-o-matic",
    )
