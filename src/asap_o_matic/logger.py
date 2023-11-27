import datetime
from pathlib import Path
from sys import modules, stdout
from typing import TextIO

from loguru import logger

parent_module = modules[".".join(__name__.split(".")[:-1]) or "__main__"]

def init_logger(verbose: int = 0, msg_format: str | None = None, save: bool = True, module_name: str=parent_module) -> None:
    if msg_format is None:
        msg_format = "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan>·-·<level>{message}</level>"

    timezone = datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo

    match verbose:
        case 3:
            log_level = "DEBUG"
        case 2:
            log_level = "INFO"
        case 1:
            log_level = "WARNING"
        case _:
            log_level = "ERROR"

    output_sink: Path | TextIO = (
        Path(f"module_name_{datetime.datetime.now(tz=timezone).strftime('%d-%m-%Y--%H-%M-%S')}.log") if save else stdout
    )

    logger.add(sink=output_sink, format=msg_format, level=log_level, backtrace=True, diagnose=True)
