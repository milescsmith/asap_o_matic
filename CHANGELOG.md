# CHANGELOG

# [2.5.0] - 2025-10-03

## Added

- An actual README.md, with updates to reflect how to actually use asap-o-matic
- Added several instances of casting from a `pathlib.Path` object to a `str`
- Added type hints to the unit tests

## Changed

- Increased minimum Python version to 3.11
- Increased Python versions tested to include 3.11 to 3.14
- Replaced use of pdm with uv
- Update noxfile to use uv
- Replaced use of mypy with ty
- Renamed Rust submodule to `_rust_o_matic`
- Updated Rust dependencies

# [2.4.0] - 2025-02-05

## Added

- Added parameter to asap_o_matic to save the log file

## Fixed

- And we're back to actually writing the output to disk.
- Logging actually works now

## Changed

- Update dependencies
- Added a note about where output is being written to

# [2.3.0] - 2025-02-03

## Changed

- Removed tqdm.  Already importing Rich, might as well use `rich.progress`

## Fixed

- If output directory does not exist, attempt to create it.

# [2.2.0] - 2024-07-29

## Added

- Unit tests

## Changed

- Rust code cleaned a little


# [2.1.1] - 2024-04-18

## Added

- Added messages to indicate rearrangement has completed and compressing the files had begun

# [2.1.0] - 2024-04-18

## Fixed

- I guess using write/append to a gzipped file doesn't work? Replaced the way output FASTQs were being written
so that they first go to a temporary file and are then bgzipped

## Changed

- Using PDM for the Python package management side of things

# [2.0.0] - 2024-04-18

## Changed

- Switch build-backend from hatch to maturin
- Rearrange module structure to account for now being a joint Rust/Python package
- Rewrote `formatRead` and part of `asap_to_kite` in Rust

# [1.4.1] - 2024-04-10

## Fix

- Changed the `bcl_source` argument to `asap-o-matic` to `fastq_source`

# [1.4.0] - 2024-04-10

## Added

- Add ability to handle FASTQs created by bcl2fastq/bcl-convert (really, the only difference in the naming but it was
still annoying)

## Changed

- Write gzipped FASTQ files
- Write output FASTQ files during the reformatting loop instead of all at once at the end

# [1.3.0] - 2023-11-28

## Changed

- Turns out `process_map` is slower than just a simple list comprehension so switched to Joblib
- Removed chunking because?

# [1.2.0] - 2023-11-27

## Changed

- Made logging a litle less overwhelming and actually reflect passing the "debug" argument
- Replaced `multiprocessing.Pool` with `tqdm.contrib.concurrent.process_map`
- Reduced the DEFAULT_MAX_READS_PER_ITERATION from 100 to 1 million

## Removed

- Removed now unused `batch_iterator`

# [1.1.0] - 2023-11-27

## Added

- Ability to set the output directory
- Some docstrings

## Changed

- Replacing custom `batch_iterator` with `more_itertools.chunked`
- Overhauled logging submodule
- Changed how the next item in an iterator is called

# [1.0.0] - 2023-10-19

## Changed

- Changed name
- Reworked previous script to act as an installble package
- Added logging, (at least some) typing, etc...

[2.5.0]: https://github.com/milescsmith/asap_o_matic/compare/2.4.0..2.5.0
[2.4.0]: https://github.com/milescsmith/asap_o_matic/compare/2.3.0..2.4.0
[2.3.0]: https://github.com/milescsmith/asap_o_matic/compare/2.2.0..2.3.0
[2.2.0]: https://github.com/milescsmith/asap_o_matic/compare/2.1.1..2.2.0
[2.1.1]: https://github.com/milescsmith/asap_o_matic/compare/2.1.0..2.1.1
[2.1.0]: https://github.com/milescsmith/asap_o_matic/compare/2.0.0..2.1.0
[2.0.0]: https://github.com/milescsmith/asap_o_matic/compare/1.4.1..2.0.0
[1.4.1]: https://github.com/milescsmith/asap_o_matic/compare/1.4.0..1.4.1
[1.4.0]: https://github.com/milescsmith/asap_o_matic/compare/1.3.0..1.4.0
[1.2.0]: https://github.com/milescsmith/asap_o_matic/compare/1.2.0..1.3.0
[1.2.0]: https://github.com/milescsmith/asap_o_matic/compare/1.1.0..1.2.0
[1.1.0]: https://github.com/milescsmith/asap_o_matic/compare/1.0.0..1.1.0
[1.0.0]: https://github.com/milescsmith/asap_o_matic/releases/tag/1.0.0
