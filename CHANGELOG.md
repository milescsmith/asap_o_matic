# CHANGELOG

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


[1.1.0]: https://github.com/milescsmith/asap_o_matic/compare/1.0.0..1.1.0
[1.0.0]: https://github.com/milescsmith/asap_o_matic/releases/tag/1.0.0