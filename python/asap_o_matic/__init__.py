"""
asap-o-matic

A *lightly* altered version of the asap_to_kite_v2.py script
written by Caleb Lareau and
found at https://github.com/caleblareau/asap_to_kite
"""

from asap_o_matic.__main__ import __version__, asap_to_kite, parse_directories, verify_sample_from_R1

# from rich.console import Console

__all__ = ["verify_sample_from_R1", "parse_directories", "asap_to_kite", "__version__"]
