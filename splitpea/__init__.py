# splitpea/__init__.py

from .main import calculate_delta_psi, combine_spliced_exon, rewire, plot, stats
from .get_background_ppi import get_background
from .get_consensus_network import get_consensus_network

__all__ = [
    "calculate_delta_psi",
    "combine_spliced_exon",
    "get_background",
    "get_consensus_network", 
    "rewire",
    "plot",
    "stats"
]
