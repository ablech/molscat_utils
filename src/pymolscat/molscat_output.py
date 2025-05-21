from dataclasses import dataclass, field
from typing import List, Tuple, Dict
import numpy as np


@dataclass
class ScatteringBlock:
    """Results of the calculation for one (J, parity) block."""
    Jtot: int
    sym: int
    energies: List
    n_channels: List[int] # List of number of channels
    n_open: List[int] # List of number of open channels
    #channels: List[List[Tuple[int, ...]]]     # e.g. [[(i,v1,v2,...), ...]]
    #open_channels: List[List[Dict]]
    #open_idx: List[int] # indices of open channels from the channel list
    S: List[np.ndarray]                # shape (Nch, Nch), dtype=complex
    pcs: List[np.ndarray] # list of partial cross sections


@dataclass
class MolscatResult:
    """Aggregate result of a full MOLSCAT run."""
    para: object
    type: str # interaction type as extracted from output
    mu: float # reduced mass
    energies: List # List of input energies
    pair_energies: List # List of threshold energies
    #quantum_numbers: List[str] # List of pair state quantum numbers
    Jtot: List[int] = field(init=False) # List of Jtot values
    symmetries: List[int] # List of symmetries
    blocks: List[ScatteringBlock]
    #exit_code: int
    runtime: float
    memory: float
    n_blocks: int = field(init=False)
    ics: List[np.ndarray] # List of total integrated cross sections
    acs: List[List[np.ndarray]] # For each energy, the list of accumulated integrated cross sections

    def __post_init__(self):
        self.n_blocks = len(self.blocks)
        self.Jtot = np.arange(self.para.JTOTL, self.para.JTOTU+1, self.para.JSTEP)


@dataclass
class MolscatInputParameters:
    """Molscat input variables extracted from the output"""
    JTOTL: int
    JTOTU: int
    JSTEP: int
    IPRINT: int
