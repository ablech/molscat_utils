MOLSCAT utilities
=================

This package is a collection of utility programs for running and processing
scattering calculations with `MOLSCAT`_.

.. _MOLSCAT: https://github.com/molscat/molscat

Key features:

- MOLSCAT wrapper with automatic input file generation
- Parsing  MOLSCAT output into MolscatResult objects
- Command line interface for wrapper and parser


Installation
------------

After cloning the repository,

    git clone https://github.com/ablech/molscat_utils.git

it can be installed with ``pip`` into your active python environment by running

    pip install -e ./molscat_utils

Publication to PyPi is planned.


Command line usage
------------------

Wrapper
^^^^^^^

Clone and compile `MOLSCAT`_.

The wrapper needs to know, where to find the molscat executables. This can be
controlled via the ``MOLSCAT_BIN`` environment variable.
For example, if molscat compiles binaries into its `source_code` directory
(the default), define

  export MOLSCAT_BIN=path_to_your_molscat_clone/source_code

The call signature of the wrapper is

  molscat [OPTIONS] SYSTEM RUNFOLDER

The wrapper calls the molscat binary ``molscat-SYSTEM`` found in ``$MOLSCAT_BIN``.
``RUNFOLDER`` needs to be a directory, from which input files are read and to which
all output is written. The runfolder should contain an input template file `molscat.in`
and optionally an additional configuration file ``config.json`` in JSON format.
The input template file needs to follow the namelist syntax of a molscat input file.
It may, however, contain variables, e.g., ``$[JTOTU]``, which are replaced by the wrapper
when writing the actual input file (``inp.dat``).
The values of these variables are inferred from the content of ``config.json``.


Output Parser
^^^^^^^^^^^^^

The call signature of the wrapper is

  molscat-parser [OPTIONS] FILE

The resulting ``MolscatResult`` object is then stored as ``FILE.parsed.pkl`` in the same
directory as the original output file.

For example, to parse a molscat output file, `./runfolder/out.dat`, run

  molscat-parser ./runfolder/out.dat

The resulting ``MolscatResult`` object is stored as is stored as
``./runfolder/out.dat.parsed.pkl``.

In a python script, the `MolscatResult` object can be load with the ``load()``
function of the ``pickle`` package.
