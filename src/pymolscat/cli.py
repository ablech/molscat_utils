import click
import os
from pathlib import Path
from .wrapper import MolscatWrapper
#from .io     import write_cross_sections, write_s_matrices

@click.command()
@click.argument('system', type=str)
@click.argument('runfolder',
                type=click.Path(
                    exists=True, file_okay=False, dir_okay=True,
                    writable=True, readable=True,
                    ),
                )
@click.option('--molscat-bin', envvar='MOLSCAT_BIN', type=click.Path(exists=True), default=None,
              help="Path to the Molscat executables. Defaults to the current directory or the MOLSCAT_BIN environment variable."
              )
@click.option('--quiet', type=bool, default=False,
              help="Don't print Molscat output to stdout."
              )
def main(system, runfolder, molscat_bin, quiet):
    """Run MOLSCAT.

    \b
    Arguments:
        SYSTEM :
            Name of the collision system. This is appended to the molscat executable name as f'molscat-{system}'
        RUNFOLDER :
            Path to the runfolder where the input files are located and output will be written.
    """
    runfolder = Path(runfolder)
    if molscat_bin is None:
        molscat_bin = Path.cwd() #if "MOLSCAT_BIN" not in os.environ else Path(os.environ["MOLSCAT_BIN"])
    else:
        molscat_bin = Path(molscat_bin)
    executable = molscat_bin / f"molscat-{system}"
    assert executable.exists(), f"Molscat executable {executable} does not exist"

    wrapper = MolscatWrapper(executable, quiet=quiet)
    wrapper.run(runfolder)

    return 0