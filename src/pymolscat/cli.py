import click
import os
import json
import pickle
from pathlib import Path
from .wrapper import MolscatWrapper
from .parser import parse_output
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
@click.option('--quiet', is_flag=True,
              help="Don't print Molscat output to stdout."
              )
def wrapper(system, runfolder, molscat_bin, quiet):
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


@click.command()
@click.argument('file',
                type=click.Path(
                    exists=True, file_okay=True, dir_okay=False,
                    writable=False, readable=True,
                    ),
                )
# NOTE: json serialization currently not possible due to numpy arrays
@click.option('--write-json', is_flag=True,
              help="Save parsed output in JSON format instead of as pickled object."
              )
def parser(file, write_json):
    """Parse MOLSCAT output file.

    \b
    Arguments:
        FILE :
            The MOLSCAT output file to parse.

    The parsed output is stored as FILE.parsed.pkl or FILE.parsed.json in the
    same directory as FILE.
    """
    file = Path(file).absolute()
    with open(file, 'r') as f:
        output = f.readlines()
    MolscatResult, output_data = parse_output(output)
    if write_json:
        with open(f"{file}.parsed.json", 'w') as f:
            json.dump(output_data, f, indent=4)
    else:
        with open(f"{file}.parsed.pkl", 'wb') as f:
            pickle.dump(output_data, f)
    return 0