from importlib import metadata
import json
from pathlib import Path
import subprocess


class MolscatWrapper(object):
    """
    Wrapper for the Molscat program.
    """

    def __init__(self, molscat_bin: str):
        """
        Initialize the MolscatWrapper class.
        """
        self.molscat = Path(molscat_bin)
        self.config = {}

    def print_header(self, log):
        """
        Print information about MolscatWrapper and the MOLSCAT binary at the
        beginning of the log file and to stdout.
        """
        version = metadata.version("pymolscat")
        header = (
            f"#   Running MolscatWrapper, version {version}\n"
            f"#   Calling MOLSCAT binary: {self.molscat.name}\n"
            f"#   from: {self.molscat.parent}\n\n"
        )
        log.write(header)
        print(header)

    def run(self, runfolder: Path, log_file: str = 'out.dat'):
        """
        Run the molscat program.

        Args:
            command (str): The command to execute.
            log_file (Path): Path to the log file where stdout will be written.
        """
        # Prepare the input template
        self.parse_config(runfolder)
        input_data = self.write_molscat_input(runfolder)

        command = [self.molscat]
        with open(runfolder / log_file, 'w') as log:
            self.print_header(log)
            process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            stdin=subprocess.PIPE,
            text=True,
            )
            # Send the input data to the process
            process.stdin.write(input_data)
            process.stdin.close()
            # Stream stdout to the log file in real time
            output = []
            for line in process.stdout:
                log.write(line)
                log.flush()  # Ensure the log file is updated immediately
                print(line, end='')  # Optionally print to the console
                output.append(line)  # Collect the output

            process.wait()  # Wait for the process to complete

        return ''.join(output)  # Return collected output as a single string

    def parse_config(self, runfolder: Path, filename: str='config.json'):
        """
        Read and parse the configuration file.
        """
        try:
            with open(runfolder / filename, 'r') as f:
                config_data = json.load(f)
        except:
            self.config = {}
            return

        MAX_LINE_LENGTH = 79
        indentation = ' ' * 13
        jlevel = ''
        line = ''
        for v, j_data in config_data['pair levels']['v'].items():
            for j in range(j_data['jmin'], j_data['jmax'] + 1, j_data['jstep']):
                tuple_str = f"{j},{v}, "
                if len(line) + len(tuple_str) > MAX_LINE_LENGTH:
                    jlevel += line.rstrip() + '\n' + indentation
                    line = tuple_str
                else:
                    line += tuple_str
        jlevel += line.rstrip()  # Add the remaining part of the line
        jlevel = jlevel.rstrip(', ')  # Remove double comma after last entry

        nlevel = len(jlevel.replace('\n', '').split(', '))

        # figure out index of reference level
        # TODO

        self.config['nlevel'] = nlevel
        self.config['jlevel'] = jlevel
        #self.config['iref'] = iref

    def write_molscat_input(self, runfolder: Path, filename: str='molscat.in'):
        """
        Write the Molscat input file.
        """
        with open(runfolder / filename, 'r') as f:
            # Read the template file
            template = f.read()

        if len(self.config) > 0:
            # Define the replacements dictionary
            replacements = {
                "$[NLEVEL]": str(self.config['nlevel']),
                "$[JLEVEL]": str(self.config['jlevel']),
                #"$[IREF]": str(self.config['iref']),
            }

            # Perform the replacements dynamically
            for placeholder, value in replacements.items():
                template = template.replace(placeholder, value)

        with open(runfolder / 'inp.dat', 'w') as f:
            # Write the input data to the file
            f.write(template)

        return template