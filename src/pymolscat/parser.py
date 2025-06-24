import re
from typing import TextIO, List, Dict
import numpy as np
import xarray as xr
from .molscat_output import MolscatInputParameters, ScatteringBlock, MolscatResult


# TODO: add sanity checks:
# - number of open channels is <= number of pair states/levels
# - dimension of S matrix corresponds to number of open channels
# - dimension of printed partial cross section matrix corresponds to number of pair levels in open channel list
# - dimension of printed accumulated cross section matrix corresponds to number of pair states/levels
# - dimension of printed total integrated cross section matrix corresponds to maximum number of pair levels of all partial cross sections,
#   i.e., zero-entries are omitted

INITIALIZATION_PATTERS = {
    'separator': re.compile(r'^s*'+59*r'==', flags=re.MULTILINE), # see mol.driver.f90, format 1060
}

CALCULATION_PATTERNS = {
    'label_line': re.compile(r'^\s*=+\s(.*?)\s=', flags=re.MULTILINE), # see mol.driver.f90, format F710
    'block_start': re.compile(r'^\s*\*+\s+ANGULAR MOMENTUM JTOT  =\s*(\d+)\s+AND SYMMETRY BLOCK  =\s*(\d+)', flags=re.MULTILINE),
    'channels': re.compile(r'CHANNEL FUNCTION LIST'),
    'energy_block_separator': re.compile(r'^\s*'+58*r'- '+'-', flags=re.MULTILINE), # see mol.driver.f90, format 2801
}

ENERGY_BLOCK_SECTION_PATTERNS = {
    'propagation': re.compile(r'^ *COUPLED EQUATIONS PROPAGATED AT *ENERGY', flags=re.MULTILINE),
    'open_channels': re.compile(r'^ *OPEN CHANNEL * WVEC', flags=re.MULTILINE),
    's_matrix': re.compile(r'^ *ROW +COL + S\*\*2'),
    'pcs': re.compile(r'^ *\* \* \* \* \* \* \* \* \* \*  STATE-TO-STATE PARTIAL CROSS SECTIONS', flags=re.MULTILINE),
    'ics': re.compile(r'^ *\* \* \* \* \* \* \* \* \* \* STATE-TO-STATE INTEGRAL CROSS SECTIONS', flags=re.MULTILINE),
}

SUMMARY_SECTION_PATTERNS = {
    'ics': re.compile(r'^\s*STATE-TO-STATE INTEGRAL CROSS SECTIONS IN ANGSTROM', flags=re.MULTILINE),
    'footer': re.compile(r'^\s*-'+5*r'--- MOLSCAT ---', flags=re.MULTILINE), # see mol.driver.f90, format 1001
}


def parse_output(stream: TextIO | List[str]) -> Dict[str, List[List[str]]]:
    stream = iter(stream)  # Ensure the stream is an iterator
    init_data = parse_initialization(stream)
    calc_data = parse_calculation(stream)
    summary_data = parse_summary(stream)
    output_data = {'init': init_data, 'calc': calc_data, 'summary': summary_data}

    debug = False

    if debug:
        try:
            MolscatResult = _process_result_data(
                output_data
            )
            return MolscatResult, output_data
        except:
            print('Error in output processing.')
            print('Returning None instead')
            return None, output_data
    else:
        MolscatResult = _process_result_data(
                output_data
            )
        return MolscatResult, output_data


def parse_initialization(stream: TextIO | List[str]) -> Dict[str, List[List[str]]]:
    data = {}
    read_energies = False
    read_pair_states = False

    skip = 0
    for line in stream:
        line = line.strip()
        if not line: continue # skip empty lines

        if skip > 0:
            skip -= 1
            continue

        if line.startswith('INITIALIZATION DONE'):
            # End of initialization: extrat time and memory usage and return.
            match = re.search(r'TIME WAS\s*(\d+\.\d+) CPU SECS\.\s*(\d+) WORDS', line)
            if match:
                runtime = float(match.group(1))
                words_used = int(match.group(2))
                data['runtime'] = runtime
                data['words_used'] = words_used
            return data

        # Read miscellaneous initialization data

        # Read header information
        match = re.search(r'\| +Run on +(.+\S) +at +(\d\d:\d\d:\d\d) +\|', line)
        if match:
            data['date'] = match.group(1)
            data['time'] = match.group(2)
            continue
        match = re.search(r'\| +Version +(\d+\.\d*) +\|', line)
        if match:
            data['version'] = match.group(1)
            continue

        match = re.search(r'PRINT LEVEL \(IPRINT\) = *(\d+)', line)
        if match:
            data['IPRINT'] = int(match.group(1))
            continue

        match = re.search(r'REDUCED MASS FOR INTERACTION = *(\d+\.\d+).+\((.+)\)', line)
        if match:
            data['reduced mass'] = float(match.group(1))
            data['mass unit'] = match.group(2)
            continue

        match = re.search(r'INTERACTION TYPE IS +(.+)\.', line)
        if match:
            data['interaction type'] = match.group(1)
            continue

        # Read pair levels
        match = re.search(r'EACH PAIR STATE IS LABELLED BY *(\d+) QUANTUM NUMBER', line)
        if match:
            num_quantum_numbers = int(match.group(1))
            data['num_quantum_numbers'] = num_quantum_numbers
            continue
        if line.startswith('PAIR STATE     PAIR STATE QUANTUM NUMBERS    PAIR LEVEL'):
            read_pair_states = True
            pair_states = []
            skip = 1 # skip header line
            continue
        if read_pair_states and INITIALIZATION_PATTERS['separator'].match(line):
            read_pair_states = False
            data['pair_states'] = pair_states
            continue
        elif read_pair_states:
            split = line.split()
            assert len(split) == 3 + num_quantum_numbers
            pair_states.append(
                [int(col) for col in split[:-1]] + [float(split[-1])]
            )
            continue

        match = re.search(r'POTENTIAL RETURNED IN UNITS OF EPSIL *= * (\d+\.\d+) +CM-1', line)
        if match:
            data['pot unit fac'] = float(match.group(1))
            continue

        # Read list of collision energies
        match = re.search(r'CALCULATIONS WILL BE PERFORMED FOR *(\d+) *ENERG', line)
        if match:
            num_energies = int(match.group(1))
            data['num_energies'] = num_energies
            data['energy_unit'] = 'CM-1'
            # Energy list is always printed in units of wavenumbers,
            # see mol.driver.f, format 1250
            energies = []
            read_energies = True
            continue
        if read_energies and len(energies) < num_energies:
            match = re.search(r'ENERGY +\d+ += +(\d+\.\d+) +CM-1', line)
            if not match:
                raise ValueError(f"Problem reading energies from line\n{line}")
            else:
                energies.append(float(match.group(1)))
            if len(energies) == num_energies:
                data['energies'] = energies
                read_energies = False
            continue

        # Read number of J and symmetry blocks
        match = re.search(r'TOTAL ANGULAR MOMENTUM JTOT RUNS FROM +(\d+) +TO +(\d+) +IN STEPS OF +(\d+)', line)
        if match:
            data['JTOTL'] = int(match.group(1))
            data['JTOTU'] = int(match.group(2))
            data['JSTEP'] = int(match.group(3))
            continue
        match = re.search(r'EACH JTOT IS SPLIT INTO A MAXIMUM OF +(\d+) +SYMMETRY BLOCKS', line)
        if match:
            data['max symmetries'] = int(match.group(1))
            continue

    # Stream ended before initialization was complete
    # TODO: handle this case
    return data


def parse_calculation(stream: TextIO | List[str]) -> Dict[str, List[List[str]]]:
    data = {'j_blocks': []}
    block: List[str] = []
    nblocks: int = 0

    calculation_started = False
    for line in stream:
        line = line.strip()
        if not line: continue # skip empty lines

        # Start of calculation
        if not calculation_started:
            # We check for the separator line containing the LABEL.
            # TODO: consider that this line is only printed for IPRINT >= 1
            match = re.search(CALCULATION_PATTERNS['label_line'], line)
            if match:
                data['label'] = match.group(1).strip()
                calculation_started = True
                separator_line = line
                continue

        # End of calculation
        elif calculation_started and line == separator_line:
            break

        # Collect and parse the calculation blocks
        elif calculation_started:

            # Start of a new JTOT/SYMMETRY block
            if CALCULATION_PATTERNS['block_start'].match(line):
                # If we already found a block, finalize it
                if nblocks > 0:
                    data = _parse_j_block(block, data)
                    block.clear()
                block.append(line)

                nblocks += 1
                continue

            # If block started, collect lines
            if nblocks > 0:
                block.append(line)
                continue

    # Finalize last block (if any)
    if nblocks > 0:
        data = _parse_j_block(block, data)
    data['nblocks'] = nblocks
    #TODO: finalize data

    return data


def _parse_j_block(lines: List[str], data: Dict[str, List[List[str]]]) -> Dict[str, List[List[str]]]:
    block_data = {'E_blocks': []}

    # From first line of block, extract J and symmetry
    match = CALCULATION_PATTERNS['block_start'].match(lines[0])
    block_data['Jtot'] = int(match.group(1))
    block_data['symmetry'] = int(match.group(2))

    # From beginning of j block, extract the channels
    channel_block = []
    i_start = None
    i_stop = None
    for i, line in enumerate(lines):
        line = line.strip()
        if not line: continue
        if i_start is None and CALCULATION_PATTERNS['channels'].match(line):
            # found channel list
            i_start = i
        if (
            i_start is not None
            and ENERGY_BLOCK_SECTION_PATTERNS['propagation'].match(line)
        ): # end of channel list
            i_stop = i - 1
            break
    if i_start is not None and i_stop is not None:
        channel_block = lines[i_start:i_stop+1]
    if channel_block:
        block_data = _parse_jblock_channels(channel_block, block_data)

    # Split j block into energy blocks
    j_block = '\n'.join(lines)
    E_blocks = re.split(CALCULATION_PATTERNS['energy_block_separator'], j_block)
    # TODO: check number of energy blocks matches number of collision energies

    for block in E_blocks:
        block_data = _parse_energy_block(block, block_data)
    data['j_blocks'].append(block_data)

    return data # TODO: finalize data


def _parse_jblock_channels(lines, data):
    channels = {}
    quantum_numbers = None
    for line in lines:
        line = line.strip()
        if not line: continue
        if quantum_numbers is None:
            match = re.search(r'^\s*-+\s+(.+)\s+-+', line)
            if match:
                quantum_numbers = match.group(1).split()
                column_header = (
                    ['channel', 'pair state']
                    + quantum_numbers
                    + ['L', 'pair level', 'pair energy']
                )
                channels['column_header'] = column_header
                channels['rows'] = []
                continue
        if quantum_numbers:
            match = re.search(r'^ *REFERENCE ENERGY IS *(\d+\.?\d*) *(\S*) *= *(\d+\.?\d*) *(\S*)', line)
            if match:
                # End of channel list
                channels['ref_energy'] = float(match.group(1))
                channels['ref_energy_unit'] = match.group(2)
                break
            else:
                # Row of channel list
                split = line.split()
                assert len(split) == len(column_header)
                channels['rows'].append(
                    [int(col) for col in split[:-1]] + [float(split[-1])]
                )
    data['channels'] = channels
    return data


def _parse_energy_block(block, data):
    block_data = {}
    sections: Dict[str, List[List[str]]] = {name: [] for name in ENERGY_BLOCK_SECTION_PATTERNS.keys()}
    current_section = None
    buffer: List[str] = []

    def flush():
        if current_section and buffer and not sections[current_section]:
            sections[current_section] = buffer.copy()
        buffer.clear()

    lines = block.split('\n')
    for line in lines:
        line = line.strip()
        if not line: continue # skip empty lines

        # check if this line starts any section
        for name, pattern in ENERGY_BLOCK_SECTION_PATTERNS.items():
            if pattern.match(line):
                # start new section
                flush()
                current_section = name
                break

        # if we have an active section, collect
        if current_section:
            buffer.append(line)

    # flush last section
    flush()

    if sections.get('propagation'):
        block_data = _parse_Eblock_propagation(sections['propagation'], block_data)
    if sections.get('open_channels'):
        block_data = _parse_Eblock_open_channels(sections['open_channels'], block_data)
    if sections.get('s_matrix'):
        block_data = _parse_Eblock_s_matrix(sections['s_matrix'], block_data)
    if sections.get('pcs'):
        block_data = _parse_Eblock_partial_cross_sections(sections['pcs'], block_data)
    if sections.get('ics'):
        block_data = _parse_Eblock_integrated_cross_sections(sections['ics'], block_data)

    # TODO: make sure, parsing is robust against unknown sections

    data['E_blocks'].append(block_data)

    return data


def _parse_Eblock_propagation(lines, data):

    # Extract collision energy from first line
    split = lines[0].split()
    data['energy'] = float(split[6])
    data['energy_unit'] = split[7]

    # Extract absolute energy from second line
    match = re.search(r'ABSOLUTE ENERGY IS +(\d+\.\d+) +CM-1', lines[1])
    if match:
        data['abs_energy'] = float(match.group(1))

    #TODO: add fallback, if not found
    return data


def _parse_Eblock_open_channels(lines, data):
    channels = {}
    column_header = None
    for line in lines:
        line = line.strip()
        if not line: continue
        if line.startswith('OPEN CHANNEL'):
            column_header = re.split(r' {2,}', line)
            channels['column_header'] = column_header
            channels['rows'] = []
        elif column_header:
            split = line.split()
            assert len(split) == len(column_header)
            channels['rows'].append([
                int(split[0]),
                float(split[1]),
                int(split[2]),
                int(split[3]),
                int(split[4]),
                float(split[5]),
            ])
        else:
            raise ValueError("No open channels found")
    data['open channels'] = channels
    return data


def _parse_Eblock_s_matrix(lines, data):
    s_matrix = {}
    column_header = None
    for line in lines:
        line = line.strip()
        if not line: continue
        if line.startswith('ROW  COL'):
            column_header = re.split(r' {2,}', line)
            s_matrix['column_header'] = column_header
            s_matrix['rows'] = []
        elif column_header:
            split = line.split()
            assert len(split) == len(column_header)
            s_matrix['rows'].append(
                [int(split[0]), int(split[1])] + [float(col) for col in split[2:]]
            )
        else:
            raise ValueError("No S matrix found")
        #TODO: make sure, number of lines corresponds to number of open channels
    data['S'] = s_matrix
    return data


def _parse_Eblock_partial_cross_sections(lines, data):
    cs_data = {}
    # From second line, extract JTOT, SYMMETRY and ENERGY
    line = lines[1]
    match = re.search(r'^ *FOR JTOT = *(\d+) *AND SYMMETRY BLOCK = *(\d+) * AT ENERGY\( *(\d+)\) = *(\d+\.\d+) *(\S+)', line)
    cs_data['jtot'] = int(match.group(1))
    cs_data['sym'] = int(match.group(2))
    cs_data['iE'] = int(match.group(3))
    cs_data['E'] = float(match.group(4))
    cs_data['E unit'] = match.group(5)

    cs_data['data'] = _parse_cross_sections(lines)
    data['pcs'] = cs_data
    return data


def _parse_Eblock_integrated_cross_sections(lines, data):
    cs_data = {}
    # From first line, extract accumulated JTOT values
    line = lines[0]
    match = re.search(r'ACCUMULATED FROM JTOT = *(\d+) *TO *(\d+)', line)
    cs_data['jtot min'] = int(match.group(1))
    cs_data['jtot max'] = int(match.group(2))

    cs_data['data'] = _parse_cross_sections(lines)
    data['acs'] = cs_data
    return data


def _parse_cross_sections(lines):
    data = {}
    # Find the header lines for columns (lines starting with 'F  I =')
    cross_section_blocks = []
    block = []
    header_indices = []
    for idx, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
        if line.startswith('F  I ='):
            if block:
                cross_section_blocks.append((header_indices[-1], block))
                block = []
            header_indices.append(idx)
            continue
        if line[0].isdigit():
            block.append(line)
    if block:
        cross_section_blocks.append((header_indices[-1], block))

    # Now parse the blocks into a matrix
    # Find the maximum F and I indices to determine the matrix size
    max_f = 0
    max_i = 0
    for header_idx, block in cross_section_blocks:
        header = lines[header_idx]
        col_indices = [int(x) for x in header.split('=')[1].split()]
        max_i = max(max_i, max(col_indices))
        for row in block:
            parts = row.split()
            f_idx = int(parts[0])
            max_f = max(max_f, f_idx)
    #TODO: max_f and max_i should probably be equal

    # Populate the matrix
    #TODO: As a sanity check, we can use that the maximum value max_f and max_i
    #      is the number of levels. This can maybe also used to simply parsing.
    cross_section = np.zeros((max_f, max_i))
    for header_idx, block in cross_section_blocks:
        header = lines[header_idx]
        col_indices = [int(x) for x in header.split('=')[1].split()]
        for row in block:
            parts = row.split()
            f_idx = int(parts[0]) - 1  # zero-based index
            for c, val in enumerate(parts[1:]):
                i_idx = col_indices[c] - 1  # zero-based index
                try:
                    cross_section[f_idx, i_idx] = float(val)
                except ValueError:
                    # This might occur because a broken Fortran exponent,
                    # e.g., 7.65245-104
                    if re.match(r'[+-]\d{3}', val[-4:]):
                        new_val = val[0:-4] + 'E' + val[-4:]
                        cross_section[f_idx, i_idx] = float(new_val)

    data['values'] = cross_section
    data['I'] = [i+1 for i in range(max_i)]
    data['F'] = [i+1 for i in range(max_f)]
    return data


def parse_summary(stream: TextIO | List[str]) -> Dict[str, List[List[str]]]:
    data = {}
    sections: Dict[str, List[List[str]]] = {name: [] for name in SUMMARY_SECTION_PATTERNS.keys()}
    current_section = None
    buffer: List[str] = []

    def flush():
        if current_section and buffer and not sections[current_section]:
            sections[current_section] = buffer.copy()
        buffer.clear()

    for line in stream:
        line = line.strip()
        if not line: continue # skip empty lines

        # check if this line starts any section
        for name, pattern in SUMMARY_SECTION_PATTERNS.items():
            if pattern.match(line):
                # start new section
                flush()
                current_section = name
                break

        # if we have an active section, collect
        if current_section:
            buffer.append(line)

    # flush last section
    flush()

    if sections.get('ics'):
        data = _parse_summary_ics(sections['ics'], data)
    if sections.get('footer'):
        data = _parse_summary_footer(sections['footer'], data)

    return data


def _parse_summary_ics(lines: List[str], data: Dict[str, List[List[str]]]) -> Dict[str, List[List[str]]]:
    re_separator = re.compile(r'ENERGY \(CM-1\)\s+JTOTL\s+JSTEP\s+JTOTU\s+F\s+I\s+SIG\(F,I\)')
    section = '\n'.join(lines)
    blocks = re.split(re_separator, section)

    # First block contains threshold energies
    block = blocks.pop(0)
    threshold_energies = np.array(
        [float(line.split()[1]) for line in block.strip().split("\n")[1:]]
        # Note: we have to skip the first line, which is the separator line
    )
    data['threshold_energies'] = threshold_energies

    # Remaining blocks contain the integral cross sections for each collision energy
    ics = []
    total_inelastic = []
    for block in blocks:
        i_inel = []
        sigma_inel = []
        f_idx = []
        i_idx = []
        sigma = []
        _lines = block.strip().split('\n')
        for line in _lines:
            split = line.split()
            if len(split) == 7:
                f_idx.append(int(split[4]))
                i_idx.append(int(split[5]))
                sigma.append(float(split[6]))
            elif len(split) == 2:
                sigma_inel.append(float(split[0]))
                i_inel.append(int(split[1]))
        ics.append([f_idx, i_idx, sigma])
        total_inelastic.append([i_inel, sigma_inel])
    data['ics'] = ics
    data['total_inelastic'] = total_inelastic

    return data


def _parse_summary_footer(lines: List[str], data: Dict[str, List[List[str]]]) -> Dict[str, List[List[str]]]:
    footer = '\n'.join(lines)
    data['time'] = float(re.search(r"This run used\s+(\d+\.\d+) cpu secs", footer).group(1))
    words = int(re.search(r"\s+([\d,]+)\s+of the allocated", footer).group(1))
    data['memory'] = words * 8 / 1e6 # convert to MB
    return data


def _process_result_data(data):

    para = MolscatInputParameters(
        IPRINT = data['init'].get('IPRINT'),
        JTOTL = data['init'].get('JTOTL'),
        JTOTU = data['init'].get('JTOTU'),
        JSTEP = data['init'].get('JSTEP'),
    )

    result = MolscatResult(
        para = para,
        type = data['init'].get('interaction type'),
        runtime = data['summary'].get('time'),
        memory = data['summary'].get('memory'),
        mu = data['init'].get('reduced mass'),
        energies = data['init'].get('energies'),
        _pair_state_list = data['init'].get('pair_states'),
        symmetries = data['init'].get('max symmetries'),
        blocks = _process_blocks(data),
        ics = _process_total_integrated_cross_sections(data),
        acs = _process_accumulated_cross_sections(data),
    )

    return result


def _process_total_integrated_cross_sections(data):
    """
    Read total integrated cross section from data dict and store as DataArray.
    """
    ics_data = data['summary'].get('ics')
    if ics_data is None:
        return None
    else:
        # All ics_data entries have the same structure: [f_idx, i_idx, sigma]
        # We'll use the first entry to get the unique indices
        f_indices = np.unique(np.concatenate([np.array(entry[0]) for entry in ics_data]))
        i_indices = np.unique(np.concatenate([np.array(entry[1]) for entry in ics_data]))
        energies = data['init'].get('energies')

        # Prepare a 3D array: (energy, f, i)
        cross_sections = np.zeros((len(ics_data), len(f_indices), len(i_indices)))

        for e_idx, entry in enumerate(ics_data):
            f_idx_arr = np.array(entry[0], dtype=int) - 1  # zero-based
            i_idx_arr = np.array(entry[1], dtype=int) - 1
            sigma_arr = np.array(entry[2], dtype=float)
            cross_sections[e_idx, f_idx_arr, i_idx_arr] = sigma_arr

        return xr.DataArray(
            cross_sections,
            dims=["E", "f", "i"],
            coords={
                "E": energies,
                "f": f_indices,
                "i": i_indices,
            },
            name="cross_section"
        )


def _process_accumulated_cross_sections(data):
    acs_list = [[]]
    for block_data in data['calc'].get('j_blocks', []):
        for E_block in block_data.get('E_blocks', []):
            acs_dict = E_block.get('acs')
            if acs_dict is not None:
                # Extract one cross section
                acs_data = acs_dict.pop('data')
                f_indices = acs_data['F']
                i_indices = acs_data['I']
                values = acs_data['values']
                acs = xr.DataArray(
                    values,
                    dims=["f", "i"],
                    coords={
                        "f": f_indices,
                        "i": i_indices,
                    },
                    name="partial cross_section"
                )
                acs.attrs = acs_dict

                # Store cross section at correct place in acs_list
                # Ensure acs_list has enough sublists for each energy
                energy = float(E_block.get('energy'))
                energy_idx = np.argmin(abs(np.asarray(data['init']['energies']) - energy))
                while len(acs_list) <= energy_idx:
                    acs_list.append([])
                acs_list[energy_idx].append(acs)
    return acs_list


def _process_blocks(data):
    blocks = []
    for block_data in data['calc'].get('j_blocks', []):

        channels = block_data.get('channels', {}).get('rows', [])

        # Extract data from energy sub-blocks
        energies = []
        n_open = []
        S = []
        pcs = []
        for E_block in block_data.get('E_blocks', []):
            energies.append(E_block.get('energy'))
            n_open.append(len(E_block.get('open channels', {}).get('rows', [])))
            S.append(_process_S_matrix(E_block, n_open[-1]))
            pcs.append(_process_partial_cross_sections(E_block, n_open[-1]))

        # Fill J-block data structure
        blocks.append(ScatteringBlock(
            Jtot = block_data.get('Jtot'),
            sym = block_data.get('symmetry'),
            energies = energies,
            channels = channels,
            n_channels = len(channels),
            n_open = n_open,
            S = S,
            pcs = pcs,
        ))
    return blocks


def _process_S_matrix(E_block, n_open):
    S_dict = E_block.get('S')
    if S_dict is not None:
        S_array = np.asarray(S_dict['rows'], dtype=object)
        S = np.zeros((n_open, n_open), dtype=complex)
        for i in range(len(S_array)):
            row, col, S_sqr, S_ang, S_re, S_im = S_array[i]
            S_val = S_re + 1j*S_im
            assert np.allclose(S_sqr, np.abs(S_val)**2), (
                f"{S_sqr} {np.abs(S_val)**2} {S_sqr-np.abs(S_val)**2}"
            )
            S[row-1, col-1] = S_val
    else:
        S = None
    return S


def _process_partial_cross_sections(E_block, n_open):
    pcs_dict = E_block.get('pcs')
    if pcs_dict is not None:
        pcs_data = pcs_dict.pop('data')
        f_indices = pcs_data['F']
        i_indices = pcs_data['I']
        values = pcs_data['values']
        pcs = xr.DataArray(
            values,
            dims=["f", "i"],
            coords={
                "f": f_indices,
                "i": i_indices,
            },
            name="partial cross_section"
        )
        pcs.attrs = pcs_dict
    else:
        pcs = None
    return pcs