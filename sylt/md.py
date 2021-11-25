import pandas as pd
import numpy as np


def load_md_data(filename, skip_rows=98):
    """load & parse md data files and return tau, turn, V arrays"""

    with open(filename, 'r') as file:
        header = [next(file).strip() for _ in range(skip_rows)]

    V_scope = pd.read_csv(filename, header=skip_rows-1)
    V_scope = np.array(V_scope)[:, 0]

    data = {
        'epoch': int(header[1]),
        'n_frames': int(header[16]),
        'n_bins': int(header[20]),
        'bin_width': float(header[22]),
        'turns_frame': float(header[24]),
        'Vg': float(header[61]),
        'h': int(header[69]),
        'B': float(header[75]),
        'dBdt': float(header[77]),
        'R': float(header[79]),
        'rho': float(header[81]),
        'gamma_t': float(header[83]),
        'E0': float(header[85]),
        'q': int(header[87]),
        'pickup_sensitivity': float(header[97]),
    }

    data['tau'] = np.linspace(-0.5, 0.5, data['n_bins'])*data['bin_width']*data['n_bins']
    data['turns'] = np.arange(data['n_frames'])*data['turns_frame']
    data['V'] = np.reshape(V_scope, (data['n_frames'], data['n_bins']))

    return data
