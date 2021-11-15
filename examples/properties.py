ring_properties = {
    'R': 100,              # machine radius (m)
    'h': 8,                # rf harmonic
    'Vg': 80e3,            # gap voltage (V)
    'gamma_t': 6.1,        # transition

    'beta': 17,            # nominal beta function
    'D': (2.66, 0),        # nominal dispersion function
    'b': (73e-3, 35e-3),   # aperture radius (m)
}

bunch_properties = {
    'E': 2.938272e9,        # beam energy (eV)
    'sig_w': 5e6,           # energy spread (eV)
    'sig_tau': 30e-9,       # bunch length (s)
    'sig_eps': 1e-6,        # transverse emittance (m)
}
