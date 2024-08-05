from pathlib import Path

cycle32 = Path('cycle 32')
input32 = Path(cycle32 / 'input_catalogs')
output32 = Path(cycle32 / 'output_catalogs')
etc32 = Path(cycle32 / 'ETC Grids')

reffolder = Path('reference_tables')
refspecfolder = Path('reference_uv_spectra')

stage1_targets = input32 / 'stage1_targets_proposed.ecsv'
preliminary_targets = input32 / 'preliminary_targets_2024-02.ecsv'
toi_false_positives = input32 / 'toi_false_positives_from_shreyas.txt'