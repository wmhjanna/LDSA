from pathlib import Path
import importlib.util, shutil
ROOT = Path('@PROJECT_ROOT@')
W = Path(__file__).parent
MAIN = W.parent
assert (MAIN / 'workbook_done').exists(), 'Final workbook analysis still running'
spec = importlib.util.spec_from_file_location('qc_table_pairs', MAIN / 'workbook/redraw.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
m.WORK = W
m.tr = m.tr[~m.tr.Site.isin(['PAR_SUB', 'PRA_SUB'])].copy()
m.ms = m.ms[~m.ms.Site.isin(['PAR_SUB', 'PRA_SUB'])].copy()
m.pair_tables()
for name in ['exact_overlap_pairs', 'classification', 'exact_overlap_breakpoints_min12', 'correlation_temporal_sensitivity']:
    shutil.copy2(W / (name + '.csv'), m.ANA / (name + '.csv'))
print('Paired-comparison tables prepared for document output')
