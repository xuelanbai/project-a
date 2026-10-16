import contextlib
import gzip
import io
import json
from pathlib import Path
import sys
import tempfile

script = Path(__file__).with_name('collect_cuttag_v2_stats.slurm').read_text(encoding='utf-8')
code = script.split("<<'PY'\n", 1)[1].rsplit('\nPY', 1)[0]
namespace = {'__name__': 'collector_tests'}
exec(compile(code, 'collector_embedded.py', 'exec'), namespace)

def rec(start, end):
    return ('1', start, end, [])

shared = namespace['shared_intervals']
support = namespace['supported_count']
assert shared([rec(0,10)], [rec(10,20)]) == []
assert shared([rec(0,10)], [rec(9,20)]) == [('1',9,10)]
assert shared([rec(0,10), rec(10,20)], [rec(9,11)]) == [('1',9,11)]
assert support([rec(0,10), rec(10,20)], [rec(9,11)]) == 2
assert support([rec(0,10)], []) == 0
assert shared([], [rec(0,10)]) == []
# Duplicate source peaks must not duplicate common genomic intervals.
assert shared([rec(0,10),rec(0,10)], [rec(2,8)]) == [('1',2,8)]

with tempfile.TemporaryDirectory(prefix='cuttag_stats_test_') as directory:
    root = Path(directory)
    (root/'config').mkdir()
    (root/'config/samplesheet.csv').write_text('group,replicate,fastq_1,fastq_2,control\n', encoding='utf-8')
    results = root/'results_pig_v2'
    for caller in ['macs2','seacr']:
        (results/'03_peak_calling/04_called_peaks'/caller).mkdir(parents=True)
    peak_dir = results/'03_peak_calling/04_called_peaks'
    for state in ['119','137']:
        for mark in ['Pol2','Tag']:
            for rep in [1,2]:
                prefix = '3D4_{}_{}_R{}'.format(state,mark,rep)
                start,end = (0,10) if rep==1 else (9,20)
                (peak_dir/'macs2'/(prefix+'.macs2_peaks.narrowPeak')).write_text(
                    '1\t{}\t{}\tpeak\t100\t.\t10\t5\t2\t0\n'.format(start,end),encoding='utf-8')
                (peak_dir/'seacr'/(prefix+'.seacr.peaks.stringent.bed')).write_text('',encoding='utf-8')
    # Include a present and a missing final pooled track.
    tracks = results/'custom_tracks/macs2_ppois'
    tracks.mkdir(parents=True)
    (tracks/'3D4_119_Pol2_rep1-rep2_pooled_vs_3D4_119_IgG_rep1.MACS2_ppois_Sscrofa11.1.bigWig').write_bytes(b'test')
    sys.argv = ['collector', str(root)]
    with contextlib.redirect_stdout(io.StringIO()):
        namespace['main']()
    with gzip.open(root/'logs/cuttag_v2_table_inputs.json.gz','rt',encoding='utf-8') as handle:
        data=json.load(handle)
    assert len(data['peak_files'])==16
    assert len(data['reproducibility'])==12
    assert not data['errors'] and not data['missing']
    for row in data['reproducibility']:
        assert row['shared_intervals'] == (0 if row['caller']=='SEACR' else 1)
    assert sum(t['exists'] for t in data['merged_tracks'])==1
    assert all(t['bytes'] is None for t in data['merged_tracks'] if not t['exists'])
print('PASS: overlap >=1 bp, adjacency excluded, empty peak sets, duplicate intervals, P/Q filtering and end-to-end JSON export.')
