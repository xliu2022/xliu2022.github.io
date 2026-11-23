"""Reproduce the exact finite result using Python 3 and a C++17 compiler."""
from pathlib import Path
import json, subprocess, tempfile, shutil

ROOT = Path(__file__).resolve().parent
compiler = shutil.which('g++') or shutil.which('clang++')
if compiler is None:
    raise SystemExit('A C++17 compiler (g++ or clang++) is required.')

subprocess.run(['python3', str(ROOT/'universal_nine_missing_orbits.py')],
               check=True, stdout=subprocess.DEVNULL)
subprocess.run(['python3', str(ROOT/'root_finite_proof_audit.py')], check=True)

outputs = []
with tempfile.TemporaryDirectory(prefix='rigid_nine_check_') as tmp:
    for stem in ['universal_nine_exact_certificate', 'universal_nine_integer_recheck']:
        exe = Path(tmp)/stem
        subprocess.run([compiler, '-std=c++17', '-O3', str(ROOT/(stem+'.cpp')),
                        '-o', str(exe)], check=True)
        p = subprocess.run([str(exe), str(ROOT/'universal_nine_missing_orbits.txt')],
                           check=True, capture_output=True, text=True)
        output = [json.loads(line) for line in p.stdout.splitlines() if line.strip()]
        outputs.append([row for row in output if not row.get('summary')])
        (ROOT/(stem+'_rerun.jsonl')).write_text(p.stdout)

a, b = outputs
assert len(a) == len(b) == 20
for first, second in zip(a, b):
    for key in ['missing_mask', 'orbit_size', 'triangle_compatible',
                'parity_compatible', 'residual_covered', 'enough_free_edges']:
        assert first[key] == second[key], (key, first, second)
    assert second['enough_free_edges'] == 0
result = {
    'status': 'PASS', 'missing_pair_orbits': 20,
    'balanced_orientations': sum(r['all_balanced_orientations'] for r in b),
    'triangle_compatible': sum(r['triangle_compatible'] for r in b),
    'parity_compatible': sum(r['parity_compatible'] for r in b),
    'residual_covered': sum(r['residual_covered'] for r in b),
    'final_survivors': sum(r['enough_free_edges'] for r in b),
    'computer_assisted_conclusion': 'f_3(9)=13',
}
assert result['balanced_orientations'] == 19427814
assert result['triangle_compatible'] == 1584297
assert result['residual_covered'] == 716
(ROOT/'finite_reproduction_result.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
