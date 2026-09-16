# usage: python dev/latex/list_unused_glossary.py
# lists items that appear in glossary but are not used in papers

import re
from pathlib import Path

key_pat = re.compile(r'\\newsym\s*\{\s*([^{}]+?)\s*\}')
ref_pat = re.compile(r'\\gls\s*\{\s*(sym-[^}]+?)\s*\}')

roots = [Path('pubs'), Path('mdocean/simulation/modules/OpenFLASH/pubs/JFM')]

keys = set()
for p in Path('pubs/shared/glossary').glob('*.tex'):
    text = p.read_text(encoding='utf-8', errors='ignore')
    keys |= {m.group(1).strip() for m in key_pat.finditer(text)}

used = set()
for root in roots:
    if not root.exists():
        continue
    for p in root.rglob('*.tex'):
        text = p.read_text(encoding='utf-8', errors='ignore')
        for m in ref_pat.finditer(text):
            ref = m.group(1)
            if ref.startswith('sym-'):
                key = ref[4:]
                if key in keys:
                    used.add(key)

unused = sorted(keys - used)
print(f'UNUSED_KEYS={len(unused)}')
for key in unused:
    print(key)