#!/usr/bin/env python3

import re
from pathlib import Path

WORDS = [
    'array', 'buckle', 'circ', 'clear', 'constr', 'heave', 'mult',
    'cost', 'struct', 'crit', 'tube', 'plate', 'top', 'edge', 'cent',
    'fixed', 'constant', 'max', 'min', 'elec', 'drag', 'up', 'down',
    'linear', 'slam', 'wave', 'oval', 'unconstrained', 'shaft',
    'harmonics', 'limits', 'shape', 'stiff', 'rel', 'damp',
    'avg', 'surf', 'sub', 'guess', 'mech', 'reac', 'storm', 
    'design', 'sens', 'spar', 'unsat',
]

word_pat = re.compile(
    r'(?<!\\text\{)\b(' + '|'.join(map(re.escape, WORDS)) + r')\b'
)

math_pat = re.compile(
    r'''
    (?<!\\)\$\$(?:\\.|[^\\])*?(?<!\\)\$\$         |   # $$ ... $$
    (?<!\\)\$(?!\$)(?:\\.|[^\\$])*?(?<!\\)\$(?!\$) |   # $ ... $
    \\\((?:\\.|[^\\)])*?\\\)                      |   # \( ... \)
    \\\[(?:\\.|[^\\\]])*?\\\]                      # \[ ... \]
    ''',
    re.DOTALL | re.VERBOSE,
)


def replace_in_subscripts(math):
    # {...} subscripts
    def brace_repl(m):
        contents = m.group(1)

        def wrap_word(w):
            if w.group(0).startswith(r'\text{'):
                return w.group(0)
            return rf'\text{{{w.group(1)}}}'

        new_contents = word_pat.sub(wrap_word, contents)
        return '_{' + new_contents + '}'

    math = re.sub(
        r'_\{([^{}]*)\}',
        brace_repl,
        math,
    )

    # Single-token subscripts: _max -> _{\text{max}}
    def bare_repl(m):
        word = m.group(1)
        return rf'_{{\text{{{word}}}}}'

    math = re.sub(
        r'_(?:(' + '|'.join(map(re.escape, WORDS)) + r'))\b',
        bare_repl,
        math,
    )

    return math


def process_text(text):
    def math_repl(m):
        return replace_in_subscripts(m.group(0))

    return math_pat.sub(math_repl, text)


for path in Path("pubs").rglob("*.tex"):
    old = path.read_text(encoding="utf-8")
    new = process_text(old)

    if old != new:
        path.write_text(new, encoding="utf-8")
        print(f"Modified {path}")