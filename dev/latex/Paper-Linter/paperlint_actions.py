import re
from pathlib import Path

import glossary_symbols as gs


REF_CMDS = ["ref", "eqref", "cref", "Cref", "autoref", "pageref", "nameref"]

WORDS = [
    "array", "buckle", "circ", "clear", "constr", "heave", "mult",
    "cost", "struct", "crit", "tube", "plate", "top", "edge", "cent",
    "fixed", "constant", "max", "min", "elec", "drag", "up", "down",
    "linear", "slam", "wave", "oval", "unconstrained", "shaft",
    "harmonics", "limits", "shape", "stiff", "rel", "damp", "opt",
    "avg", "surf", "sub", "guess", "mech", "reac", "storm",
    "design", "sens", "spar", "unsat", "in", "rad"
]

word_pat = re.compile(r'(?<!\\text\{)\b(' + '|'.join(map(re.escape, WORDS)) + r')\b')

math_pat = re.compile(
    r'''
    (?<!\\)\$\$(?:\\.|[^\\])*?(?<!\\)\$\$         |   # $$ ... $$
    (?<!\\)\$(?!\$)(?:\\.|[^\\$])*?(?<!\\)\$(?!\$) |   # $ ... $
    \\\((?:\\.|[^\\)])*?\\\)                      |   # \( ... \)
    \\\[(?:\\.|[^\\\]])*?\\\]                      # \[ ... \]
    ''',
    re.DOTALL | re.VERBOSE,
)

WORD_STYLE_REFERENCE_TEXT = ""

TERMINOLOGY_REPLACEMENTS = [
    (re.compile(r"\bexternal region(s)?\b", re.IGNORECASE), "exterior region"),
    (re.compile(r"\binternal region(s)?\b", re.IGNORECASE), "interior region"),
    (re.compile(r"\breaction plate(s)?\b", re.IGNORECASE), "damping plate"),
]


def _read_text(path):
    return Path(path).read_text(encoding="utf-8")


def _write_text(path, text):
    Path(path).write_text(text, encoding="utf-8")


def _rewrite_lines(path, transform):
    path = Path(path)
    old = _read_text(path)
    new = "".join(transform(line) for line in old.splitlines(keepends=True))
    if old != new:
        _write_text(path, new)
        return 1
    return 0


def replace_in_subscripts(math):
    def brace_repl(match):
        contents = match.group(1)

        def wrap_word(word_match):
            if word_match.group(0).startswith(r'\text{'):
                return word_match.group(0)
            return rf'\text{{{word_match.group(1)}}}'

        return '_{' + word_pat.sub(wrap_word, contents) + '}'

    math = re.sub(r'_\{([^{}]*)\}', brace_repl, math)

    def bare_repl(m):
        word = m.group(1)
        return rf'_{{\text{{{word}}}}}'

    return re.sub(r'_(?:(' + '|'.join(map(re.escape, WORDS)) + r'))\b', bare_repl, math)


def replace_mathmode_subscripts(text):
    return math_pat.sub(lambda match: replace_in_subscripts(match.group(0)), text)


def fix_mathmode_subscripts_in_file(path):
    path = Path(path)
    old = _read_text(path)
    new = replace_mathmode_subscripts(old)
    if old != new:
        _write_text(path, new)
        return 1
    return 0


def space_before_cite_line(line):
    match = re.search(r'[^ ~]\\cite(?![A-Za-z])', line)
    if not match or '\\etal\\cite' in line:
        return line, None
    start = match.start(0)
    return line[:start + 1] + ' ' + line[start + 1:], match.span(0)


def fix_space_before_cite_in_file(path):
    return _rewrite_lines(path, lambda line: space_before_cite_line(line)[0])


def preferred_terminology_line(line):
    updated = line
    changed = False

    for pattern, replacement in TERMINOLOGY_REPLACEMENTS:
        def replace(match, replacement=replacement):
            plural_suffix = match.group(1) or ""
            term = replacement + plural_suffix
            if match.group(0) and match.group(0)[0].isupper():
                term = term[0].upper() + term[1:]
            return term

        new_updated = pattern.sub(replace, updated)
        changed |= new_updated != updated
        updated = new_updated

    if not changed:
        return line, None
    return updated, (0, len(line))


def fix_preferred_terminology_in_file(path):
    return _rewrite_lines(path, lambda line: preferred_terminology_line(line)[0])


def collect_unused_macro_definitions(text):
    definitions = []
    patterns = [
        re.compile(r'\\newcommand\*?\{\\([A-Za-z@]+)\}'),
        re.compile(r'\\providecommand\*?\{\\([A-Za-z@]+)\}'),
        re.compile(r'\\def\\([A-Za-z@]+)'),
        re.compile(r'\\gdef\\([A-Za-z@]+)'),
        re.compile(r'\\xdef\\([A-Za-z@]+)'),
    ]
    for line_no, line in enumerate(text.splitlines()):
        for pattern in patterns:
            match = pattern.search(line)
            if not match:
                continue
            name = match.group(1)
            if len(re.findall(r'\\%s\b' % re.escape(name), text)) <= 1:
                definitions.append((name, line_no, match.span(1)))
            break
    return definitions


def fix_unused_macros_in_file(path):
    path = Path(path)
    old = _read_text(path)
    unused = {name for name, _, _ in collect_unused_macro_definitions(old)}
    if not unused:
        return 0

    def transform(line):
        updated = line
        for name in unused:
            updated = re.sub(
                r'\\(?:newcommand|providecommand)\*?\{\\%s\}(?:\[[^\]]*\])*(?:\{[^{}]*\})?' % re.escape(name),
                '',
                updated,
            )
            updated = re.sub(r'\\(?:gdef|xdef|def)\\%s\b[^\n]*' % re.escape(name), '', updated)
        return updated

    return _rewrite_lines(path, transform)


def duplicate_word_line(line):
    match = re.search(r'\b(\w+)\s+\1\b', line, re.IGNORECASE)
    if not match:
        return line, None
    word = match.group(1)
    start, end = match.span()
    return line[:start] + word + line[end:], match.span()


def fix_duplicate_words_in_file(path):
    return _rewrite_lines(path, lambda line: duplicate_word_line(line)[0])


def comment_has_space_line(line):
    ls = line.strip()
    if '%' not in ls or ls.endswith('%') or ls.startswith('%'):
        return line, None
    match = re.search(r'[^\s\\\}\{%]+%', line)
    if not match:
        return line, None
    return re.sub(r'([^\s\\\}\{%])%', r'\1 %', line, count=1), match.span()


def fix_comment_has_space_in_file(path):
    return _rewrite_lines(path, lambda line: comment_has_space_line(line)[0])


def section_capitalization_line(line):
    match = re.search(r'(\\section\{)([^}]*)\}(?=\s|$)', line)
    if not match:
        return line, None
    title = match.group(2)

    def capitalize_word(word_match):
        word = word_match.group(0)
        return word[0].upper() + word[1:]

    updated_title = re.sub(r'\b[a-z][A-Za-z\'\-]{4,}\b', capitalize_word, title)
    if updated_title == title:
        return line, None
    return line[:match.start(2)] + updated_title + line[match.end(2):], match.span(2)


def fix_section_capitalization_in_file(path):
    return _rewrite_lines(path, lambda line: section_capitalization_line(line)[0])


def quotation_marks_line(line):
    if re.match(r'\s*\\(?:g|e)?def\b', line):
        return line, None
    comment = re.search(r'(?<!\\)%', line)
    text = line[:comment.start()] if comment else line
    tail = line[comment.start():] if comment else ''
    if '"' not in text:
        return line, None
    result = []
    open_quote = True
    changed = False
    for idx, char in enumerate(text):
        if char == '"' and (idx == 0 or text[idx - 1] != '\\'):
            result.append('``' if open_quote else "''")
            open_quote = not open_quote
            changed = True
        else:
            result.append(char)
    if not changed:
        return line, None
    first_quote = text.find('"')
    return ''.join(result) + tail, (first_quote, first_quote + 1)


def fix_quotation_marks_in_file(path):
    return _rewrite_lines(path, lambda line: quotation_marks_line(line)[0])


def space_before_punctuation_line(line):
    match = re.search(r'\s+([,.!?:;])', line)
    if not match:
        return line, None
    return re.sub(r'\s+([,.!?:;])', r'\1', line, count=1), match.span()


def fix_space_before_punctuation_in_file(path):
    return _rewrite_lines(path, lambda line: space_before_punctuation_line(line)[0])


def multicite_line(line):
    match = re.search(r'(?P<cmd>\\citeA?)\{(?P<left>[^}]*)\}\s*(?P=cmd)\{(?P<right>[^}]*)\}', line)
    if not match:
        return line, None
    cmd = match.group('cmd')
    left = match.group('left').strip()
    right = match.group('right').strip()
    return line[:match.start()] + f'{cmd}{{{left},{right}}}' + line[match.end():], match.span()


def fix_multicite_in_file(path):
    return _rewrite_lines(path, lambda line: multicite_line(line)[0])


def cite_duplicate_line(line):
    match = re.search(r'(?P<cmd>\\(?:no)?citeA?)\{(?P<keys>[^}]*)\}', line)
    if not match:
        return line, None
    keys = [key.strip() for key in match.group('keys').split(',') if key.strip()]
    deduped = []
    seen = set()
    for key in keys:
        if key in seen:
            continue
        seen.add(key)
        deduped.append(key)
    updated = line[:match.start()] + f"{match.group('cmd')}{{{','.join(deduped)}}}" + line[match.end():]
    if updated == line:
        return line, None
    return updated, match.span()


def fix_cite_duplicate_in_file(path):
    return _rewrite_lines(path, lambda line: cite_duplicate_line(line)[0])


def brackets_space_line(line):
    updated = line
    changed = False
    new_updated = re.sub(r'(?<![\\\s\{\[~])\((?=[A-Za-z\\])', ' (', updated)
    changed |= new_updated != updated
    updated = new_updated
    new_updated = re.sub(r'\(\s+', '(', updated)
    changed |= new_updated != updated
    updated = new_updated
    new_updated = re.sub(r'\s+\)', ')', updated)
    changed |= new_updated != updated
    updated = new_updated
    if not changed:
        return line, None
    return updated, (0, len(line))


def fix_brackets_space_in_file(path):
    return _rewrite_lines(path, lambda line: brackets_space_line(line)[0])


def numerals_line(line):
    replacement_pairs = [
        (r'\bthree\b', '3'),
        (r'\bfour\b', '4'),
        (r'\bfive\b', '5'),
        (r'\bsix\b', '6'),
        (r'\bseven\b', '7'),
        (r'\beight\b', '8'),
        (r'\bnine\b', '9'),
        (r'\bten\b', '10'),
        (r'\beleven\b', '11'),
        (r'\btwelve\b', '12'),
    ]
    updated = line
    changed = False
    for pattern, numeral in replacement_pairs:
        new_updated = re.sub(pattern, numeral, updated)
        changed |= new_updated != updated
        updated = new_updated
    if not changed:
        return line, None
    return updated, (0, len(line))


def fix_numerals_in_file(path):
    return _rewrite_lines(path, lambda line: numerals_line(line)[0])


def collect_word_style_data(text):
    text = text + ("\n" + WORD_STYLE_REFERENCE_TEXT if WORD_STYLE_REFERENCE_TEXT else "")
    styled_word_pattern = re.compile(r'\\text([^\{]*)\{([^{}]+)\}')
    style_for_word = {}
    style_counts_by_word = {}
    style_counts = {}
    for styled in styled_word_pattern.finditer(text):
        style = styled.group(1)
        word = styled.group(2)
        style_counts[word] = style_counts.get(word, 0) + 1
        if word not in style_counts_by_word:
            style_counts_by_word[word] = {}
        style_counts_by_word[word][style] = style_counts_by_word[word].get(style, 0) + 1

    for word, counts in style_counts_by_word.items():
        # Prefer the most common style seen for this word.
        style_for_word[word] = max(counts.items(), key=lambda item: item[1])[0]
    return style_for_word, style_counts, styled_word_pattern


def missing_word_style_line(line, style_for_word, style_counts, styled_word_pattern):
    placeholders = []

    def mask_gls(match):
        placeholders.append(match.group(0))
        return f'__PAPERLINT_GLS_{len(placeholders) - 1}__'

    def mask_styled(match):
        placeholders.append(match.group(0))
        return f'__PAPERLINT_STYLE_{len(placeholders) - 1}__'

    updated = re.sub(r'\\gls[a-zA-Z]*\*?\{[^{}]*\}', mask_gls, line)
    updated = styled_word_pattern.sub(mask_styled, updated)
    changed = updated != line
    for word, style in style_for_word.items():
        new_updated = re.sub(
            r'(?<!\\)\b%s\b' % re.escape(word),
            lambda _match, style=style, word=word: f'\\text{style}{{{word}}}',
            updated,
        )
        changed |= new_updated != updated
        updated = new_updated
    for index, original in enumerate(placeholders):
        updated = updated.replace(f'__PAPERLINT_STYLE_{index}__', original)
    if not changed:
        return line, None
    return updated, (0, len(line))


def fix_missing_word_style_in_file(path):
    path = Path(path)
    old = _read_text(path)
    style_for_word, style_counts, styled_word_pattern = collect_word_style_data(old)
    if not style_for_word:
        return 0
    return _rewrite_lines(path, lambda line: missing_word_style_line(line, style_for_word, style_counts, styled_word_pattern)[0])


def fix_glossary_refs_in_file(path, glossary_file):
    path = Path(path)
    replacements = gs.read_glossary_replacements(glossary_file)
    if not replacements:
        return 0
    return gs.replace_glossary_refs_in_file(str(path), replacements)


def _iter_tex_files(root):
    root = Path(root)
    for path in root.rglob('*.tex'):
        if '/aux/' in str(path) or path.name.endswith('.bak'):
            continue
        yield path


def _apply_prefix(path, prefix, recursive=False, source_prefix=None, write=False):
    path = Path(path)
    source_prefix = source_prefix or prefix
    labels = set()
    tex_files = []
    paths = _iter_tex_files(path) if recursive else [path]
    for tex_path in paths:
        tex_files.append(tex_path)
        text = _read_text(tex_path)
        for label in re.findall(rf'\\label\{{{re.escape(source_prefix)}([^}}]+)\}}', text):
            labels.add(label)

    if not labels:
        return 0

    changed = 0
    pattern = re.compile(r'\\(' + '|'.join(REF_CMDS) + r')\{([^}]*)\}')
    for tex_path in tex_files:
        original = _read_text(tex_path)

        def repl(match):
            command = match.group(1)
            body = match.group(2)
            parts = [part.strip() for part in body.split(',')]
            updated_parts = []
            for part in parts:
                if part.startswith(prefix):
                    updated_parts.append(part)
                elif part in labels:
                    updated_parts.append(f'{prefix}{part}')
                else:
                    updated_parts.append(part)
            return f'\\{command}{{{",".join(updated_parts)}}}'

        updated = pattern.sub(repl, original)
        if updated != original:
            if write:
                _write_text(tex_path, updated)
            changed += 1
    return changed


def preview_prefix(path, prefix, recursive=False, source_prefix=None):
    return _apply_prefix(path, prefix, recursive=recursive, source_prefix=source_prefix, write=False)


def fix_prefix(path, prefix, recursive=False, source_prefix=None):
    return _apply_prefix(path, prefix, recursive=recursive, source_prefix=source_prefix, write=True)
