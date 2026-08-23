import os
import pathlib
import re

import tex_structure as ts

from math_symbol_parser import (
    MATH_SYMBOL_STOPWORDS,
    SYMBOL_DECORATOR_COMMANDS,
    _normalize_scripts_in_expression as shared_normalize_scripts_in_expression,
    _normalize_symbol as shared_normalize_symbol,
    _read_braced_group as shared_read_braced_group,
    _skip_ws as shared_skip_ws,
    _unwrap_braces as shared_unwrap_braces,
    canonicalize_symbol as shared_canonicalize_symbol,
    extract_math_symbols as shared_extract_math_symbols,
    find_math_symbol_matches as shared_find_math_symbol_matches,
    symbol_key as shared_symbol_key,
)

EXPLICIT_MATH_TEXT_COMMANDS = [
    "text",
    "textrm",
    "textit",
    "textbf",
    "textsf",
    "texttt",
    "textsc",
    "mbox",
    "operatorname",
    "gls",
]
GLOSSARY_LABEL_STRIP_COMMANDS = {
    "text",
    "textrm",
    "textit",
    "textbf",
    "textsf",
    "texttt",
    "textsc",
    "mbox",
    "operatorname",
}

GLOSSARY_SORT_STRIP_COMMANDS = SYMBOL_DECORATOR_COMMANDS | GLOSSARY_LABEL_STRIP_COMMANDS
GLOSSARY_IGNORED_SYMBOLS = set(MATH_SYMBOL_STOPWORDS)
GLS_COMMAND_RE = re.compile(r"\\gls[a-zA-Z]*\{([^{}]+)\}")
LABEL_COMMAND_RE = re.compile(r"\\label\*?\{[^{}]*\}")
INCLUDEGRAPHICS_COMMAND_RE = re.compile(r"\\includegraphics\*?(?:\[[^\]]*\])?\{[^{}]*\}")
DOCUMENTCLASS_COMMAND_RE = re.compile(r"\\documentclass\*?(?:\[[^\]]*\])?\{[^{}]*\}")
TITLE_COMMAND_RE = re.compile(r"\\title\*?(?:\[[^\]]*\])?\{[^{}]*\}")
SHORTTITLE_COMMAND_RE = re.compile(r"\\shorttitle\*?(?:\[[^\]]*\])?\{[^{}]*\}")
AUTHOR_COMMAND_RE = re.compile(r"\\author\*?(?:\[[^\]]*\])?\{[^{}]*\}")
AFFILIATION_COMMAND_RE = re.compile(r"\\affiliation\*?(?:\[[^\]]*\])?\{[^{}]*\}")
REF_COMMANDS_RE = re.compile(r"\\(?:ref|eqref|cref|Cref|autoref|pageref|nameref)\*?\{[^{}]*\}")
CITE_COMMANDS_RE = re.compile(r"\\(?:cite|citet|citep)\*?(?:\[[^\]]*\]){0,2}\{[^{}]*\}")
DECLARE_FIGURE_OPTIONS_RE = re.compile(r"\\DeclareFigureOptions\*?\{[^{}]*\}\{[^{}]*\}\{[^{}]*\}")

_skip_ws = shared_skip_ws
_read_braced_group = shared_read_braced_group
_normalize_symbol = shared_normalize_symbol
_normalize_scripts_in_expression = shared_normalize_scripts_in_expression
_unwrap_braces = shared_unwrap_braces
_canonicalize_symbol = shared_canonicalize_symbol
_symbol_key = shared_symbol_key


def _collapse_simple_script_braces(text):
    return re.sub(r"([_^])\{([A-Za-z0-9]+)\}", r"\1\2", text)


def _normalize_glossary_equivalence_symbol(symbol):
    stripped = _normalize_scripts_in_expression(symbol)
    stripped = _strip_wrapped_commands(stripped, GLOSSARY_SORT_STRIP_COMMANDS)
    stripped = _collapse_simple_script_braces(stripped)
    return _normalize_symbol(stripped)


def _prefer_glossary_symbol(existing_symbol, candidate_symbol):
    existing_has_text = "\\text{" in existing_symbol
    candidate_has_text = "\\text{" in candidate_symbol
    if candidate_has_text and not existing_has_text:
        return candidate_symbol
    if existing_has_text and not candidate_has_text:
        return existing_symbol
    if len(candidate_symbol) < len(existing_symbol):
        return candidate_symbol
    return existing_symbol


def read_glossary_replacements(glossary_dir):
    replacements = {}
    if not glossary_dir or not os.path.isdir(glossary_dir):
        return replacements
    for path in sorted(pathlib.Path(glossary_dir).glob("*.tex")):
        for key, symbol, _ in parse_newsym_entries(path):
            label = "sym-%s" % key
            replacements[shared_symbol_key(shared_canonicalize_symbol(symbol))] = label
    return replacements


def parse_acronym_entries(path):
    entries = []
    if not path.exists():
        return entries

    for raw_line in path.read_text().splitlines():
        line = ts.strip_comment_from_line(raw_line).strip()
        if not line.startswith(r"\newabbreviation"):
            continue

        idx = len(r"\newabbreviation")
        idx = _skip_ws(line, idx)

        if idx < len(line) and line[idx] == "[":
            end_idx = idx + 1
            brace_depth = 0
            while end_idx < len(line):
                char = line[end_idx]
                if char == "{":
                    brace_depth += 1
                elif char == "}":
                    if brace_depth > 0:
                        brace_depth -= 1
                elif char == "]" and brace_depth == 0:
                    break
                end_idx += 1
            if end_idx >= len(line) or line[end_idx] != "]":
                continue
            idx = end_idx + 1

        idx = _skip_ws(line, idx)
        label_group, idx = _read_braced_group(line, idx)
        idx = _skip_ws(line, idx)
        short_group, idx = _read_braced_group(line, idx)
        idx = _skip_ws(line, idx)
        long_group, idx = _read_braced_group(line, idx)
        if label_group is None or short_group is None or long_group is None:
            continue

        label = _unwrap_braces(label_group).strip()
        short = _unwrap_braces(short_group).strip()
        long_form = _unwrap_braces(long_group).strip()
        if label and short:
            entries.append((label, short, long_form))
    return entries


def read_acronym_replacements(path):
    replacements = []
    if not path:
        return replacements
    for label, short, long_form in parse_acronym_entries(pathlib.Path(path)):
        if long_form:
            replacements.append((f"{long_form} ({short})", rf"\gls{{{label}}}"))
            if not long_form.endswith("s") and short.endswith("s") is False:
                replacements.append((f"{long_form}s ({short}s)", rf"\glspl{{{label}}}"))
        replacements.append((short, rf"\gls{{{label}}}"))
        replacements.append((f"{short}s", rf"\glspl{{{label}}}"))
        replacements.append((f"{short}'s", rf"\gls{{{label}}}'s"))
        if long_form:
            replacements.append((long_form, rf"\gls{{{label}}}"))
            if not long_form.endswith("s"):
                replacements.append((f"{long_form}s", rf"\glspl{{{label}}}"))

    # Replace longer phrases first so singular does not consume plural forms.
    replacements.sort(key=lambda item: len(item[0]), reverse=True)
    return replacements


def replace_glossary_refs_in_span(span_text, replacements):
    def has_spaced_backslash_prefix(text, idx):
        probe = idx - 1
        saw_space = False
        while probe >= 0 and text[probe].isspace():
            saw_space = True
            probe -= 1
        return saw_space and probe >= 0 and text[probe] == "\\"

    updates = []
    for start, end, symbol in shared_find_math_symbol_matches(span_text):
        if has_spaced_backslash_prefix(span_text, start):
            continue
        label = replacements.get(shared_symbol_key(shared_canonicalize_symbol(symbol)))
        if label is None:
            continue
        updates.append((start, end, label))

    if not updates:
        return span_text, 0

    updated = span_text
    for start, end, label in sorted(updates, key=lambda item: item[0], reverse=True):
        updated = updated[:start] + r"\gls{%s}" % label + updated[end:]
    updated = re.sub(r"([_^])\s*(\\gls[a-zA-Z]*\{[^{}]+\})", r"\1{\2}", updated)
    return updated, len(updates)


def _mask_label_commands(text):
    placeholders = []

    def repl(match):
        placeholders.append(match.group(0))
        return f"__PAPERLINT_LABEL_{len(placeholders) - 1}__"

    masked = LABEL_COMMAND_RE.sub(repl, text)
    masked = INCLUDEGRAPHICS_COMMAND_RE.sub(repl, masked)
    masked = DOCUMENTCLASS_COMMAND_RE.sub(repl, masked)
    masked = TITLE_COMMAND_RE.sub(repl, masked)
    masked = SHORTTITLE_COMMAND_RE.sub(repl, masked)
    masked = AUTHOR_COMMAND_RE.sub(repl, masked)
    masked = AFFILIATION_COMMAND_RE.sub(repl, masked)
    masked = REF_COMMANDS_RE.sub(repl, masked)
    masked = CITE_COMMANDS_RE.sub(repl, masked)
    masked = DECLARE_FIGURE_OPTIONS_RE.sub(repl, masked)
    return masked, placeholders


def _restore_label_commands(text, placeholders):
    restored = text
    for idx, original in enumerate(placeholders):
        restored = restored.replace(f"__PAPERLINT_LABEL_{idx}__", original)
    return restored


def replace_acronym_refs_in_span(span_text, acronym_replacements):
    if not acronym_replacements:
        return span_text, 0

    updates = []
    for source, replacement in acronym_replacements:
        if not source:
            continue
        pattern = re.compile(r"(?<![A-Za-z0-9\\])%s(?![A-Za-z0-9])" % re.escape(source))
        for match in pattern.finditer(span_text):
            updates.append((match.start(), match.end(), replacement))

    if not updates:
        return span_text, 0

    # Keep the longest replacement for overlapping matches.
    updates.sort(key=lambda item: (item[0], -(item[1] - item[0])))
    selected = []
    for start, end, replacement in updates:
        if selected and start < selected[-1][1]:
            continue
        selected.append((start, end, replacement))

    updated = span_text
    for start, end, replacement in reversed(selected):
        updated = updated[:start] + replacement + updated[end:]
    return updated, len(selected)


def replace_glossary_refs_in_line(line, replacements, acronym_replacements=None, equation_line=False):
    content, comment = ts.split_comment_from_line(line)
    if not content:
        return line, 0

    content, label_placeholders = _mask_label_commands(content)

    stripped = content.strip()
    if equation_line and stripped.startswith((r"\begin", r"\end", r"\label", r"\tag", r"\nonumber", r"\notag")):
        return _restore_label_commands(content, label_placeholders) + comment, 0

    total = 0
    if equation_line:
        updated, count = replace_glossary_refs_in_span(content, replacements)
        total += count
        return _restore_label_commands(updated, label_placeholders) + comment, total

    spans = ts.get_math_spans(content)
    if not spans:
        updated, count = replace_acronym_refs_in_span(content, acronym_replacements)
        total += count
        return _restore_label_commands(updated, label_placeholders) + comment, total

    pieces = []
    last_index = 0
    for start, end, _ in spans:
        prose_updated, prose_count = replace_acronym_refs_in_span(content[last_index:start], acronym_replacements)
        pieces.append(prose_updated)
        total += prose_count
        updated, count = replace_glossary_refs_in_span(content[start:end], replacements)
        pieces.append(updated)
        total += count
        last_index = end

    prose_updated, prose_count = replace_acronym_refs_in_span(content[last_index:], acronym_replacements)
    pieces.append(prose_updated)
    total += prose_count
    updated = "".join(pieces)
    return _restore_label_commands(updated, label_placeholders) + comment, total


def replace_glossary_refs_in_file(path, replacements, acronym_replacements=None):
    try:
        with open(path) as handle:
            lines = handle.readlines()
    except Exception:
        return 0

    updated_lines = []
    total = 0
    for i, line in enumerate(lines):
        updated_line, count = replace_glossary_refs_in_line(
            line,
            replacements,
            acronym_replacements=acronym_replacements,
            equation_line=ts.in_equation(i),
        )
        updated_lines.append(updated_line)
        total += count

    if total > 0:
        with open(path, "w") as handle:
            handle.writelines(updated_lines)
    return total


def strip_explicit_math_text(content):
    stripped = content
    pattern = re.compile(r"\\(?:%s)\*?\{[^{}]*\}" % "|".join(EXPLICIT_MATH_TEXT_COMMANDS))
    while True:
        updated = pattern.sub(" ", stripped)
        if updated == stripped:
            return stripped
        stripped = updated


def extract_math_symbols(content):
    return shared_extract_math_symbols(content)


def symbol_to_glossary_label(symbol):
    stripped = symbol[1:] if symbol.startswith("\\") else symbol
    stripped = stripped.replace(r"\prime", "prime")
    stripped = stripped.replace("'", "prime")
    stripped = stripped.replace("*", "star")
    safe = re.sub(r"[^A-Za-z0-9]+", "-", stripped).strip("-")
    if not safe:
        safe = "symbol"
    return "sym-%s" % safe


def _parse_newsym_line(line):
    """Parse a \\newsym{key}{symbol}{description} definition line, if present."""
    stripped = line.strip()
    comment_idx = re.search(r"(?<!\\)%", stripped)
    if comment_idx:
        stripped = stripped[:comment_idx.start()].rstrip()
    if not stripped.startswith(r"\newsym"):
        return None

    idx = len(r"\newsym")
    idx = _skip_ws(stripped, idx)
    key_group, idx = _read_braced_group(stripped, idx)
    if key_group is None:
        return None
    idx = _skip_ws(stripped, idx)
    symbol_group, idx = _read_braced_group(stripped, idx)
    if symbol_group is None:
        return None
    idx = _skip_ws(stripped, idx)
    description_group, idx = _read_braced_group(stripped, idx)
    if description_group is None:
        return None
    if _skip_ws(stripped, idx) != len(stripped):
        return None

    key = _unwrap_braces(key_group).strip()
    symbol = _unwrap_braces(symbol_group)
    description = _unwrap_braces(description_group).strip()
    return key, symbol, description


def _strip_wrapped_commands(text, commands):
    updated = text
    while True:
        changed = False
        for command in commands:
            pattern = re.compile(r"\\%s\s*\{([^{}]*)\}" % re.escape(command))
            replaced = pattern.sub(r"\1", updated)
            if replaced != updated:
                updated = replaced
                changed = True
        if not changed:
            return updated


def _normalize_glossary_label_symbol(symbol):
    stripped = _normalize_scripts_in_expression(symbol)
    stripped = _strip_wrapped_commands(stripped, GLOSSARY_LABEL_STRIP_COMMANDS)
    return _normalize_symbol(stripped)


def _normalize_glossary_sort_symbol(symbol):
    stripped = _normalize_scripts_in_expression(symbol)
    stripped = _strip_wrapped_commands(stripped, GLOSSARY_SORT_STRIP_COMMANDS)
    return _normalize_symbol(stripped)


def parse_newsym_entries(path):
    """Read all \\newsym{key}{symbol}{description} entries from a single file."""
    entries = []
    if not os.path.exists(path):
        return entries
    try:
        with open(path) as handle:
            lines = handle.readlines()
    except Exception:
        return entries
    for line in lines:
        parsed = _parse_newsym_line(line)
        if parsed is None:
            continue
        entries.append(parsed)
    return entries


def read_newsym_symbols_from_dir(glossary_dir, exclude_path=None):
    """Collect the normalized symbols already defined by \\newsym in a glossary folder."""
    symbols = set()
    exclude_resolved = pathlib.Path(exclude_path).resolve() if exclude_path else None
    if not os.path.isdir(glossary_dir):
        return symbols
    for path in sorted(pathlib.Path(glossary_dir).glob("*.tex")):
        if exclude_resolved is not None and path.resolve() == exclude_resolved:
            continue
        for _, symbol, _ in parse_newsym_entries(path):
            symbols.add(_normalize_glossary_equivalence_symbol(_canonicalize_symbol(symbol.strip())))
    return symbols


def collect_glossary_labels_from_files(paths):
    labels = set()
    for path in paths:
        try:
            with open(path) as handle:
                content = ts.strip_tex_comments(handle.read())
        except Exception:
            continue
        for match in GLS_COMMAND_RE.finditer(content):
            label = match.group(1).strip()
            if label:
                labels.add(label)
    return labels


def symbol_to_newsym_key(symbol):
    label = symbol_to_glossary_label(symbol)
    return label[len("sym-"):] if label.startswith("sym-") else label


def build_newsym_lines(symbols, metadata_by_symbol, used_labels=None):
    lines = []
    used_labels = set(used_labels or [])
    filtered_symbols = [symbol for symbol in symbols if _symbol_key(_normalize_glossary_equivalence_symbol(symbol)) not in GLOSSARY_IGNORED_SYMBOLS]
    for symbol in sorted(filtered_symbols, key=lambda item: (
        symbol_to_newsym_key(_normalize_glossary_sort_symbol(item)).lower(),
        symbol_to_newsym_key(_normalize_glossary_label_symbol(item)),
        item,
    )):
        key_base = symbol_to_newsym_key(_normalize_glossary_label_symbol(symbol))
        lookup_key = _symbol_key(_normalize_glossary_equivalence_symbol(symbol))
        key = key_base
        suffix = 2
        while ("sym-%s" % key) in used_labels:
            key = "%s-%d" % (key_base, suffix)
            suffix += 1
        used_labels.add("sym-%s" % key)
        description = metadata_by_symbol.get(lookup_key, {}).get("description", "")
        lines.append(f"\\newsym{{{key}}}{{{symbol}}}{{{description}}}")
    return lines


def write_symbol_glossary(path, symbols, glossary_dir=None, used_labels=None):
    """Write newly discovered symbols not already documented in glossary_dir to path.

    Entries already present in path from a previous run are carried over (with
    their description) unless a symbol has since been documented elsewhere in
    glossary_dir, in which case it is dropped from the generated file.
    """
    known_symbols = read_newsym_symbols_from_dir(glossary_dir, exclude_path=path) if glossary_dir else set()

    grouped_symbols = {}
    metadata_by_symbol = {}
    for key, symbol, description in parse_newsym_entries(path):
        normalized_key = _symbol_key(_normalize_glossary_equivalence_symbol(symbol))
        if normalized_key in known_symbols:
            continue
        grouped_symbols[normalized_key] = symbol
        metadata_by_symbol[normalized_key] = {"description": description}

    for symbol in set(symbols):
        normalized_key = _symbol_key(_normalize_glossary_equivalence_symbol(symbol))
        if normalized_key in known_symbols:
            continue
        existing_symbol = grouped_symbols.get(normalized_key)
        if existing_symbol is None:
            grouped_symbols[normalized_key] = symbol
        else:
            grouped_symbols[normalized_key] = _prefer_glossary_symbol(existing_symbol, symbol)

    lines = build_newsym_lines(set(grouped_symbols.values()), metadata_by_symbol, used_labels=used_labels)

    with open(path, "w") as handle:
        for line in lines:
            handle.write(line + "\n")

