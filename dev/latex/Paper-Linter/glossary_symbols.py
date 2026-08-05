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
    parse_glossary_entries,
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

_skip_ws = shared_skip_ws
_read_braced_group = shared_read_braced_group
_normalize_symbol = shared_normalize_symbol
_normalize_scripts_in_expression = shared_normalize_scripts_in_expression
_unwrap_braces = shared_unwrap_braces
_canonicalize_symbol = shared_canonicalize_symbol
_symbol_key = shared_symbol_key


def read_glossary_replacements(path):
    replacements = {}
    if not os.path.exists(path):
        return replacements
    for label, symbol, _ in parse_glossary_entries(pathlib.Path(path)):
        replacements[shared_symbol_key(shared_canonicalize_symbol(symbol))] = label
    return replacements


def replace_glossary_refs_in_span(span_text, replacements):
    updates = []
    for start, end, symbol in shared_find_math_symbol_matches(span_text):
        label = replacements.get(shared_symbol_key(shared_canonicalize_symbol(symbol)))
        if label is None:
            continue
        updates.append((start, end, label))

    if not updates:
        return span_text, 0

    updated = span_text
    for start, end, label in sorted(updates, key=lambda item: item[0], reverse=True):
        updated = updated[:start] + r"\gls{%s}" % label + updated[end:]
    return updated, len(updates)


def replace_glossary_refs_in_line(line, replacements, equation_line=False):
    content, comment = ts.split_comment_from_line(line)
    if not content:
        return line, 0

    stripped = content.strip()
    if equation_line and stripped.startswith((r"\begin", r"\end", r"\label", r"\tag", r"\nonumber", r"\notag")):
        return line, 0

    total = 0
    if equation_line:
        updated, count = replace_glossary_refs_in_span(content, replacements)
        total += count
        return updated + comment, total

    spans = ts.get_math_spans(content)
    if not spans:
        return line, 0

    pieces = []
    last_index = 0
    for start, end, _ in spans:
        pieces.append(content[last_index:start])
        updated, count = replace_glossary_refs_in_span(content[start:end], replacements)
        pieces.append(updated)
        total += count
        last_index = end

    pieces.append(content[last_index:])
    return "".join(pieces) + comment, total


def replace_glossary_refs_in_file(path, replacements):
    try:
        with open(path) as handle:
            lines = handle.readlines()
    except Exception:
        return 0

    updated_lines = []
    total = 0
    for i, line in enumerate(lines):
        updated_line, count = replace_glossary_refs_in_line(line, replacements, equation_line=ts.in_equation(i))
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
    safe = re.sub(r"[^A-Za-z0-9]+", "-", stripped).strip("-")
    if not safe:
        safe = "symbol"
    return "sym-%s" % safe


def _split_top_level_options(options):
    parts = []
    current = ""
    depth = 0
    for ch in options:
        if ch == "{":
            depth += 1
        elif ch == "}":
            if depth > 0:
                depth -= 1
        if ch == "," and depth == 0:
            if current.strip():
                parts.append(current.strip())
            current = ""
            continue
        current += ch
    if current.strip():
        parts.append(current.strip())
    return parts


def _extract_option_value(options, wanted_key):
    for option in _split_top_level_options(options):
        if "=" not in option:
            continue
        key, value = option.split("=", 1)
        if key.strip() != wanted_key:
            continue
        value = value.strip()
        group, end_idx = _read_braced_group(value, 0)
        if group is not None and end_idx == len(value):
            return _unwrap_braces(group)
        return value
    return None


def _extract_description_option(options):
    return _extract_option_value(options, "description")


def _extract_sort_option(options):
    return _extract_option_value(options, "sort")


def _parse_glsxtrnewsymbol_line(line):
    stripped = line.strip()
    comment_idx = re.search(r"(?<!\\)%", stripped)
    if comment_idx:
        stripped = stripped[:comment_idx.start()].rstrip()
    if not stripped.startswith(r"\glsxtrnewsymbol"):
        return None, None

    idx = len(r"\glsxtrnewsymbol")
    idx = _skip_ws(stripped, idx)
    options = ""
    if idx < len(stripped) and stripped[idx] == "[":
        end_idx = idx + 1
        brace_depth = 0
        while end_idx < len(stripped):
            char = stripped[end_idx]
            if char == "{":
                brace_depth += 1
            elif char == "}":
                if brace_depth > 0:
                    brace_depth -= 1
            elif char == "]" and brace_depth == 0:
                break
            end_idx += 1
        if end_idx >= len(stripped) or stripped[end_idx] != "]":
            return None, None
        options = stripped[idx + 1:end_idx]
        idx = end_idx + 1

    idx = _skip_ws(stripped, idx)
    label_group, idx = _read_braced_group(stripped, idx)
    if label_group is None:
        return None, None
    idx = _skip_ws(stripped, idx)
    ensuremath_group, idx = _read_braced_group(stripped, idx)
    if ensuremath_group is None:
        return None, None
    if _skip_ws(stripped, idx) != len(stripped):
        return None, None

    ensuremath_content = _unwrap_braces(ensuremath_group).strip()
    if not ensuremath_content.startswith(r"\ensuremath"):
        return None, None
    symbol_idx = _skip_ws(ensuremath_content, len(r"\ensuremath"))
    symbol_group, symbol_end = _read_braced_group(ensuremath_content, symbol_idx)
    if symbol_group is None:
        return None, None
    if _skip_ws(ensuremath_content, symbol_end) != len(ensuremath_content):
        return None, None

    symbol = _unwrap_braces(symbol_group)
    description = _extract_description_option(options)
    sort = _extract_sort_option(options)
    metadata = {}
    if description is not None:
        metadata["description"] = description
    if sort is not None:
        metadata["sort"] = sort
    return symbol, metadata


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


def read_existing_symbol_metadata(path):
    metadata_by_symbol = {}
    if not os.path.exists(path):
        return metadata_by_symbol

    try:
        with open(path) as handle:
            lines = handle.readlines()
    except Exception:
        return metadata_by_symbol

    for line in lines:
        symbol, metadata = _parse_glsxtrnewsymbol_line(line)
        if symbol is None or not metadata:
            continue
        symbol = _canonicalize_symbol(symbol.strip())
        existing = metadata_by_symbol.setdefault(_symbol_key(symbol), {})
        for key, value in metadata.items():
            if value or key not in existing:
                existing[key] = value
    return metadata_by_symbol


def read_existing_symbols(path):
    symbols = set()
    if not os.path.exists(path):
        return symbols
    try:
        with open(path) as handle:
            lines = handle.readlines()
    except Exception:
        return symbols
    for line in lines:
        symbol, _ = _parse_glsxtrnewsymbol_line(line)
        if symbol is None:
            continue
        symbols.add(_canonicalize_symbol(symbol.strip()))
    return symbols


def build_symbol_glossary_lines(symbols, metadata_by_symbol):
    lines = []
    used_labels = set()
    filtered_symbols = [symbol for symbol in symbols if _symbol_key(symbol) not in GLOSSARY_IGNORED_SYMBOLS]
    for symbol in sorted(filtered_symbols, key=lambda item: (
        metadata_by_symbol.get(_symbol_key(item), {}).get("sort", symbol_to_glossary_label(_normalize_glossary_sort_symbol(item))).lower(),
        symbol_to_glossary_label(_normalize_glossary_label_symbol(item)),
        item,
    )):
        label_base = symbol_to_glossary_label(_normalize_glossary_label_symbol(symbol))
        sort_value = metadata_by_symbol.get(_symbol_key(symbol), {}).get("sort")
        if sort_value is None:
            sort_value = symbol_to_glossary_label(_normalize_glossary_sort_symbol(symbol))
        label = label_base
        suffix = 2
        while label in used_labels:
            label = "%s-%d" % (label_base, suffix)
            suffix += 1
        used_labels.add(label)
        metadata = metadata_by_symbol.get(_symbol_key(symbol), {})
        description = metadata.get("description", "")
        sort_option = ""
        if metadata.get("sort") or sort_value != label_base:
            sort_option = f",sort={sort_value}"
        lines.append(f"\\glsxtrnewsymbol[description={{{description}}}{sort_option}]{{{label}}}{{\\ensuremath{{{symbol}}}}}")
    return lines


def write_symbol_glossary(path, symbols, seed_path=None):
    source_path = seed_path if seed_path else path
    symbols = set(symbols) | read_existing_symbols(source_path)
    metadata_by_symbol = read_existing_symbol_metadata(source_path)
    lines = build_symbol_glossary_lines(symbols, metadata_by_symbol)
    with open(path, "w") as handle:
        for line in lines:
            handle.write(line + "\n")
