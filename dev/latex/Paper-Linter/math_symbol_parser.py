#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path

ENTRY_RE = re.compile(
    r"\\glsxtrnewsymbol\[(?P<attrs>[^\]]*)\]\{(?P<label>[^}]+)\}\{\\ensuremath\{(?P<body>.*)\}\}"
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
GREEK_MATH_SYMBOLS = {
    "alpha", "beta", "gamma", "delta", "epsilon", "varepsilon", "zeta", "eta", "theta", "vartheta",
    "iota", "kappa", "lambda", "mu", "nu", "xi", "pi", "varpi", "rho", "varrho", "sigma", "varsigma",
    "tau", "upsilon", "phi", "varphi", "chi", "psi", "omega", "Gamma", "Delta", "Theta", "Lambda",
    "Xi", "Pi", "Sigma", "Upsilon", "Phi", "Psi", "Omega",
}
MATH_SYMBOL_STOPWORDS = {
    "sin", "cos", "tan", "cot", "sec", "csc", "log", "ln", "exp", "min", "max", "arg", "det",
    "dim", "gcd", "hom", "ker", "rank", "span", "supp", "mod", "Pr", "Re", "Im", "sgn", "diag",
    "lim", "limsup", "liminf", "sup", "inf", "left", "right", "big", "Big", "bigg", "Bigg", "e", "i", "d",
}
SYMBOL_DECORATOR_COMMANDS = {
    "vec", "mathbf", "boldsymbol", "bm", "hat", "dot", "ddot", "bar", "tilde",
    "mathit", "mathsf", "mathtt", "mathcal", "mathbb", "mathfrak", "mathrm",
}
SPECIAL_MATHBB_BASES = {r"\mathbb{R}", r"\mathbb{C}"}


def _skip_ws(text: str, idx: int) -> int:
    while idx < len(text) and text[idx].isspace():
        idx += 1
    return idx


def _read_command(text: str, idx: int):
    if idx >= len(text) or text[idx] != "\\":
        return None, idx
    end = idx + 1
    while end < len(text) and (text[end].isalpha() or text[end] == "@"):
        end += 1
    return text[idx:end], end


def _read_braced_group(text: str, idx: int):
    if idx >= len(text) or text[idx] != "{":
        return None, idx
    depth = 0
    end = idx
    while end < len(text):
        if text[end] == "{":
            depth += 1
        elif text[end] == "}":
            depth -= 1
            if depth == 0:
                return text[idx:end + 1], end + 1
        end += 1
    return None, idx


def _read_script_value(text: str, idx: int):
    idx = _skip_ws(text, idx)
    if idx >= len(text):
        return None, idx
    if text[idx] == "{":
        return _read_braced_group(text, idx)
    if text[idx] == "\\":
        command, command_end = _read_command(text, idx)
        if command is None:
            return None, idx
        next_idx = _skip_ws(text, command_end)
        if next_idx < len(text) and text[next_idx] == "{":
            group, group_end = _read_braced_group(text, next_idx)
            if group is not None:
                return command + group, group_end
        return command, command_end
    if text[idx].isalnum():
        # TeX takes a single token for unbraced scripts, not an alnum run.
        return text[idx], idx + 1
    return text[idx], idx + 1


def _normalize_symbol(symbol: str) -> str:
    return re.sub(r"\s+", " ", symbol).strip()


def symbol_key(symbol: str) -> str:
    return re.sub(r"\s+", "", symbol)


def _is_simple_script_token(value: str) -> bool:
    return re.fullmatch(r"[A-Za-z0-9]+", value) is not None or re.fullmatch(r"\\[A-Za-z@]+(?:\{[^{}]*\})?", value) is not None


def _normalize_script_value(script_value: str):
    normalized_value = _normalize_symbol(script_value)
    if not normalized_value:
        return None
    if normalized_value.startswith("{") and normalized_value.endswith("}"):
        inner_value = _normalize_script_value(normalized_value[1:-1])
        if inner_value is None:
            return None
        if _is_simple_script_token(inner_value):
            return inner_value
        return "{%s}" % inner_value
    normalized_value = re.sub(r",+\Z", "", normalized_value).strip()
    if not normalized_value:
        return None
    return normalized_value


def _unwrap_braces(value: str) -> str:
    if len(value) >= 2 and value[0] == "{" and value[-1] == "}":
        return value[1:-1]
    return value


def _script_is_ignored_superscript(script_value: str) -> bool:
    raw = _unwrap_braces(script_value)
    if raw in {"T", "*"}:
        return True
    if re.fullmatch(r"[+-]?\d+", raw):
        return True
    if re.fullmatch(r"[+-]?\d+\\times[+-]?\d+", raw):
        return True
    return False


def _canonicalize_script(script_op: str, script_value: str):
    original_value = _normalize_symbol(script_value)
    had_outer_braces = original_value.startswith("{") and original_value.endswith("}")
    original_inner_value = _normalize_symbol(_unwrap_braces(original_value)) if had_outer_braces else original_value
    normalized_value = _normalize_script_value(script_value)
    if normalized_value is None:
        return None
    if script_op == "^" and _script_is_ignored_superscript(normalized_value):
        return None
    # Preserve explicitly braced control-sequence-like scripts with spacing,
    # e.g., ^{\ e}, so they do not collapse into undefined commands like ^\e.
    if had_outer_braces and re.search(r"\\\s+[A-Za-z@]", original_inner_value):
        return script_op + "{" + original_inner_value + "}"
    # Preserve explicit braces around control sequences when they were present
    # in the source script.
    if had_outer_braces and original_inner_value.startswith("\\"):
        normalized_inner = _unwrap_braces(normalized_value)
        return script_op + "{" + normalized_inner + "}"
    # Preserve braces when they were explicitly present around a multi-character
    # simple script value (e.g., a_{in}, x^{ij}).
    if had_outer_braces and re.fullmatch(r"[A-Za-z0-9]{2,}", normalized_value):
        normalized_value = "{%s}" % normalized_value
    return script_op + normalized_value


def _normalize_scripts_in_expression(expr: str) -> str:
    canonical = ""
    idx = 0
    while idx < len(expr):
        if expr[idx] not in "_^":
            canonical += expr[idx]
            idx += 1
            continue
        script_op = expr[idx]
        script_value, next_idx = _read_script_value(expr, idx + 1)
        if script_value is None:
            canonical += expr[idx:]
            break
        updated = _canonicalize_script(script_op, script_value)
        if updated is not None:
            canonical += updated
        idx = next_idx
    return canonical


def canonicalize_symbol(symbol: str) -> str:
    return _normalize_scripts_in_expression(symbol)


def _group_contains_gls_command(group: str) -> bool:
    return re.search(r"\\gls\s*\{", group) is not None


def _decorator_group_has_symbol(group: str) -> bool:
    inner = _unwrap_braces(group)
    return len(extract_math_symbols(inner)) > 0


def strip_explicit_math_text(content: str) -> str:
    stripped = content
    pattern = re.compile(r"\\(?:%s)\*?\{[^{}]*\}" % "|".join(EXPLICIT_MATH_TEXT_COMMANDS))
    while True:
        updated = pattern.sub(lambda match: " " * (match.end() - match.start()), stripped)
        if updated == stripped:
            return stripped
        stripped = updated


def _collect_math_symbol_matches(content: str):
    cleaned = content
    cleaned = re.sub(r"\\(?:begin|end)\{[^{}]+\}", lambda match: " " * (match.end() - match.start()), cleaned)
    cleaned = re.sub(
        r"\\(?:label|tag|nonumber|notag)\*?(?:\{[^{}]*\})?",
        lambda match: " " * (match.end() - match.start()),
        cleaned,
    )
    cleaned = re.sub(
        r"\\(?:left|right)(?![A-Za-z@])|\\(?:,|;|:|!|quad|qquad|medspace|thinspace|enspace)",
        lambda match: " " * (match.end() - match.start()),
        cleaned,
    )
    matches = []
    idx = 0
    while idx < len(cleaned):
        base = None
        next_idx = idx + 1
        char = cleaned[idx]

        if char.isalpha():
            prev_alpha = idx > 0 and cleaned[idx - 1].isalpha()
            next_alpha = idx + 1 < len(cleaned) and cleaned[idx + 1].isalpha()
            if prev_alpha or next_alpha:
                idx += 1
                continue
            base = char
            next_idx = idx + 1
        elif char == "\\":
            command, command_end = _read_command(cleaned, idx)
            if command is None:
                idx += 1
                continue
            name = command[1:]
            if name in GREEK_MATH_SYMBOLS:
                base = command
                next_idx = command_end
            elif name == "gls":
                skip_idx = _skip_ws(cleaned, command_end)
                if skip_idx < len(cleaned) and cleaned[skip_idx] == "{":
                    _, skip_end = _read_braced_group(cleaned, skip_idx)
                    if skip_end > skip_idx:
                        idx = skip_end
                    else:
                        idx = command_end
                else:
                    idx = command_end
                continue
            elif name in SYMBOL_DECORATOR_COMMANDS:
                arg_start = _skip_ws(cleaned, command_end)
                if arg_start < len(cleaned) and cleaned[arg_start] == "{":
                    group, group_end = _read_braced_group(cleaned, arg_start)
                    if (
                        group is not None
                        and _normalize_symbol(_unwrap_braces(group))
                        and not _group_contains_gls_command(group)
                        and _decorator_group_has_symbol(group)
                    ):
                        base = command + group
                        next_idx = group_end
                if base is None:
                    idx = command_end
                    continue
            else:
                idx = command_end
                continue
        else:
            idx += 1
            continue

        symbol = base
        base_is_special_mathbb = _normalize_symbol(base) in SPECIAL_MATHBB_BASES
        script_idx = next_idx
        replace_end = next_idx
        while True:
            script_idx = _skip_ws(cleaned, script_idx)
            if script_idx >= len(cleaned) or cleaned[script_idx] not in "_^":
                break
            script_op = cleaned[script_idx]
            script_val, script_end = _read_script_value(cleaned, script_idx + 1)
            if script_val is None:
                break
            if base == "e" and script_op == "^":
                idx = next_idx
                break
            if base_is_special_mathbb and script_op == "^":
                # Keep the superscript out of the base symbol and let its
                # contents be parsed as independent symbols.
                break
            canonical_script = _canonicalize_script(script_op, script_val)
            if canonical_script is not None:
                symbol += canonical_script
                replace_end = script_end
            script_idx = script_end
        if base == "e" and idx == next_idx:
            continue

        normalized = _normalize_symbol(canonicalize_symbol(symbol.strip()))
        normalized_key = symbol_key(normalized)
        if normalized_key.startswith("\\"):
            base_name = normalized_key[1:].split("{", 1)[0]
            if base_name in MATH_SYMBOL_STOPWORDS:
                idx = script_idx
                continue
        elif normalized_key in MATH_SYMBOL_STOPWORDS:
            idx = script_idx
            continue

        matches.append((idx, replace_end, normalized))
        idx = script_idx
    return matches


def extract_math_symbols(content: str):
    return {item[2] for item in _collect_math_symbol_matches(content)}


def find_math_symbol_matches(content: str):
    return _collect_math_symbol_matches(content)


def parse_glossary_entries(path: Path):
    entries = []
    for line in path.read_text().splitlines():
        match = ENTRY_RE.fullmatch(line.strip())
        if not match:
            continue
        symbol_set = extract_math_symbols(match.group("body"))
        if not symbol_set:
            continue
        symbol = sorted(symbol_set, key=len, reverse=True)[0]
        entries.append((match.group("label"), symbol, match.group("body")))
    return entries


def load_glossary_symbol_map(path: Path):
    symbol_map = {}
    if not path.exists():
        return symbol_map
    for label, symbol, _ in parse_glossary_entries(path):
        symbol_map[symbol_key(canonicalize_symbol(symbol))] = label
    return symbol_map