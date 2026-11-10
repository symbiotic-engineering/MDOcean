#!/usr/bin/env python3

import re
import sys


GLS_RE = re.compile(r"\\gls\s*\{([^{}]+)\}")
GLSPL_RE = re.compile(r"\\glspl\s*\{([^{}]+)\}")

NEWSYM_START_RE = re.compile(r"\\newsym\s*\{")
NEWABBREVIATION_START_RE = re.compile(
    r"\\newabbreviation(?:\s*\[[^\]]*\])?\s*\{"
)

DIF_INLINE_RE = re.compile(
    r"\\DIF(?:add|del)(?:FL)?\{([^{}]*)\}"
)


# ---------------------------------------------------------------------------
# Balanced-brace parsing
# ---------------------------------------------------------------------------

def read_braced_argument(text, start):
    """
    Read a LaTeX {...} argument starting at `start`.

    Returns:
        (contents, position_after_closing_brace)

    Handles arbitrary nested braces, e.g.

        {\\ddot{\\xi}}

    and:

        {\\gls{acr-pto} force}
    """

    if start >= len(text) or text[start] != "{":
        raise ValueError(f"expected '{{' at position {start}")

    depth = 0
    i = start

    while i < len(text):
        char = text[i]

        # Ignore escaped characters.
        if char == "\\":
            i += 2
            continue

        if char == "{":
            depth += 1

        elif char == "}":
            depth -= 1

            if depth == 0:
                return text[start + 1:i], i + 1

        i += 1

    raise ValueError("unmatched '{'")


def skip_whitespace(text, pos):
    while pos < len(text) and text[pos].isspace():
        pos += 1
    return pos
def comparison_normalize(text, glossary_map):
    """
    Normalize a diff hunk for equivalence testing.

    Nested DIFadd/DIFdel wrappers are removed, glossary references are
    expanded, and whitespace is ignored.

    The original source is never modified by this function.
    """

    # Remove simple nested \DIFadd{...} / \DIFdel{...} wrappers.
    while True:
        new_text, count = DIF_INLINE_RE.subn(r"\1", text)

        if count == 0:
            break

        text = new_text

    text = normalize(text, glossary_map)

    # Ignore whitespace for equivalence testing.
    return re.sub(r"\s+", "", text)
# ---------------------------------------------------------------------------
# Extract glossary definitions
# ---------------------------------------------------------------------------

def extract_glossary_map(text):
    """
    Extract:

        \\newsym{D-d}{D_d}{Description}

    as:

        sym-D-d -> D_d

    and:

        \\newabbreviation[type=econ]{acr-lcoe}{LCOE}{Levelized Cost of Energy}

    as:

        acr-lcoe -> LCOE
    """

    mappings = {}

    # -----------------------------------------------------------------------
    # \newsym
    # -----------------------------------------------------------------------

    for match in NEWSYM_START_RE.finditer(text):
        pos = match.end() - 1

        try:
            key, pos = read_braced_argument(text, pos)

            pos = skip_whitespace(text, pos)
            symbol, pos = read_braced_argument(text, pos)

            pos = skip_whitespace(text, pos)
            _, pos = read_braced_argument(text, pos)

        except ValueError as exc:
            print(
                f"WARNING: could not parse \\newsym near position "
                f"{match.start()}: {exc}",
                file=sys.stderr,
            )
            continue

        mappings[f"sym-{key}"] = symbol

    # -----------------------------------------------------------------------
    # \newabbreviation
    # -----------------------------------------------------------------------

    for match in NEWABBREVIATION_START_RE.finditer(text):
        pos = match.end() - 1

        try:
            key, pos = read_braced_argument(text, pos)

            pos = skip_whitespace(text, pos)
            abbreviation, pos = read_braced_argument(text, pos)

            pos = skip_whitespace(text, pos)
            _, pos = read_braced_argument(text, pos)

        except ValueError as exc:
            print(
                f"WARNING: could not parse \\newabbreviation near position "
                f"{match.start()}: {exc}",
                file=sys.stderr,
            )
            continue

        mappings[key] = abbreviation

    return mappings


# ---------------------------------------------------------------------------
# Normalize glossary references
# ---------------------------------------------------------------------------
def normalize(text, glossary_map):
    """
    Replace \\gls{...} and \\glspl{...} with their literal forms.

    Used only for equivalence comparison.
    """

    def replace_glspl(match):
        key = match.group(1)

        if key in glossary_map:
            return glossary_map[key] + "s"

        print(
            f"WARNING: unknown glossary key: {key}",
            file=sys.stderr,
        )

        return match.group(0)

    def replace_gls(match):
        key = match.group(1)

        if key in glossary_map:
            return glossary_map[key]

        print(
            f"WARNING: unknown glossary key: {key}",
            file=sys.stderr,
        )

        return match.group(0)

    text = GLSPL_RE.sub(replace_glspl, text)
    text = GLS_RE.sub(replace_gls, text)

    return text


# ---------------------------------------------------------------------------
# Parse one latexdiff deletion/addition pair
# ---------------------------------------------------------------------------
def parse_diff_block(text, pos, kind, fl):
    """
    Parse one latexdiff block.

    Supports both:

        \\DIFaddbegin \\DIFadd{...} \\DIFaddend

    and:

        \\DIFaddbegin ... \\DIFaddend

    The content is parsed with balanced braces, so nested \\gls{...}
    commands are handled correctly.

    Returns:
        (content, end_position)

    or None if the block cannot be parsed.
    """

    if fl:
        begin = rf"\DIF{kind}beginFL"
        command = rf"\DIF{kind}FL"
        end = rf"\DIF{kind}endFL"
    else:
        begin = rf"\DIF{kind}begin"
        command = rf"\DIF{kind}"
        end = rf"\DIF{kind}end"

    if not text.startswith(begin, pos):
        return None

    pos += len(begin)
    pos = skip_whitespace(text, pos)

    # ---------------------------------------------------------------
    # Explicit \DIFadd{...} / \DIFdel{...}
    # ---------------------------------------------------------------

    if text.startswith(command, pos):
        command_end = pos + len(command)

        # Make sure this really is the command and not something like
        # \DIFaddend.
        if command_end < len(text) and text[command_end] == "{":
            content, pos = read_braced_argument(text, command_end)

            pos = skip_whitespace(text, pos)

            if not text.startswith(end, pos):
                return None

            pos += len(end)

            return content, pos

    # ---------------------------------------------------------------
    # Unwrapped form:
    #
    #   \DIFaddbegin WEC \DIFaddend
    #
    # Find the matching end marker. Since the content itself can contain
    # braces, this does not use a naive closing-brace search.
    # ---------------------------------------------------------------

    end_match = re.search(re.escape(end), text[pos:])

    if end_match is None:
        return None

    end_pos = pos + end_match.start()

    content = text[pos:end_pos]

    return content.rstrip(), end_pos + len(end)


def parse_diff_pair(text, start):
    """
    Parse a complete adjacent latexdiff deletion/addition pair.

    Both braced and unbraced DIF blocks are supported.
    """

    pos = start

    # ---------------------------------------------------------------
    # Determine FL vs non-FL
    # ---------------------------------------------------------------

    if text.startswith(r"\DIFdelbeginFL", pos):
        fl = True
    elif text.startswith(r"\DIFdelbegin", pos):
        fl = False
    else:
        return None

    # ---------------------------------------------------------------
    # Parse deletion
    # ---------------------------------------------------------------

    result = parse_diff_block(text, pos, "del", fl)

    if result is None:
        return None

    old, pos = result

    pos = skip_whitespace(text, pos)

    # ---------------------------------------------------------------
    # Parse addition
    # ---------------------------------------------------------------

    result = parse_diff_block(text, pos, "add", fl)

    if result is None:
        return None

    new, pos = result

    return old, new, pos
# ---------------------------------------------------------------------------
# Remove equivalent diffs
# ---------------------------------------------------------------------------

def remove_equivalent_diffs(text, glossary_map):
    """
    Find adjacent latexdiff deletion/addition pairs.

    If their contents become identical after replacing \\gls{...} with
    the corresponding literal symbol/abbreviation, remove the diff markup
    and keep the normalized new text.
    """

    output = []
    pos = 0

    while pos < len(text):

        # Look for the next possible deletion block.
        match = re.search(r"\\DIFdelbegin(?:FL)?", text[pos:])

        if match is None:
            output.append(text[pos:])
            break

        start = pos + match.start()

        # Copy everything before the candidate diff unchanged.
        output.append(text[pos:start])

        try:
            parsed = parse_diff_pair(text, start)
        except ValueError as exc:
            print(
                f"WARNING: could not parse diff near position "
                f"{start}: {exc}",
                file=sys.stderr,
            )
            output.append(text[start])
            pos = start + 1
            continue

        if parsed is None:
            # This wasn't a complete deletion/addition pair.
            output.append(text[start])
            pos = start + 1
            continue

        old, new, end_pos = parsed

        old_normalized = comparison_normalize(old, glossary_map)
        new_normalized = comparison_normalize(new, glossary_map)

        if old_normalized == new_normalized:
            # The only difference was glossary markup.
            #
            # Return the normalized NEW source, so the output contains
            # the literal symbol/abbreviation rather than \gls{...}.
            output.append(new)
        else:
            # Keep the original diff unchanged.
            output.append(text[start:end_pos])

        pos = end_pos

    return "".join(output)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    text = sys.stdin.read()

    glossary_map = extract_glossary_map(text)

    text = remove_equivalent_diffs(text, glossary_map)

    sys.stdout.write(text)

if __name__ == "__main__":
    main()