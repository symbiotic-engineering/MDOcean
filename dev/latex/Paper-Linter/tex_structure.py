import os
import re

tex = None
tex_lines = None
tex_lines_clean = None
tex_lines_math_masked = None
tex_lines_math_texttt_masked = None
in_env = None
envs = None
is_document_root = False
current_file_path = None

FLOAT_ENVS = ["figure", "listing", "table"]
SECTION_COMMAND_RE = re.compile("\\\\(section|subsection|subsubsection|paragraph)\\*?\\{([^\\}]+)\\}")
SENTENCE_END_RE = re.compile(r"(?<!\\)[.!?](?=(?:['\")\]}]*)(?:\s|$))")
COMMON_SENTENCE_ABBREVIATIONS = ("lin", "nonlin", "ineq", "homog", "partic")


def file_has_document_environment(path):
    try:
        with open(path) as handle:
            content = handle.read()
        return re.search(r"\\begin\{document\}", content) is not None
    except Exception:
        return False


def strip_tex_comments(text):
    lines = []
    for line in text.split("\n"):
        match = re.search(r"(?<!\\)%", line)
        if match:
            lines.append(line[:match.start()])
        else:
            lines.append(line)
    return "\n".join(lines)


def resolve_included_file(base_file, target):
    target = target.strip().strip('"').strip("'")
    if not target:
        return None

    base_dir = os.path.dirname(base_file)
    candidates = []

    if os.path.isabs(target):
        candidates.append(target)
    else:
        candidates.append(os.path.join(base_dir, target))

    if not os.path.splitext(target)[1]:
        extensions = [".tex", ".tikz", ".ltx"]
        for candidate in list(candidates):
            for extension in extensions:
                candidates.append(candidate + extension)

    for candidate in candidates:
        if os.path.isfile(candidate):
            return os.path.normpath(candidate)
    return None


def should_ignore_file(file_path, ignore_specs):
    normalized_path = os.path.normpath(file_path)
    file_name = os.path.basename(normalized_path)
    for spec in ignore_specs:
        normalized_spec = os.path.normpath(spec)
        if normalized_path == normalized_spec:
            return True
        if file_name == os.path.basename(normalized_spec):
            return True
        if normalized_path.endswith(normalized_spec):
            return True
        if normalized_path.endswith(spec):
            return True
    return False


def collect_tex_files(root_path, ignore_specs=None):
    if ignore_specs is None:
        ignore_specs = []

    if root_path.endswith(".tex"):
        initial_files = [os.path.normpath(root_path)]
    else:
        all_tex_files = []
        root_tex_files = []
        for path, subdirs, files in os.walk(root_path):
            for f in files:
                if f.endswith(".tex"):
                    file_path = os.path.join(path, f)
                    file_path = os.path.normpath(file_path)
                    all_tex_files.append(file_path)
                    if file_has_document_environment(file_path):
                        root_tex_files.append(file_path)
        initial_files = root_tex_files if len(root_tex_files) > 0 else all_tex_files

    initial_files = [file_path for file_path in initial_files if not should_ignore_file(file_path, ignore_specs)]

    included_pattern = re.compile(r"\\(?:input|include)\*?\{([^\}]+)\}")
    collected = []
    seen = set()

    def visit(file_path):
        file_path = os.path.normpath(file_path)
        if should_ignore_file(file_path, ignore_specs):
            return
        if file_path in seen:
            return
        seen.add(file_path)
        collected.append(file_path)

        try:
            with open(file_path) as handle:
                content = handle.read()
        except Exception:
            return

        content = strip_tex_comments(content)
        for match in included_pattern.finditer(content):
            included_file = resolve_included_file(file_path, match.group(1))
            if included_file is not None:
                visit(included_file)

    for file_path in initial_files:
        visit(file_path)

    return collected


def strip_comment_from_line(line):
    match = re.search(r"(?<!\\)%", line)
    if match:
        return line[:match.start()]
    return line


def split_comment_from_line(line):
    match = re.search(r"(?<!\\)%", line)
    if match:
        return line[:match.start()], line[match.start():]
    return line, ""


def mask_false_sentence_ends(text):
    masked = text
    for abbreviation in COMMON_SENTENCE_ABBREVIATIONS:
        masked = re.sub(r"\b%s\.(?=\s|$)" % re.escape(abbreviation), abbreviation + "§", masked, flags=re.IGNORECASE)
    masked = re.sub(r"(?<!\w)[A-Za-z]\.(?=\s|$)", lambda match: match.group(0)[:-1] + "§", masked)
    return masked


def next_file(file):
    global tex, tex_lines, tex_lines_clean, tex_lines_math_masked, tex_lines_math_texttt_masked, in_env, envs, is_document_root, current_file_path
    try:
        tex = open(file).read()
    except Exception:
        raise RuntimeError("Could not open '%s'" % file)

    current_file_path = file
    is_document_root = re.search(r"\\begin\{document\}", tex) is not None

    tex_lines = tex.split("\n")
    tex_lines_clean = tex.split("\n")
    tex_lines_math_masked = tex.split("\n")
    tex_lines_math_texttt_masked = tex.split("\n")
    in_env = {}
    envs = {}


def preprocess():
    global tex_lines_math_masked, tex_lines_math_texttt_masked
    env = list(set(re.findall(r"\\begin\{(\w+)\*?\}", tex)))
    for e in env:
        in_env[e] = []
        envs[e] = []
        current_start = -1
        begin_pattern = re.compile(r"\\begin\{%s\*?\}" % re.escape(e))
        end_pattern = re.compile(r"\\end\{%s\*?\}" % re.escape(e))
        for i, l in enumerate(tex_lines):
            if begin_pattern.search(l):
                in_env[e].append(True)
                current_start = i
            elif end_pattern.search(l):
                in_env[e].append(False)
                envs[e].append((current_start, i))
            else:
                if len(in_env[e]) == 0:
                    in_env[e].append(False)
                else:
                    in_env[e].append(in_env[e][-1])
            if "%" in tex_lines[i]:
                idx = tex_lines[i].index("%")
                if idx > 0 and tex_lines[i][idx - 1] != "\\":
                    tex_lines_clean[i] = tex_lines[i][0:max(0, (tex_lines[i].index("%") - 1))]
                    if tex_lines_clean[i].startswith("%"):
                        tex_lines_clean[i] = ""
                else:
                    tex_lines_clean[i] = tex_lines[i]
            else:
                tex_lines_clean[i] = tex_lines[i]
    if "comment" in in_env:
        for i in range(len(tex_lines_clean)):
            if in_env["comment"][i]:
                tex_lines_clean[i] = ""
    tex_lines_math_masked = mask_math_in_lines(tex_lines_clean)
    tex_lines_math_texttt_masked = [
        mask_spans(line, get_command_brace_spans(line, "texttt"))
        for line in tex_lines_math_masked
    ]


def in_any_env(line):
    for e in in_env:
        if in_env[e][line]:
            return True
    return False


def in_any_float(line):
    for f in FLOAT_ENVS:
        if f in in_env and in_env[f][line]:
            return True
    return False


def in_table(line):
    return (
        ("table" in in_env and in_env["table"][line])
        or ("tabular" in in_env and in_env["tabular"][line])
        or ("longtable" in in_env and in_env["longtable"][line])
    )


def in_code(line):
    if "lstlisting" in in_env:
        return in_env["lstlisting"][line]
    return False


def in_comment(line):
    if "comment" in in_env:
        return in_env["comment"][line]
    return False


def in_list(line):
    for env in ["itemize", "enumerate", "compactitem", "compactenum", "description"]:
        if env in in_env and in_env[env][line]:
            return True
    return False


def in_equation(line):
    if "equation" in in_env and in_env["equation"][line]:
        return True
    if "align" in in_env and in_env["align"][line]:
        return True
    if "align*" in in_env and in_env["align*"][line]:
        return True
    if "eqnarray" in in_env and in_env["eqnarray"][line]:
        return True
    if "theorem" in in_env and in_env["theorem"][line]:
        return True
    if "proof" in in_env and in_env["proof"][line]:
        return True
    if "proposition" in in_env and in_env["proposition"][line]:
        return True

    return False


def normalize_title(title):
    return re.sub("\\s+", " ", title.strip()).lower()


def get_section_entries():
    entries = []
    for i, l in enumerate(tex_lines_clean):
        match = SECTION_COMMAND_RE.search(l)
        if match:
            entries.append({
                "line": i,
                "level": match.group(1),
                "title": match.group(2).strip(),
                "span": match.span(),
                "starred": "*" in match.group(0)
            })
    return entries


def extract_section_ranges():
    entries = [e for e in get_section_entries() if e["level"] == "section"]
    appendix_line = tex.find("\\appendix")
    appendix_line = tex[:appendix_line].count("\n") if appendix_line != -1 else None
    ranges = []
    for idx, entry in enumerate(entries):
        start = entry["line"] + 1
        end = len(tex_lines)
        for nxt in entries[idx + 1:]:
            end = nxt["line"]
            break
        ranges.append({
            "title": entry["title"],
            "line": entry["line"],
            "span": entry["span"],
            "starred": entry["starred"],
            "in_appendix": appendix_line is not None and entry["line"] > appendix_line,
            "start": start,
            "end": end
        })
    return ranges


def count_sentences(text):
    stripped = text.strip()
    if not stripped:
        return 0
    cleaned = re.sub(r"\\[A-Za-z@]+\*?(?:\[[^\]]*\])?(?:\{[^{}]*\})?", " ", stripped)
    cleaned = mask_false_sentence_ends(cleaned)
    return len(SENTENCE_END_RE.findall(cleaned))


def is_paragraph_content(line):
    stripped = line.strip()
    if not stripped:
        return False
    if stripped.startswith("%"):
        return False
    if stripped in ["{", "}"]:
        return False
    if re.match("\\\\(section|subsection|subsubsection|paragraph|begin|end|item|label|caption|centering|bibliography|bibliographystyle|appendix|FloatBarrier|clearpage|cleardoublepage|maketitle|title|author|date|thanks|if|else|fi|newcommand|renewcommand|providecommand|definecolor|newlength|setlength|AtBeginDocument)\\b", stripped):
        return False
    return True


def get_paragraphs(start=0, end=None):
    if end is None:
        end = len(tex_lines_clean)
    paragraphs = []
    para_start = None
    para_lines = []
    for idx in range(start, end):
        if "document" in in_env and not in_env["document"][idx]:
            continue
        stripped = tex_lines_clean[idx].strip()
        if stripped.startswith("%"):
            continue
        if in_code(idx) or in_equation(idx) or in_any_float(idx):
            if para_lines:
                paragraphs.append((para_start, idx, " ".join(para_lines).strip()))
                para_start = None
                para_lines = []
            continue
        item_match = re.match(r"\\item(?:\s*\[[^\]]+\])?\s*(.*)", stripped)
        if item_match:
            if para_lines:
                paragraphs.append((para_start, idx, " ".join(para_lines).strip()))
            para_start = idx
            para_lines = [item_match.group(1).strip()]
            continue
        if is_paragraph_content(tex_lines_clean[idx]):
            if para_start is None:
                para_start = idx
            para_lines.append(stripped)
        else:
            if para_lines:
                paragraphs.append((para_start, idx, " ".join(para_lines).strip()))
                para_start = None
                para_lines = []
    if para_lines:
        paragraphs.append((para_start, end, " ".join(para_lines).strip()))
    return paragraphs


def get_math_spans(line):
    spans = []
    for pattern in [r"(?<!\\)(?<!\$)\$(?!\$)(.+?)(?<!\\)(?<!\$)\$(?!\$)", r"\\\((.+?)\\\)"]:
        for match in re.finditer(pattern, line):
            spans.append((match.start(), match.end(), match.group(1)))
    return spans


def mask_spans(line, spans):
    if not spans:
        return line
    chars = list(line)
    for start, end in spans:
        for pos in range(max(0, start), min(len(chars), end)):
            chars[pos] = " "
    return "".join(chars)


def get_inline_math_spans(line):
    spans = []
    i = 0
    math_start = None
    while i < len(line):
        if line.startswith(r"\(", i):
            if math_start is None:
                math_start = i
            else:
                spans.append((math_start, i + 2))
                math_start = None
            i += 2
            continue
        if line[i] == "$" and (i == 0 or line[i - 1] != "\\"):
            if math_start is None:
                math_start = i
            else:
                spans.append((math_start, i + 1))
                math_start = None
        i += 1
    return spans


def get_command_brace_spans(line, command):
    spans = []
    pattern = re.compile(r"\\%s\*?\{" % re.escape(command))
    for match in pattern.finditer(line):
        start = match.start()
        pos = match.end() - 1
        depth = 0
        while pos < len(line):
            char = line[pos]
            if char == "{":
                depth += 1
            elif char == "}":
                depth -= 1
                if depth == 0:
                    spans.append((start, pos + 1))
                    break
            pos += 1
    return spans


def mask_math_in_lines(lines):
    masked_lines = []
    state = None
    for line in lines:
        chars = list(line)
        i = 0
        while i < len(line):
            if state is None:
                if line.startswith(r"\(", i):
                    chars[i:i + 2] = [" ", " "]
                    state = r"\("
                    i += 2
                    continue
                if line.startswith(r"\[", i):
                    chars[i:i + 2] = [" ", " "]
                    state = r"\["
                    i += 2
                    continue
                if line.startswith("$$", i) and (i == 0 or line[i - 1] != "\\"):
                    chars[i:i + 2] = [" ", " "]
                    state = "$$"
                    i += 2
                    continue
                if line[i] == "$" and (i == 0 or line[i - 1] != "\\"):
                    chars[i] = " "
                    state = "$"
            else:
                if state == r"\(" and line.startswith(r"\)", i):
                    chars[i:i + 2] = [" ", " "]
                    state = None
                    i += 2
                    continue
                if state == r"\[" and line.startswith(r"\]", i):
                    chars[i:i + 2] = [" ", " "]
                    state = None
                    i += 2
                    continue
                if state == "$$" and line.startswith("$$", i) and (i == 0 or line[i - 1] != "\\"):
                    chars[i:i + 2] = [" ", " "]
                    state = None
                    i += 2
                    continue
                if state == "$" and line[i] == "$" and (i == 0 or line[i - 1] != "\\"):
                    chars[i] = " "
                    state = None
                    i += 1
                    continue
                chars[i] = " "
            i += 1
        masked_lines.append("".join(chars))
    return masked_lines


def mask_inline_math(line):
    return mask_spans(line, get_inline_math_spans(line))


def mask_inline_math_and_texttt(line):
    spans = get_inline_math_spans(line)
    spans.extend(get_command_brace_spans(line, "texttt"))
    return mask_spans(line, spans)
