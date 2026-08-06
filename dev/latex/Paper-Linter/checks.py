import re

import glossary_symbols as gs
import paperlint_actions as actions
import tex_structure as ts

CODE_WARNING_EXCEPTIONS = {
    "listing-alignment",
    "listing-label",
    "listing-caption",
    "listing-caption-order",
    "listing-float",
}

SENTENCE_END_RE = ts.SENTENCE_END_RE

math_text_mix_strings = set()
current_file_equation_symbols = set()
GLOSSARY_FILE = None
PREFIX_CHECK_CONFIGS = {}

def split_sentences(text):
    sentences = []
    start = 0
    for match in SENTENCE_END_RE.finditer(text):
        end = match.end()
        sentence = text[start:end].strip()
        if sentence:
            sentences.append(sentence)
        start = end
    tail = text[start:].strip()
    if tail:
        sentences.append(tail)
    return sentences


def build_sentence_sequence():
    sentences = []
    for start, end, text in ts.get_paragraphs():
        for sentence in split_sentences(text):
            sentences.append({"start": start, "end": end, "text": sentence})
    return sentences


def sentence_mentions_symbol(sentence_text, symbol):
    for _, _, content in ts.get_math_spans(sentence_text):
        if symbol in gs.extract_math_symbols(content):
            return True
    return False


def collect_equation_symbol_first_use():
    symbol_first_use = {}
    macro_block_depth = 0
    for i, line in enumerate(ts.tex_lines_clean):
        if ts.in_code(i):
            continue
        stripped = ts.strip_comment_from_line(line)
        if not stripped.strip():
            continue
        stripped_no_ws = stripped.strip()
        if macro_block_depth > 0:
            macro_block_depth += stripped_no_ws.count("{") - stripped_no_ws.count("}")
            if macro_block_depth <= 0:
                macro_block_depth = 0
            continue
        if stripped_no_ws.startswith(("\\newcommand", "\\providecommand", "\\def", "\\gdef", "\\xdef", "\\makeatletter")):
            macro_block_depth = stripped_no_ws.count("{") - stripped_no_ws.count("}")
            if macro_block_depth <= 0:
                macro_block_depth = 0
            continue
        math_contents = [content for _, _, content in ts.get_math_spans(stripped)]
        if ts.in_equation(i):
            equation_line = re.sub(r"\\(?:begin|end)\{[^{}]+\}", " ", stripped)
            equation_line = re.sub(r"\\(?:label|tag|nonumber|notag)\*?(?:\{[^{}]*\})?", " ", equation_line)
            math_contents.append(equation_line)
        for content in math_contents:
            for symbol in gs.extract_math_symbols(content):
                if symbol not in symbol_first_use:
                    symbol_first_use[symbol] = i
    return symbol_first_use


def check_equation_symbols_defined():
    global current_file_equation_symbols
    warns = []
    symbol_first_use = collect_equation_symbol_first_use()
    current_file_equation_symbols = set(symbol_first_use.keys())

    if not symbol_first_use:
        return warns

    sentences = build_sentence_sequence()
    if not sentences:
        return warns

    for symbol, first_line in symbol_first_use.items():
        sentence_index = None
        for idx, sentence in enumerate(sentences):
            if sentence["start"] <= first_line <= sentence["end"]:
                sentence_index = idx
                break
        if sentence_index is None:
            continue

        nearby = sentences[max(0, sentence_index - 3):min(len(sentences), sentence_index + 4)]
        if any(sentence_mentions_symbol(sentence["text"], symbol) for sentence in nearby):
            continue

        if symbol.startswith("\\"):
            symbol_text = symbol
        else:
            symbol_text = "$%s$" % symbol
        warns.append((first_line, "Symbol %s is used in an equation but is not mentioned inline in the surrounding 3 sentences" % symbol_text, (0, 0)))
    return warns


def _count_gls_commands(text):
    return len(re.findall(r"\\gls\{[^{}]+\}", text))


def check_math_glossary_reference_coverage():
    warns = []

    for start, end in ts.envs.get("equation", []):
        equation_text = "\n".join(ts.strip_comment_from_line(ts.tex_lines_clean[i]) for i in range(start, end + 1))
        symbol_count = len(gs.extract_math_symbols(equation_text))
        if symbol_count >= 2 and _count_gls_commands(equation_text) < 2:
            warns.append((start, "Equation environment should include at least two \\gls references", (0, 0)))
        elif symbol_count == 1 and _count_gls_commands(equation_text) < 1:
            warns.append((start, "Equation environment should include at least one \\gls reference", (0, 0)))

    for i, line in enumerate(ts.tex_lines_clean):
        for match in re.finditer(r"(?<!\\)\$\$(.+?)(?<!\\)\$\$", line):
            content = match.group(1)
            if gs.extract_math_symbols(content) and _count_gls_commands(content) < 1:
                warns.append((i, "Display math $$...$$ should include at least one \\gls reference", match.span()))

    return warns


def check_prefix():
    config = PREFIX_CHECK_CONFIGS.get("prefix", {})
    path = config.get("path")
    prefix = config.get("prefix")
    source_prefix = config.get("source_prefix")
    recursive = config.get("recursive", False)
    if not path or not prefix:
        return []
    if actions.preview_prefix(path, prefix, recursive=recursive, source_prefix=source_prefix) > 0:
        return [(-1, "Prefixing labels and references is required")]
    return []


def check_mathmode_subscripts():
    warns = []
    updated = actions.replace_mathmode_subscripts(ts.tex)
    if updated == ts.tex:
        return warns
    for i, (old_line, new_line) in enumerate(zip(ts.tex_lines, updated.split("\n"))):
        if old_line != new_line:
            warns.append((i, "Math subscripts should wrap prose words in \\text{...}", (0, len(old_line))))
            break
    return warns


def check_glossary_refs():
    warns = []
    if GLOSSARY_FILE is None:
        return warns
    replacements = gs.read_glossary_replacements(GLOSSARY_FILE)
    if not replacements:
        return warns
    for i, line in enumerate(ts.tex_lines_clean):
        _, count = gs.replace_glossary_refs_in_line(line, replacements, equation_line=ts.in_equation(i))
        if count > 0:
            warns.append((i, "Math symbols should use glossary references where available", (0, len(line))))
            break
    return warns


def check_space_before_cite():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        updated, span = actions.space_before_cite_line(l)
        if span is not None:
            warns.append((i, "No space before \\cite", span))
    return warns

def check_float_alignment(env):
    warns = []
    for i, l in enumerate(ts.tex_lines):
        b = re.search(r"\\begin\{%s\}" % env, l)
        if b:
            if not re.search(r"%s}\[[^\]]*[htbH][^\]]*\]" % env, l):
                warns.append((i, "%s without alignment: %s" % (env, l.strip()), b.span()))
    return warns

def check_figure_alignment():
    return check_float_alignment("figure")

def check_table_alignment():
    return check_float_alignment("table")

def check_listing_alignment():
    return check_float_alignment("listing")

def _check_float_has_command(env, pattern, message):
    warns = []
    if env not in ts.envs:
        return warns
    for r in ts.envs[env]:
        found = False
        for i in range(*r):
            if re.search(pattern, ts.tex_lines[i]):
                found = True
        if not found:
            warns.append((r[0], message % env))
    return warns

def check_float_has_label(env):
    return _check_float_has_command(env, r"\\label\{", "%s without a label")


def check_float_has_caption(env):
    return _check_float_has_command(env, r"\\caption\{", "%s without a caption")

def check_float_caption_label_order(env):
    warns = []
    if env not in ts.envs: return warns
    for r in ts.envs[env]:
        label = -1
        caption = -1
        for i in range(*r):
            b = re.search(r"\\caption\{", ts.tex_lines[i])
            if b:
                caption = i
            b = re.search(r"\\label\{", ts.tex_lines[i])
            if b:
                label = i
        if label > -1 and caption > -1 and label < caption:
            warns.append((r[0], "label before caption in %s, swap for correct references" % env))
    return warns


def check_no_resizebox_for_tables():
    warns = []
    if "table" not in ts.envs: return warns
    for r in ts.envs["table"]:
        rb = False
        b = None
        for i in range(*r):
            b = re.search(r"\\resizebox\{", ts.tex_lines[i])
            if b:
                rb = True
                break
        if rb:
            warns.append((r[0], "table with resizebox -> use adjustbox instead"))
    return warns


def check_weird_units():
    warns = []
    block = ["\\textwidth", "\\linewidth"]
    for i, l in enumerate(ts.tex_lines):
        for b in block:
            if b in l:
                warns.append((i, "use \\hsize instead of %s" % b, (l.index(b), l.index(b) + len(b))))
    return warns

def check_figure_has_label():
    return check_float_has_label("figure")

def check_table_has_label():
    return check_float_has_label("table")

def check_listing_has_label():
    return check_float_has_label("listing")

def check_figure_has_caption():
    return check_float_has_caption("figure")

def check_table_has_caption():
    return check_float_has_caption("table")

def check_listing_has_caption():
    return check_float_has_caption("listing")

def check_figure_caption_label_order():
    return check_float_caption_label_order("figure")

def check_table_caption_label_order():
    return check_float_caption_label_order("table")

def check_listing_caption_label_order():
    return check_float_caption_label_order("listing")

def check_todos():
    warns = []
    for i, l in enumerate(ts.tex_lines_clean):
        if "TODO" in l:
            warns.append((i, "TODO found", (l.index("TODO"), l.index("TODO") + 4)))
    return warns


def check_notes():
    warns = []
    for i, l in enumerate(ts.tex_lines_clean):
        if "\\note" in l:
            warns.append((i, "\\note found", (l.index("\\note"), l.index("\\note") + 5)))
        if "\\todo" in l:
            warns.append((i, "\\todo found", (l.index("\\todo"), l.index("\\todo") + 5)))
    return warns


def check_math_numbers():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        if ts.in_code(i):
            continue
        n = re.search("\\$\\d+\\$", ts.tex_lines[i]) 
        if n and not ts.in_any_float(i):
            warns.append((i, "Number in math mode, consider using siunit instead", n.span()))
    return warns


def check_math_text_mix():
    global math_text_mix_strings
    warns = []
    for i, l in enumerate(ts.tex_lines_clean):
        if ts.in_code(i):
            continue
        for start, end, content in ts.get_math_spans(l):
            normalized = gs.strip_explicit_math_text(content)
            normalized = re.sub(r"\\(?:begin|end)\*?\{[^{}]*\}", " ", normalized)
            normalized = re.sub("\\\\[A-Za-z]+", " ", normalized)
            normalized = re.sub("[\\^_{}&=+\\-*/(),.;:0-9\\s]", " ", normalized)
            words = re.findall("[A-Za-z]{2,}", normalized)
            prose_words = [word for word in words if word.islower()]
            prose_like_words = [word for word in prose_words if len(word) >= 4]
            if prose_like_words:
                for word in prose_like_words:
                    math_text_mix_strings.add(word.lower())
                warns.append((i, "Text-like content in math mode, mark it explicitly with \\text{...}", (start, end)))
                break
    return warns


def check_units():
    warns = []
    unit_pattern = re.compile(r"(?<!\\)(\b\d+(?:\.\d+)?)(\s?)(ms|ns|us|s|min|h|Hz|kHz|MHz|GHz|B|KB|MB|GB|TB|bit|m|cm|mm|km|g|kg|mg|V|A|W|kW|J|N|K|mol|cd)\b")
    for i, l in enumerate(ts.tex_lines_clean):
        if ts.in_code(i):
            continue
        if any(cmd in l for cmd in ["\\SI{", "\\si{", "\\qty{", "\\unit{"]):
            continue
        masked = ts.tex_lines_math_masked[i]
        match = unit_pattern.search(masked)
        if match:
            before = masked[match.start() - 1] if match.start() > 0 else ""
            after = masked[match.end()] if match.end() < len(masked) else ""
            if before.isalnum() or before in ["-", "_"]:
                continue
            if after.isalnum() or after in ["-", "_"]:
                continue
            prefix = masked[:match.start()].rstrip()
            if re.search(r"\\[A-Za-z@]+\*?(?:\{[^{}]*\})*\{$", prefix):
                continue
            if re.search(r"(?:^|[,\[])\s*[A-Za-z@][A-Za-z0-9@ _-]*=\s*$", prefix):
                continue
            if match.group(2) == "":
                warns.append((i, "Unit attached to number, separate it and prefer siunitx", match.span()))
            else:
                warns.append((i, "Number with unit should use siunitx", match.span()))
    return warns


def check_float_placement():
    warns = []
    for i, l in enumerate(ts.tex_lines_clean):
        match = re.search("\\\\begin\\{(%s)\\}\\[([^\\]]*[hH][^\\]]*)\\]" % "|".join(ts.FLOAT_ENVS), l)
        if match:
            warns.append((i, "Float uses discouraged placement specifier [%s]" % match.group(2), match.span(2)))
    return warns


def check_unused_macro():
    warns = []
    for name, line, span in actions.collect_unused_macro_definitions(ts.tex):
        if len(re.findall("\\\\%s\\b" % re.escape(name), ts.tex)) <= 1:
            warns.append((line, "Macro \\%s is defined but never used" % name, span))
    return warns


def check_required_sections():
    warns = []
    if not ts.is_document_root:
        return warns
    sections = {ts.normalize_title(section["title"]) for section in ts.extract_section_ranges()}
    if "\\begin{abstract}" in ts.tex:
        sections.add("abstract")
    required = ["abstract", "introduction", "conclusion"]
    missing = [section.title() for section in required if section not in sections]
    if missing:
        warns.append((-1, "Missing required paper components: %s" % ", ".join(missing)))
    return warns


def check_section_balance():
    warns = []
    sections = ts.extract_section_ranges()
    if len(sections) < 2:
        return warns
    word_counts = []
    filtered_sections = [section for section in sections if not section["starred"] and not section["in_appendix"]]
    if len(filtered_sections) < 2:
        return warns
    for section in filtered_sections:
        text = " ".join(ts.tex_lines_clean[section["start"]:section["end"]])
        word_counts.append(len(re.findall("\\b\\w+\\b", text)))
    sorted_counts = sorted(word_counts)
    median = sorted_counts[len(sorted_counts) // 2]
    if median == 0:
        return warns
    for section, count in zip(filtered_sections, word_counts):
        if count < max(40, int(median * 0.3)):
            warns.append((section["line"], "Section '%s' is extremely short compared to the rest of the paper" % section["title"], section["span"]))
        elif count > max(600, int(median * 2.5)):
            warns.append((section["line"], "Section '%s' is extremely long compared to the rest of the paper" % section["title"], section["span"]))
    return warns


def check_deprecated():
    warns = []
    patterns = [
        (re.compile("\\{\\\\bf\\b"), "Deprecated font switch \\bf, use \\textbf{...}"),
        (re.compile("\\{\\\\it\\b"), "Deprecated font switch \\it, use \\textit{...}"),
        (re.compile("\\{\\\\rm\\b"), "Deprecated font switch \\rm, use \\textrm{...}"),
        (re.compile("\\{\\\\sc\\b"), "Deprecated font switch \\sc, use \\textsc{...}"),
        (re.compile("\\{\\\\tt\\b"), "Deprecated font switch \\tt, use \\texttt{...}"),
        (re.compile("\\$\\$"), "Deprecated display math $$...$$, use \\[...\\] or an equation environment"),
        (re.compile("\\\\centerline\\b"), "Deprecated \\centerline command"),
        (re.compile("\\\\epsfig\\b"), "Deprecated \\epsfig command, use \\includegraphics")
    ]
    for i, l in enumerate(ts.tex_lines_clean):
        for pattern, message in patterns:
            match = pattern.search(l)
            if match:
                warns.append((i, message, match.span()))
    return warns


def check_package_conflict():
    warns = []
    packages = {}
    for i, l in enumerate(ts.tex_lines_clean):
        for match in re.finditer("\\\\usepackage(?:\\[[^\\]]+\\])?\\{([^\\}]+)\\}", l):
            for package in [p.strip() for p in match.group(1).split(",")]:
                packages[package] = i
    conflicts = [
        ("subfigure", "subcaption"),
        ("subfig", "subcaption"),
        ("natbib", "biblatex")
    ]
    for left, right in conflicts:
        if left in packages and right in packages:
            warns.append((packages[right], "Conflicting packages loaded together: %s and %s" % (left, right)))
    return warns


def check_duplicate_word():
    warns = []
    for i, l in enumerate(ts.tex_lines_clean):
        if ts.in_code(i):
            continue
        updated, span = actions.duplicate_word_line(l)
        if span is not None:
            match = re.search(r"\b(\w+)\s+\1\b", l, re.IGNORECASE)
            if match is not None:
                warns.append((i, "Duplicated adjacent word '%s'" % match.group(1), span))
            else:
                warns.append((i, "Duplicated adjacent word", span))
    return warns


def check_paragraph_length():
    warns = []
    for start, end, text in ts.get_paragraphs():
        if ts.in_list(start):
            continue
        if not re.search("[.!?]", text):
            continue
        sentence_count = ts.count_sentences(text)
        if sentence_count == 2:
            warns.append((start, "Paragraph has only %d sentence%s; expected at least 3" % (sentence_count, "" if sentence_count == 1 else "s"), (0, len(ts.tex_lines[start]))))
    return warns


def check_section_length():
    warns = []
    for section in ts.extract_section_ranges():
        if section["starred"] or section["in_appendix"]:
            continue
        paragraphs = [p for p in ts.get_paragraphs(section["start"], section["end"]) if p[0] >= section["start"] and p[0] < section["end"]]
        if len(paragraphs) < 3:
            warns.append((section["line"], "Section '%s' has only %d paragraph%s; expected at least 3" % (section["title"], len(paragraphs), "" if len(paragraphs) == 1 else "s"), section["span"]))
    return warns


def check_subsection_count_sanity():
    warns = []
    sections = ts.extract_section_ranges()
    for idx, section in enumerate(sections):
        if section["starred"] or section["in_appendix"]:
            continue
        next_line = sections[idx + 1]["line"] if idx + 1 < len(sections) else len(ts.tex_lines)
        count = 0
        for line_no in range(section["start"], next_line):
            if re.search("\\\\subsection\\*?\\{", ts.tex_lines_clean[line_no]):
                count += 1
        if count == 1 or count > 10:
            warns.append((section["line"], "Section '%s' has %d subsection%s; expected 0 or between 2 and 10" % (section["title"], count, "" if count == 1 else "s"), section["span"]))
    return warns


def check_large_numbers_without_si():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        n = re.search(r"[\s(]\d{5,}[\s),\.]", ts.tex_lines[i]) 
        if n and not ts.in_any_float(i):
            warns.append((i, "Large number without formating, consider using siunit", n.span()))
    return warns

def check_env_not_in_float(env, float_env):
    warns = []
    if env in ts.envs:
        for e in ts.envs[env]:
            if (float_env not in ts.in_env) or (not ts.in_env[float_env][e[0]]):
                warns.append((e[0], "%s not within %s environment" % (env, float_env)))
    return warns
    

def check_listing_in_correct_float():
    return check_env_not_in_float("lstlisting", "listing")

def check_tabular_in_correct_float():
    return check_env_not_in_float("tabular", "table")

def check_tikz_in_correct_float():
    return check_env_not_in_float("tikzpicture", "figure")


def check_comment_has_space():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        updated, span = actions.comment_has_space_line(l)
        if span is not None and not ts.in_code(i):
            warns.append((i, "Comment without a whitespace before", span))
    return warns


def check_percent_without_siunix():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        n = re.search("\\d+\\s*\\\\%", l)
        if n:
            warns.append((i, "Number with percent without siunit", n.span(0)))
    return warns


def check_short_form():
    warns = []
    for i, l in enumerate(ts.tex_lines_clean):
        if ts.in_comment(i) or l.strip().startswith("%"):
            continue
        n = re.search("[^`%]\\w+'[a-rt-z]", l)
        if n:
            warns.append((i, "Contracted form used", n.span()))
    return warns


def check_labels_referenced():
    warns = []
    labels = [] #re.findall("\\\\label\{([^\\}]+)\}", ts.tex)
    for i, l in enumerate(ts.tex_lines_clean):
        lab = re.search(r"\\label\{([^\}]+)\}", l)
        if lab:
            labels.append((lab.group(1), i, lab.span()))

    referenced_labels = set()
    reference_pattern = re.compile(
        r"\\(?:[cC]ref|[vV]ref|ref|eqref|autoref|nameref|labelcref|namecref|lcnamecref|fref|Fref|fullref)\*?\{([^\}]*)\}"
    )
    range_pattern = re.compile(
        r"\\(?:[cC]refrange|[cC]labelcrefrange)\*?\{([^\}]*)\}\{([^\}]*)\}"
    )

    for line in ts.tex_lines_clean:
        for match in reference_pattern.finditer(line):
            for ref_label in match.group(1).split(","):
                ref_label = ref_label.strip()
                if ref_label:
                    referenced_labels.add(ref_label)
        for match in range_pattern.finditer(line):
            start_label = match.group(1).strip()
            end_label = match.group(2).strip()
            if start_label:
                referenced_labels.add(start_label)
            if end_label:
                referenced_labels.add(end_label)

    def subfigure_parent_referenced(label_line):
        if "subfigure" not in ts.in_env or not ts.in_env["subfigure"][label_line]:
            return False
        if "figure" not in ts.envs:
            return False
        parent_range = None
        for fig_range in ts.envs["figure"]:
            if fig_range[0] <= label_line <= fig_range[1]:
                parent_range = fig_range
                break
        if parent_range is None:
            return False
        for label_name, other_line, _ in labels:
            if other_line < parent_range[0] or other_line > parent_range[1]:
                continue
            if "subfigure" in ts.in_env and ts.in_env["subfigure"][other_line]:
                continue
            if label_name in referenced_labels:
                return True
        return False

    for lab in labels:
        if "sec:" in lab[0] or "eq:" in lab[0]:
            continue
        if lab[0] not in referenced_labels:
            if subfigure_parent_referenced(lab[1]):
                continue
            if not (lab[0].startswith("sec") or lab[0].startswith("subsec")):
                warns.append((lab[1], "Label %s is not referenced" % lab[0], lab[2]))
    return warns


def check_section_capitalization():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        updated, span = actions.section_capitalization_line(l)
        if span is not None:
            warns.append((i, "Wrong capitalization of header", span))
    return warns


def check_quotation():
    warns = []
    for i, l in enumerate(ts.tex_lines_clean):
        if ts.in_code(i):
            continue
        updated, span = actions.quotation_marks_line(l)
        if span is not None:
            warns.append((i, "Wrong quotation, use `` and '' instead of \"", span))
    return warns


def check_hline_in_table():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        hl = re.search("\\\\hline", l)
        if hl:
            if "tabular" in ts.in_env and ts.in_env["tabular"]:
                warns.append((i, "\\hline in table, consider using \\toprule, \\midrule, \\bottomrule.", hl.span()))
    return warns


def check_space_before_punctuation():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        updated, span = actions.space_before_punctuation_line(l)
        if span is not None and not ts.in_any_env(i):
            warns.append((i, "Spacing before punctuation", span))
    return warns


def check_headers_without_text():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        n = re.search("(section|paragraph)\\{([^\\}]+)\\}", l)
        if n:
            nx = i
            while (nx + 1) < len(ts.tex_lines):
                nx += 1
                if len(ts.tex_lines[nx].strip()) == 0: continue
                if ts.tex_lines[nx].strip().startswith("%"): continue
                nn = re.search("(section|paragraph)\\{([^\\}]+)\\}", ts.tex_lines[nx])
                if nn:
                    warns.append((i, "Section header without text before next header", n.span()))
                break
    return warns


def check_one_sentence_paragraphs():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        if i > 0 and i < len(ts.tex_lines) - 1:
            if len(ts.tex_lines[i - 1].strip()) == 0 and len(ts.tex_lines[i + 1].strip()) == 0 and len(ts.tex_lines[i].strip()) > 0:
                if ts.tex_lines[i].strip().startswith("\\"): continue
                if ". " in ts.tex_lines[i]: continue
                warns.append((i, "One-sentence paragraph", (0, len(ts.tex_lines[i]))))
    return warns


def check_multiple_sentences_per_line():
    warns = []
    for i, l in enumerate(ts.tex_lines_clean):
        if ts.in_code(i) or ts.in_comment(i) or l.strip().startswith("%"):
            continue
        if "vs." in l.rstrip():
            continue
        masked = ts.mask_false_sentence_ends(l.rstrip())
        if ts.count_sentences(masked) > 1:
            p = SENTENCE_END_RE.search(masked)
            warns.append((i, "Multiple sentences in one line", p.span() if p else (0, len(l))))
    return warns


def check_unbalanced_brackets():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        if ts.in_code(i):
            continue
        masked = ts.tex_lines_math_masked[i]
        if masked.count("(") != masked.count(")"):
            first = min(masked.index("(") if masked.count("(") > 0 else len(masked), masked.index(")") if masked.count(")") > 0 else len(masked))
            last = max(masked.rindex("(") if masked.count("(") > 0 else len(masked), masked.rindex(")") if masked.count(")") > 0 else len(masked))
            warns.append((i, "Mismatch of opening and closing parenthesis", (first, last)))
    return warns


def check_and_or():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        ao = re.search("and/or", l)
        if ao:
            warns.append((i, "And/or discouraged in academic writing", ao.span()))
    return warns


def check_ellipsis():
    warns = []
    for i, l in enumerate(ts.tex_lines_clean):
        el = re.search("\\w+\\.\\.\\.", l)
        if el:
            warns.append((i, "Ellipsis \"...\" discouraged in academic writing", el.span()))
    return warns


def check_etc():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        el = re.search("\\s+etc[\\.\\w]", l)
        if el:
            warns.append((i, "Unspecific \"etc\" discouraged in academic writing", el.span()))
    return warns


def check_footnote():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        fn = re.search("\\s*\\\\footnote\\{[^\\}]+\\}\\.", l)
        if fn:
            warns.append((i, "Footnote must be after the full stop", fn.span()))
    return warns


def check_table_top_caption():
    warns = []
    if "table" in ts.envs:
        for table in ts.envs["table"]:
            caption = -1
            tab = -1
            for intab in range(*table):
                if re.search("\\\\caption\\{", ts.tex_lines[intab]):
                   caption = intab
                if re.search("\\\\begin\\{tabular", ts.tex_lines[intab]):
                    tab = intab
            if tab != -1 and caption != -1 and tab < caption:
                warns.append((table[0], "Table caption must be above table"))
    return warns



def check_punctuation_end_of_line():
    warns = []
    for i, l in enumerate(ts.tex_lines_clean):
        sl = l.strip()
        if len(sl) < 10: continue
        if len(sl.split(" ")) < 8: continue
        if ts.in_any_float(i): continue
        if "lstlisting" in ts.in_env and ts.in_env["lstlisting"][i]: continue
        if sl.startswith("\\") or sl.startswith("%"): continue
        if sl.endswith("\\\\") or sl.endswith("}"): continue
        if sl.endswith(".") or sl.endswith("!") or sl.endswith("?") or sl.endswith(":") or sl.endswith(";"): continue
        p = re.search("\\s*[\\w})$]+[\\.!?}{:;\\\\]\\s*$", l.rstrip())
        if not p:
            warns.append((i, "Line ends without punctuation", (len(l) - 2, len(l))))
    return warns


def check_table_vertical_lines():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        t = re.search("\\\\begin\\{tabular\\}\\{([^\\}]+)\\}", l)
        if t and "|" in t.group(1):
            warns.append((i, "Vertical lines in tables are discouraged", t.span()))
    return warns


def check_will():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        w = re.search("\\s+will\\s+", l)
        if w:
            warns.append((i, "Usage of \"will\" is discouraged.", w.span()))
    return warns


def check_subsection_count():
    warns = []
    last_section = -1
    subsections = []
    for i, l in enumerate(ts.tex_lines):
        if re.search("\\\\section{", l):
            if last_section != -1 and len(subsections) == 1:
                warns.append((last_section, "Section only has one subsection", re.search("\\\\section{", ts.tex_lines[last_section]).span()))
            last_section = i
            subsections = []
        if re.search("\\\\subsection{", l):
            subsections.append(i)
    return warns


def check_mixed_compact_and_item():
    warns = []
    if "\\begin{compactenum}" in ts.tex:
        for i, l in enumerate(ts.tex_lines):
            it = re.search(r"\\begin\{enumerate\}", l)
            if it:
                warns.append((i, "compactenum mixed with enumerate", it.span()))
    if "\\begin{compactitem}" in ts.tex:
        for i, l in enumerate(ts.tex_lines):
            it = re.search(r"\\begin\{itemize\}", l)
            if it:
                warns.append((i, "compactitem mixed with itemize", it.span()))
    return warns


def check_center_in_float():
    warns = []
    if "center" in ts.envs:
        for c in ts.envs["center"]:
            if ts.in_any_float(c[0]):
                warns.append((c[0], "Use \\centering instead of \\begin{center} inside floats", re.search(r"\\begin\{center\}", ts.tex_lines[c[0]]).span()))
    return warns


def check_appendix():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        ap = re.search(r"\\begin\{appendix\}", l)
        if ap:
            warns.append((i, "Use \\appendix instead of \\begin{appendix}", ap.span()))
    return warns


def check_eqnarray():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        ap = re.search(r"\\begin\{eqnarray\}", l)
        if ap:
            warns.append((i, "Use \\begin{align} instead of \\begin{eqnarray}", ap.span()))
    return warns


def check_acm_pc():
    # based on https://www.acm.org/diversity-inclusion/words-matter
    warns = []
    replace = [
        ("\\bsupremacy\\b", "advantage"),
        ("\\bmaster\\b", "main/primary/leader/parent/host"),
        ("\\bslave\\b", "secondary/replica/follower/child/worker/client"),
        ("\\bhe\\b", "they"),
        ("\\bshe\\b", "they"),
        ("\\bhis\\b", "their"),
        ("\\bhers?\\b", "their/them"),
        ("\\bhim\\b", "them"),
        ("\\bmale\\bconnector\\b", "plug"),
        ("\\bfemale\\bconnector\\b", "socket"),
        ("\\bblind\\b", "anonymous"),
        ("\\bblack\\-?\\s?list\\b", "blocklist/unapprovedlist"),
        ("\\bwhite\\-?\\s?list\\b", "allowlist/approvedlist"),
        ("\\bblack\\-?\\s?hat\\b", "unethical attacker/hostile force"),
        ("\\bwhite\\-?\\s?hat\\b", "ethical attacker/friendly force"),
        ("\\bblack\\-?\\s?box\\b", "opaque box"),
        ("\\bwhite\\-?\\s?box\\b", "clear box"),
        ("\\baverage\\s?user\\b", "common/standard/typical user"),
        ("\\babort\\s?child\\b", "cancel/force quit/stop/end/finalize"),
        ("\\bterminate\\s?child\\b", "cancel/force quit/stop/end/finalize"),
        ("\\bdark\\-?\\s?pattern\\b", "deceptive design"),
        ("\\bdummy\\-?\\s?head\\b", "temporary head"),
        ("\\bgender\\-?\\s?bender\\b", "plug-socket adapter"),
        ("\\borphaned\\-?\\s?object\\b", "unreferenced/unlinked object"),
        ("\\bsanity\\-?\\s?check", "coherence/quick/well-formedness check")
    ]
    for i, l in enumerate(ts.tex_lines):
        for r in replace:
            w = re.search(r[0], l)
            if w:
                warns.append((i, "Discouraged term \"%s\", consider replacing with \"%s\"" % (w.group(), r[1]), w.span()))
    return warns


def check_cite_noun():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        if re.search(r"\\citet\*?\{", l):
            continue
        ap = re.search("\\b(in|from|by|and|or)[\\s~]\\\\cite", l.lower())
        if ap:
            warns.append((i, "Citation is used as noun", ap.span()))
        ap = re.search("^\\s*\\\\cite", l)
        if ap and not ts.in_table(i):
            warns.append((i, "Citation at the beginning of a sentence (probably as noun)", ap.span()))
    return warns


def check_cite_duplicate():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        updated, span = actions.cite_duplicate_line(l)
        if span is not None:
            warns.append((i, "Duplicate citation key", span))
    return warns


def check_multicite():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        updated, span = actions.multicite_line(l)
        if span is not None:
            warns.append((i, "Multiple \\cite commands, use multiple citation keys in one \\cite instead", span))
    return warns


def check_emptycite():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        cites = re.search("\\\\citeA?\\{\\s*\\}", l)
        if cites:
            warns.append((i, "Empty citation key", cites.span()))
    return warns

def check_conjunction_start():
    warns = []
    for i, l in enumerate(ts.tex_lines_clean):
        p = re.search("[\\.!?]\\s+(And|Or|But)[\\s,]", l.rstrip())
        if p:
            warns.append((i, "Starting a sentence with a conjunction is discouraged", p.span()))
        p = re.search("^(And|Or|But)[\\s,]", l.rstrip())
        if p:
            warns.append((i, "Starting a sentence with a conjunction is discouraged", p.span()))
    return warns


def check_brackets_space():
    warns = []
    for i, l in enumerate(ts.tex_lines_clean):
        if ts.in_code(i) or ts.in_equation(i) or (len(l.strip()) > 0 and l.strip()[0] in ["\\", "%"]): continue
        updated, span = actions.brackets_space_line(l)
        if span is not None:
            warns.append((i, "There must be a space before an opening parenthesis", span))
    return warns  


def check_acronym_capitalization():
    warns = []
    acronyms = set()
    acronym_first = {}
    for i, l in enumerate(ts.tex_lines_clean):
        if ts.in_code(i): continue
        for p in re.finditer("\\b[A-Z]{3,}\\b", l):
            acronym = p.group()
            pos = p.span()[0]
            if pos > 0 and l[pos - 1] == '\\':
                continue
            if acronym not in acronyms:
                acronyms.add(acronym)
                acronym_first[acronym] = i
    for i, l in enumerate(ts.tex_lines_clean):
        if ts.in_code(i): continue
        if "@" in l:
            continue
        for p in re.finditer("\\b[A-Za-z]{3,}\\b", l):
            found = p.group()
            if p.span()[0] > 0 and l[p.span()[0] - 1] == '\\':
                continue # probably a macro
            if l[:p.span()[0]].count("{") != l[:p.span()[0]].count("}"):
                continue # probably inside a reference or label
            canonical = found.upper()
            if found.endswith("s") and found[:-1].upper() in acronyms:
                found = found[:-1]
                canonical = found.upper()
            if canonical not in acronyms:
                continue
            if acronym_first[canonical] >= i:
                continue
            if found.islower():
                continue
            if not found.isupper():
                warns.append((i, "(Potential) acronym with wrong capitalization (first defined in Line %d)" % (acronym_first[canonical] + 1), p.span()))
    return warns  

def check_numeral():
    warns = []
    for i, l in enumerate(ts.tex_lines):
        updated, span = actions.numerals_line(l)
        if span is not None:
            warns.append((i, "Numeral should be replaced with digits", span))
    return warns


def check_colors():
    warns = []
    cols = [
        "\\bred\\b",
        "\\bgreen\\b",
        "\\bblue\\b",
        "\\byellow\\b",
        "\\borange\\b",
        "\\bmagenta\\b",
        "\\bcyan\\b",
        "\\bbrown\\b",
        "\\bpink\\b"        
    ]
    modifiers = [
        "\\bdott?(ed)?\\b",
        "\\bdash(ed)?\\b",
        "\\bthick\\b",
        "\\bthin\\b",
        "\\bdash-?dotted\\b",
        "\\bhatch",
        "\\bcross",
        "\\bcheck",
        "\\bpattern"
    ]
    for i, l in enumerate(ts.tex_lines):
        for c in cols:
            w = re.search(c, l)
            if w:
                prefix = l[:w.span()[0]].rstrip()
                # check for = or { in front of color
                if w.span()[0] > 0 and (l[w.span()[0] - 1] == "=" or l[w.span()[0] - 1] == "{"): continue
                if re.search(r"\\[A-Za-z@]+\*?(?:\{[^{}]*\})*\[[^\]]*$", prefix):
                    continue
                if prefix.endswith("[") or prefix.endswith(","):
                    continue
                # reduce false positives by looking for modifiers
                mod = False
                for m in modifiers:
                    if re.search(m, l):
                        mod = True
                        break
                if not mod:
                    warns.append((i, "Colors (\"%s\") without a modifier such as dashed/dotted/... should be avoided." % (w[0]), w.span()))
    return warns


def check_inconsistent_word_style():
    warns = []
    word_style = {}
    for i, l in enumerate(ts.tex_lines_clean):
        styled = re.search(r"\\text([^\{]+)\{([^\}]+)\}", l)
        if styled and "newcommand" not in l:
            if styled[2] in word_style:
                if styled[1] != word_style[styled[2]][1][1]:
                    warns.append((i, "Word '%s' is styled inconsistently, used with \\text%s before at line %d" % (styled[2], word_style[styled[2]][1][1], word_style[styled[2]][0] + 1), styled.span()))
            else:
                word_style[styled[2]] = (i, styled)
    return warns


def check_missing_word_style():
    warns = []
    style_for_word, style_counts, styled_word_pattern = actions.collect_word_style_data(ts.tex)
    for i, l in enumerate(ts.tex_lines_clean):
        if ts.in_code(i):
            continue
        updated, span = actions.missing_word_style_line(l, style_for_word, style_counts, styled_word_pattern)
        if span is not None:
            warns.append((i, "Word used without a style", span))
    return warns


def print_warnings(warn, writer, use_color=True, output=True):
    warnings = 0
    sorted_warn = sorted(warn, key=lambda tup: tup[0][0])
    for cw in sorted_warn:
        w = cw[0]
        if w[0] != -1 and ts.tex_lines[w[0]].strip().startswith("%"):
            continue
        if w[0] != -1 and ts.in_comment(w[0]):
            continue
        if w[0] != -1 and ts.in_code(w[0]) and cw[1] not in CODE_WARNING_EXCEPTIONS:
            continue

        if output:
            if use_color:
                writer("\033[33mWarning %d\033[0m: " % (warnings + 1), end="")
            else:
                writer("Warning %d: " % (warnings + 1), end="")
        warnings += 1
        if w[0] != -1:
            if output:
                writer("Line %d: %s" % (w[0] + 1, w[1]), end="")
        else:
            if output:
                writer(w[1], end="")
        
        if output:
            if use_color:
                writer("  \033[90m[%s]\033[0m" % cw[1], end="")
            else:
                writer("  [%s]" % cw[1], end="")
            writer("")

        if len(w) > 2:
            if output:
                writer("    %s" % ts.tex_lines[w[0]].replace("\t", " "))
            if output:
                if use_color:
                    writer("    %s\033[33m%s\033[0m" % (" " * w[2][0], "^" * (w[2][1] - w[2][0])))
                else:
                    writer("    %s%s" % (" " * w[2][0], "^" * (w[2][1] - w[2][0])))
    return warnings


CATEGORY_GENERAL = 1
CATEGORY_TYPOGRAPHY = 2
CATEGORY_VISUAL = 4
CATEGORY_STYLE = 8
CATEGORY_REFERENCE = 16

checks = [
    (check_prefix,                      CATEGORY_REFERENCE,  "prefix", actions.fix_prefix),
    (check_mathmode_subscripts,         CATEGORY_TYPOGRAPHY, "mathmode-subscripts", actions.fix_mathmode_subscripts_in_file),
    (check_glossary_refs,               CATEGORY_REFERENCE,  "glossary-refs", actions.fix_glossary_refs_in_file),
    (check_space_before_cite,           CATEGORY_TYPOGRAPHY, "cite-space", actions.fix_space_before_cite_in_file),
    (check_math_text_mix,               CATEGORY_TYPOGRAPHY, "math-text-mix"),
    (check_units,                       CATEGORY_TYPOGRAPHY, "units"),
    (check_figure_alignment,            CATEGORY_STYLE,      "figure-alignment"),
    (check_table_alignment,             CATEGORY_STYLE,      "table-alignment"),
    (check_listing_alignment,           CATEGORY_STYLE,      "listing-alignment"),
    (check_float_placement,             CATEGORY_STYLE,      "float-placement"),
    (check_figure_has_label,            CATEGORY_REFERENCE,  "figure-label"),
    (check_table_has_label,             CATEGORY_REFERENCE,  "table-label"),
    (check_listing_has_label,           CATEGORY_REFERENCE,  "listing-label"),
    (check_figure_has_caption,          CATEGORY_STYLE,      "figure-caption"),
    (check_table_has_caption,           CATEGORY_STYLE,      "table-caption"),
    (check_listing_has_caption,         CATEGORY_STYLE,      "listing-caption"),
    (check_no_resizebox_for_tables,     CATEGORY_STYLE,      "resize-table"),
    (check_weird_units,                 CATEGORY_STYLE,      "dimensions"),
    (check_figure_caption_label_order,  CATEGORY_REFERENCE,  "figure-caption-order"),
    (check_table_caption_label_order,   CATEGORY_REFERENCE,  "table-caption-order"),
    (check_listing_caption_label_order, CATEGORY_REFERENCE,  "listing-caption-order"),
    (check_todos,                       CATEGORY_GENERAL,    "todo"),
    (check_notes,                       CATEGORY_GENERAL,    "note"),
    (check_unused_macro,                CATEGORY_GENERAL,    "unused-macro", actions.fix_unused_macros_in_file),
    (check_required_sections,           CATEGORY_GENERAL,    "required-sections"),
    (check_math_numbers,                CATEGORY_TYPOGRAPHY, "math-numbers"),
    (check_equation_symbols_defined,    CATEGORY_REFERENCE,  "symbol-mention"),
    (check_math_glossary_reference_coverage, CATEGORY_REFERENCE, "math-gls-coverage"),
    (check_large_numbers_without_si,    CATEGORY_TYPOGRAPHY, "si"),
    (check_listing_in_correct_float,    CATEGORY_REFERENCE,  "listing-float"),
    (check_tabular_in_correct_float,    CATEGORY_REFERENCE,  "tabular-float"),
    (check_tikz_in_correct_float,       CATEGORY_REFERENCE,  "tikz-float"),
    (check_comment_has_space,           CATEGORY_TYPOGRAPHY, "comment-space", actions.fix_comment_has_space_in_file),
    (check_percent_without_siunix,      CATEGORY_TYPOGRAPHY, "percentage"),
    (check_short_form,                  CATEGORY_GENERAL,    "short-form"),
    (check_labels_referenced,           CATEGORY_REFERENCE,  "label-referenced"),
    (check_section_capitalization,      CATEGORY_VISUAL,     "capitalization", actions.fix_section_capitalization_in_file),
    (check_quotation,                   CATEGORY_TYPOGRAPHY, "quotes", actions.fix_quotation_marks_in_file),
    (check_hline_in_table,              CATEGORY_VISUAL,     "hline"),
    (check_duplicate_word,              CATEGORY_TYPOGRAPHY, "duplicate-word", actions.fix_duplicate_words_in_file),
    (check_space_before_punctuation,    CATEGORY_TYPOGRAPHY, "punctuation-space", actions.fix_space_before_punctuation_in_file),
    (check_headers_without_text,        CATEGORY_VISUAL,     "two-header"),
    (check_one_sentence_paragraphs,     CATEGORY_VISUAL,     "single-sentence"),
    (check_paragraph_length,            CATEGORY_VISUAL,     "paragraph-length"),
    (check_multiple_sentences_per_line, CATEGORY_GENERAL,    "multiple-sentences"),
    (check_unbalanced_brackets,         CATEGORY_TYPOGRAPHY, "unbalanced-brackets"),
    (check_and_or,                      CATEGORY_TYPOGRAPHY, "and-or"),
    (check_ellipsis,                    CATEGORY_TYPOGRAPHY, "ellipsis"),
    (check_etc,                         CATEGORY_STYLE,      "etc"),
    (check_punctuation_end_of_line,     CATEGORY_TYPOGRAPHY, "punctuation"),
    (check_footnote,                    CATEGORY_TYPOGRAPHY, "footnote"),
    (check_table_vertical_lines,        CATEGORY_VISUAL,     "vline"),
    (check_section_balance,             CATEGORY_VISUAL,     "section-balance"),
    (check_section_length,              CATEGORY_VISUAL,     "section-length"),
    (check_table_top_caption,           CATEGORY_STYLE,      "table-top-caption"),
    (check_will,                        CATEGORY_GENERAL,    "will"),
    (check_subsection_count,            CATEGORY_VISUAL,     "single-subsection"),
    (check_subsection_count_sanity,     CATEGORY_VISUAL,     "subsection-count"),
    (check_mixed_compact_and_item,      CATEGORY_VISUAL,     "mixed-compact"),
    (check_center_in_float,             CATEGORY_VISUAL,     "float-center"),
    (check_appendix,                    CATEGORY_STYLE,      "appendix"),
    (check_eqnarray,                    CATEGORY_VISUAL,     "eqnarray"),
    (check_deprecated,                  CATEGORY_STYLE,      "deprecated"),
    (check_package_conflict,            CATEGORY_STYLE,      "package-conflict"),
    (check_acm_pc,                      CATEGORY_STYLE,      "inclusion"),
    (check_cite_noun,                   CATEGORY_STYLE,      "cite-noun"),
    (check_cite_duplicate,              CATEGORY_REFERENCE,  "cite-duplicate", actions.fix_cite_duplicate_in_file),
    (check_conjunction_start,           CATEGORY_STYLE,      "conjunction-start"),
    (check_brackets_space,              CATEGORY_TYPOGRAPHY, "bracket-spacing", actions.fix_brackets_space_in_file),
    (check_acronym_capitalization,      CATEGORY_TYPOGRAPHY, "acronym-capitalization"),
    (check_numeral,                     CATEGORY_GENERAL,    "numeral", actions.fix_numerals_in_file),
    (check_multicite,                   CATEGORY_STYLE,      "multiple-cites", actions.fix_multicite_in_file),
    (check_emptycite,                   CATEGORY_REFERENCE,  "cite-empty"),
    (check_colors,                      CATEGORY_VISUAL,     "colors"),
    (check_inconsistent_word_style,     CATEGORY_TYPOGRAPHY, "inconsistent-textstyle"),
    (check_missing_word_style,          CATEGORY_TYPOGRAPHY, "missing-textstyle", actions.fix_missing_word_style_in_file)
]

category_switches = [
    ("all",        CATEGORY_GENERAL | CATEGORY_REFERENCE | CATEGORY_STYLE | CATEGORY_TYPOGRAPHY | CATEGORY_VISUAL),
    ("general",    CATEGORY_GENERAL),
    ("reference",  CATEGORY_REFERENCE),
    ("style",      CATEGORY_STYLE),
    ("typography", CATEGORY_TYPOGRAPHY),
    ("visual",     CATEGORY_VISUAL)
]


def switch_exists(s):
    switches = [x[0] for x in category_switches] + [x[2] for x in checks]
    return s in switches


def expand_switch_to_checks(s):
    full_cat = [x[0] for x in category_switches]
    if s in full_cat:
        idx = full_cat.index(s)
        cat = category_switches[idx][1]
        return [check[2] for check in checks if check[1] & cat]
    return [s]


def add_categories(cat, new_cat):
    if type(new_cat) is str:
        full_cat = [x[0] for x in category_switches]
        if new_cat in full_cat:
            # full category, add everythingt that is not already there
            idx = full_cat.index(new_cat)
            new_cat = category_switches[idx][1]
        else:
            cat.add(new_cat)
    if type(new_cat) is int:
        for cats in checks:
            if new_cat & cats[1]:
                cat.add(cats[2])

        
def remove_categories(cat, rem_cat):
    if type(rem_cat) is str:
        full_cat = [x[0] for x in category_switches]
        if rem_cat in full_cat:
            # full category, add everythingt that is not already there
            idx = full_cat.index(rem_cat)
            rem_cat = category_switches[idx][1]
        else:
            if rem_cat in cat:
                cat.remove(rem_cat)
    if type(rem_cat) is int:
        for cats in checks:
            if (rem_cat & cats[1]) and cats[2] in cat:
                cat.remove(cats[2])


