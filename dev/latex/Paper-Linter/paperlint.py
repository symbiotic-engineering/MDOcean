#!/usr/bin/env python3
import os
import sys

LATEX_DEV_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
if LATEX_DEV_DIR not in sys.path:
    sys.path.append(LATEX_DEV_DIR)

import checks
import glossary_symbols as gs
import options_settings as opts
import tex_structure as ts

GLOSSARY_FILE = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "pubs", "shared", "symbol-glossary-shared.tex"))

output_handle = sys.stdout
use_color = True


def write_output(text="", end="\n"):
    print(text, end=end, file=output_handle)


def main():
    global output_handle, use_color

    nr_warnings = 0
    nr_suppressed = 0
    autofix_counts = {}
    autofix_files = {}
    checks.math_text_mix_strings = set()
    all_equation_symbols = set()

    options = opts.parse_cli_args(
        sys.argv,
        checks.switch_exists,
        checks.add_categories,
        checks.remove_categories,
    )

    checks.GLOSSARY_FILE = GLOSSARY_FILE
    checks.ACRONYM_GLOSSARY_FILE = options["acronym_glossary_seed_file"]
    checks.PREFIX_CHECK_CONFIGS = options["check_params"]
    if options["replace_glossary_refs"]:
        options["check_states"]["glossary-refs"] = 2

    if options["symbol_glossary_seed_file"] is not None:
        try:
            with open(options["symbol_glossary_seed_file"]) as handle:
                checks.actions.WORD_STYLE_REFERENCE_TEXT = handle.read()
        except Exception:
            checks.actions.WORD_STYLE_REFERENCE_TEXT = ""

    if options["output_file"] is not None:
        output_handle = open(options["output_file"], "w")
        use_color = False
    else:
        output_handle = sys.stdout
        use_color = output_handle.isatty()

    tex_files = ts.collect_tex_files(options["root_path"], options["ignored_files"])

    missing_word_style_check = None
    missing_word_style_prepass = options["check_states"].get("missing-textstyle", 1) == 2
    if missing_word_style_prepass:
        for c in checks.checks:
            if c[2] == "missing-textstyle":
                missing_word_style_check = c
                break

    try:
        for file in tex_files:
            try:
                ts.next_file(file)
            except RuntimeError:
                print("Could not open '%s'" % sys.argv[1])
                sys.exit(1)

            if use_color:
                write_output("Inspecting file \033[94m'%s'\033[0m" % file)
            else:
                write_output("Inspecting file '%s'" % file)

            ts.preprocess()

            if missing_word_style_prepass and missing_word_style_check is not None and len(missing_word_style_check) > 3 and missing_word_style_check[3] is not None:
                missing_word_style_check[3](file)
                try:
                    ts.next_file(file)
                except RuntimeError:
                    print("Could not open '%s'" % sys.argv[1])
                    sys.exit(1)
                ts.preprocess()

            warnings = []
            suppressed = []
            for c in checks.checks:
                state = options["check_states"].get(c[2], 1)
                config = options["check_params"].get(c[2], {})
                if state == 0:
                    add_warn = []
                else:
                    add_warn = c[0]()
                    pre_fix_warning_count = len(add_warn)
                    if c[2] == "missing-textstyle" and missing_word_style_prepass:
                        pass
                    elif state == 2 and len(c) > 3 and c[3] is not None and len(add_warn) > 0:
                        if c[2] == "glossary-refs":
                            c[3](file, GLOSSARY_FILE, options["acronym_glossary_seed_file"])
                        elif c[2] == "prefix":
                            c[3](config.get("path"), config.get("prefix"), recursive=config.get("recursive", False), source_prefix=config.get("source_prefix"))
                        else:
                            c[3](file)

                        try:
                            ts.next_file(file)
                        except RuntimeError:
                            print("Could not open '%s'" % sys.argv[1])
                            sys.exit(1)
                        ts.preprocess()

                        add_warn = c[0]()
                        fixed_count = pre_fix_warning_count - len(add_warn)
                        if fixed_count > 0:
                            autofix_counts[c[2]] = autofix_counts.get(c[2], 0) + fixed_count
                            if c[2] not in autofix_files:
                                autofix_files[c[2]] = set()
                            autofix_files[c[2]].add(file)
                if c[2] == "symbol-mention":
                    all_equation_symbols.update(checks.current_file_equation_symbols)
                if state > 0:
                    warnings += [(x, c[2]) for x in add_warn]
                else:
                    suppressed += [(x, c[2]) for x in add_warn]

            nr_warnings += checks.print_warnings(warnings, write_output, use_color=use_color)
            nr_suppressed += checks.print_warnings(suppressed, write_output, use_color=use_color, output=False)

        write_output("")
        write_output("%d warnings printed; %d suppressed warnings" % (nr_warnings, nr_suppressed))
        if autofix_counts:
            write_output("")
            write_output("Auto-fix summary:")
            for check_name in sorted(autofix_counts.keys()):
                write_output("- [%s] fixed %d issue(s) in %d file(s)" % (check_name, autofix_counts[check_name], len(autofix_files.get(check_name, set()))))
        if options["check_states"].get("math-text-mix", 0) > 0 and len(checks.math_text_mix_strings) > 0:
            write_output("")
            write_output("Unique text strings in math environments [math-text-mix]:")
            for text_string in sorted(checks.math_text_mix_strings):
                write_output("- %s" % text_string)
        if options["symbol_glossary_file"] is not None:
            used_glossary_labels = gs.collect_glossary_labels_from_files(tex_files)
            gs.write_symbol_glossary(
                options["symbol_glossary_file"],
                all_equation_symbols,
                options["symbol_glossary_seed_file"],
                used_labels=used_glossary_labels,
            )
            write_output("Wrote %d symbols to '%s'" % (len(all_equation_symbols), options["symbol_glossary_file"]))
    finally:
        if output_handle is not sys.stdout:
            output_handle.close()

    if options["exit_code"]:
        sys.exit(1 if nr_warnings > 0 else 0)


if __name__ == "__main__":
    main()
