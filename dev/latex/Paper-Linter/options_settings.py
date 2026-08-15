import sys
from pathlib import Path

import checks
import tex_structure as ts


def usage(program_name):
    print("%s <file.tex/path> [-x <excluded-switch1>] [-i <included-switch1>] [--ignore <file-or-name>] [--settings <settings-file>] [--params <params-file>] [--output <output-file>] [--symbol-glossary-output <output-file>] [--symbol-glossary-seed <seed-file>] [--replace-glossary-refs] [-i/x <switch n, evaluated in order of specification>] [--error]" % program_name)
    sys.exit(1)


def apply_settings_file(settings_file, check_states, check_params, switch_exists):
    try:
        with open(settings_file) as handle:
            lines = handle.readlines()
    except Exception:
        print("Could not open settings file '%s'" % settings_file)
        sys.exit(1)

    for line_no, raw_line in enumerate(lines, start=1):
        line = ts.strip_comment_from_line(raw_line).strip()
        if not line:
            continue

        parts = line.split()
        if len(parts) < 2:
            print("Malformed settings line %d in '%s': expected '0|1|2 check-name'" % (line_no, settings_file))
            sys.exit(1)

        state, check_name = parts[0], parts[1]
        if state not in ["0", "1", "2"]:
            print("Malformed settings line %d in '%s': expected 0, 1, or 2, got '%s'" % (line_no, settings_file, state))
            sys.exit(1)

        if not switch_exists(check_name):
            print("Unknown switch '%s' in settings file '%s' on line %d" % (check_name, settings_file, line_no))
            sys.exit(1)

        params = {}
        for token in parts[2:]:
            if "=" not in token:
                print("Malformed settings line %d in '%s': expected key=value parameters" % (line_no, settings_file))
                sys.exit(1)
            key, value = token.split("=", 1)
            if key == "recursive":
                value = value.lower() in ["1", "true", "yes", "on"]
            params[key] = value

        for name in checks.expand_switch_to_checks(check_name):
            check_states[name] = int(state)
            if params:
                check_params[name] = params


def apply_params_file(params_file, options):
    params_path = Path(params_file).resolve()

    def resolve_path(value):
        path = Path(value)
        if path.is_absolute():
            return str(path)
        return str(path.resolve())

    try:
        with open(params_path) as handle:
            lines = handle.readlines()
    except Exception:
        print("Could not open params file '%s'" % params_file)
        sys.exit(1)

    prefix_params = options["check_params"].setdefault("prefix", {})

    for line_no, raw_line in enumerate(lines, start=1):
        line = ts.strip_comment_from_line(raw_line).strip()
        if not line:
            continue

        if "=" not in line:
            print("Malformed params line %d in '%s': expected key=value" % (line_no, params_file))
            sys.exit(1)

        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()

        if key == "recursive":
            value = value.lower() in ["1", "true", "yes", "on"]

        if key in ["path", "prefix", "recursive", "source_prefix"]:
            prefix_params[key] = resolve_path(value) if key == "path" else value
        elif key in ["ignore", "ignored_file"]:
            options["ignored_files"].append(resolve_path(value))
        elif key in ["symbol_glossary_seed_file", "symbol_glossary_seed"]:
            options["symbol_glossary_seed_file"] = resolve_path(value)
        elif key in ["acronym_glossary_seed_file", "acronym_glossary_seed"]:
            options["acronym_glossary_seed_file"] = resolve_path(value)
        elif key in ["symbol_glossary_file", "symbol_glossary_output"]:
            options["symbol_glossary_file"] = resolve_path(value)
        elif key in ["output_file", "output"]:
            options["output_file"] = resolve_path(value)
        else:
            print("Unknown params key '%s' in params file '%s' on line %d" % (key, params_file, line_no))
            sys.exit(1)


def parse_cli_args(argv, switch_exists, add_categories, remove_categories):
    if len(argv) < 2:
        usage(argv[0])

    options = {
        "ignored_files": [],
        "output_file": None,
        "symbol_glossary_file": None,
        "symbol_glossary_seed_file": None,
        "acronym_glossary_seed_file": None,
        "replace_glossary_refs": False,
        "exit_code": False,
        "check_states": {},
        "check_params": {},
        "has_rules": False,
        "root_path": argv[1],
    }

    for check_name in [check[2] for check in checks.checks]:
        options["check_states"][check_name] = 1

    idx = 1
    while idx < len(argv):
        arg = argv[idx]
        if arg == "-x":
            if idx + 1 < len(argv):
                if switch_exists(argv[idx + 1]):
                    for name in checks.expand_switch_to_checks(argv[idx + 1]):
                        options["check_states"][name] = 0
                    idx += 1
                    options["has_rules"] = True
                else:
                    print("Unknown switch '%s'" % argv[idx + 1])
                    usage(argv[0])
            else:
                print("Missing switch after -x")
                usage(argv[0])

        if arg == "-i":
            if idx + 1 < len(argv):
                if switch_exists(argv[idx + 1]):
                    for name in checks.expand_switch_to_checks(argv[idx + 1]):
                        options["check_states"][name] = 1
                    idx += 1
                    options["has_rules"] = True
                else:
                    print("Unknown switch '%s'" % argv[idx + 1])
                    usage(argv[0])
            else:
                print("Missing switch after -i")
                usage(argv[0])

        if arg == "--settings":
            if idx + 1 < len(argv):
                apply_settings_file(
                    argv[idx + 1],
                    options["check_states"],
                    options["check_params"],
                    switch_exists,
                )
                idx += 1
                options["has_rules"] = True
            else:
                print("Missing file after --settings")
                usage(argv[0])

        if arg == "--params":
            if idx + 1 < len(argv):
                apply_params_file(argv[idx + 1], options)
                idx += 1
                options["has_rules"] = True
            else:
                print("Missing file after --params")
                usage(argv[0])

        if arg == "--ignore":
            if idx + 1 < len(argv):
                options["ignored_files"].append(argv[idx + 1])
                idx += 1
                options["has_rules"] = True
            else:
                print("Missing file after --ignore")
                usage(argv[0])

        if arg == "--output":
            if idx + 1 < len(argv):
                options["output_file"] = argv[idx + 1]
                idx += 1
            else:
                print("Missing file after --output")
                usage(argv[0])

        if arg == "--symbol-glossary-output":
            if idx + 1 < len(argv):
                options["symbol_glossary_file"] = argv[idx + 1]
                idx += 1
            else:
                print("Missing file after --symbol-glossary-output")
                usage(argv[0])

        if arg == "--symbol-glossary-seed":
            if idx + 1 < len(argv):
                options["symbol_glossary_seed_file"] = argv[idx + 1]
                idx += 1
            else:
                print("Missing file after --symbol-glossary-seed")
                usage(argv[0])

        if arg == "--replace-glossary-refs":
            options["replace_glossary_refs"] = True

        if arg == "--error":
            options["exit_code"] = True

        idx += 1

    return options
