from sphinx.ext.autosummary import Autosummary
from pathlib import Path

class NoTableAutosummary(Autosummary):
    # prevent table generation which causes ref warnings with matlab domain
    def get_table(self, *args, **kwargs):
        return []

def generate_matlab_modulelist(app,config):
    srcdir = Path(app.srcdir)
    toolbox_root = Path(config.matlab_src_dir)
    if toolbox_root is None:
        raise RuntimeError("matlab_src_dir is not set in conf.py")
    rst_file = srcdir / "api_generate.rst"

    entries = []

    for path in toolbox_root.rglob("*"):
        if not path.is_dir():
            continue
        if path.name.startswith((".", "+private", "@")):
            continue
        if not any(path.rglob("*.m")):
            continue
        if "generated" in path.parts:
            continue
        if "SAFE" in path.parts:
            continue

        dotted = ".".join(path.relative_to(toolbox_root).parts)
        entries.append(dotted)

    entries = sorted(set(entries))

    # Replace a marked block
    begin = ".. BEGIN MATLAB AUTO API"
    end = ".. END MATLAB AUTO API"

    text = rst_file.read_text(encoding="utf-8")

    block = "\n".join(f"    {e}" for e in entries)

    new_text = (
        text.split(begin)[0]
        + begin + "\n"
        + block + "\n"
        + end
        + text.split(end)[1]
    )

    rst_file.write_text(new_text, encoding="utf-8")
    print(f"Updated {rst_file}")

def setup(app):
    app.add_directive("autosummary", NoTableAutosummary, override=True)
    app.connect("config-inited", generate_matlab_modulelist)
    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }

