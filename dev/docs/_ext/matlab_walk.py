#!/usr/bin/env python3

import os
from pathlib import Path

RST_TEMPLATE = '''\
{modname}
{underline}

.. mat:automodule:: {modname}
   :members:
   :undoc-members:
   :show-inheritance:

{autosummary_block}
'''

TOCTREE_TEMPLATE = '''\
.. toctree::
   :maxdepth: 2

{entries}
'''

def matlab_name_from_path(path_parts):
    # Convert directory parts into MATLAB module names (dot-separated)
    return '.'.join(path_parts)

def class_name_from_path(path):
    return path.name[1:] if path.name.startswith('@') else path.name

def generate_rst_for_package(pkg_path, output_dir, rel_parts):
    modname = matlab_name_from_path(rel_parts)
    title = f"{modname} package"
    underline = '=' * len(title)

    children = []
    for entry in sorted(pkg_path.iterdir()):
        if entry.is_dir():
            children.append(entry.name.lstrip('@'))
        # Note: this script skips .m files at the package level for now

    entries = '\n'.join(f'   {modname}.{child}' for child in children)
    toctree = TOCTREE_TEMPLATE.format(entries=entries) if children else ''

    rst_text = RST_TEMPLATE.format(
        title=title, underline=underline,
        modname=modname, toctree=toctree
    )

    out_path = output_dir / f"{modname}.rst"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(rst_text)
    print(f"Generated: {out_path}")

def generate_rst_for_class(class_path, output_dir, rel_parts):
    pkg_name = matlab_name_from_path(rel_parts[:-1])
    class_name = class_name_from_path(class_path)
    full_name = f"{pkg_name}.{class_name}"
    title = f"{class_name} class"
    underline = '=' * len(title)

    rst_text = RST_TEMPLATE.format(
        title=title, underline=underline,
        modname=full_name, toctree=''
    )

    out_path = output_dir / f"{full_name}.rst"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(rst_text)
    print(f"Generated: {out_path}")

def scan_matlab_package(pkg_path, rel_parts):
    modname = matlab_name_from_path(rel_parts)
    submodules = []

    for entry in sorted(pkg_path.iterdir()):
        if entry.is_dir():
            if entry.name.startswith('@'):
                submodules.append(f"{modname}.{entry.name[1:]}")
            else:
                submodules.append(f"{modname}.{entry.name}")

    return submodules



def walk_matlab_toolbox(pkg_path, output_dir, rel_parts):
    modname = matlab_name_from_path(rel_parts)
    submodules = []

    for entry in sorted(pkg_path.iterdir()):
        if entry.is_dir():
            if entry.name.startswith('@'):
                submodules.append(f"{modname}.{entry.name[1:]}")
            else:
                submodules.append(f"{modname}.{entry.name}")
    
    autosummary_block = ''
    if submodules:
        autosummary_entries = '\n'.join(f'   {s}' for s in submodules)
        autosummary_block = f'''.. rubric:: Submodules

.. autosummary::
   :toctree: .
   :template: matmod.rst
   :recursive:

{autosummary_entries}
'''

    rst_text = RST_TEMPLATE.format(
        modname=modname,
        underline='=' * len(modname),
        autosummary_block=autosummary_block
    )

    out_path = output_dir / f"{modname}.rst"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(rst_text)
    print(f"Generated: {out_path}")

    # Recurse into non-class subdirectories
    for entry in sorted(pkg_path.iterdir()):
        if entry.is_dir() and not entry.name.startswith('@'):
            walk_matlab_toolbox(entry, output_dir, rel_parts + [entry.name])

    # Generate RST for @class directories
    for entry in sorted(pkg_path.iterdir()):
        if entry.is_dir() and entry.name.startswith('@'):
            generate_rst_for_class(entry, output_dir, rel_parts + [entry.name])


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("toolbox_dir", help="Path to folder containing your MATLAB packages")
    parser.add_argument("-o", "--output", default="./generated", help="Output folder for .rst files")
    args = parser.parse_args()

    toolbox_dir = Path(args.toolbox_dir).resolve()
    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    walk_matlab_toolbox(toolbox_dir, output_dir, [toolbox_dir.name])


if __name__ == "__main__":
    main()
