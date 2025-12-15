# from https://github.com/sphinx-contrib/matlabdomain/issues/149#issuecomment-1435385185

#!/usr/bin/env python
"""
Replicate behavior of sphinx-apidoc.

Run with --help to see options.
"""
import argparse
import sys
from pathlib import Path
import collections

import logging
_log = logging.getLogger(sys.argv[0])


def parse_args():
    """Parse the command line arguments using argparse."""
    parser = argparse.ArgumentParser(
        description="Command Line Interface to generate matlab apidoc",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg = parser.add_argument
    arg('target')
    arg('--verbose', '-v', action='count', default=1,
        help='-v for info -vv for info')
    arg('--force', action='store_true', default=False,
        help='overwrite files')
    arg('--no-toc', action='store_true', default=False,
        help='Skip TOC file')
    arg('--output', '-o', default='./',
        help='Output folder')

    args = parser.parse_args()

    args.log_verbose = 40 - (10 * args.verbose) if args.verbose > 0 else 0
    return args


def _sanitize(s):
    return s.replace('_', r'\_')


def main(args):
    """
    Replicate behavior of sphinx-apidoc.

    :param  args:  arguments parsed from command line options (see :py:func:`parse_args`)
    """
    logging.basicConfig(level=args.log_verbose, format='%(levelname)s@%(filename)s:%(funcName)s:%(lineno)d:\n\t%(message)s')
    _log.info(f'CLI arguments: {vars(args)}')

    target_path = Path(args.target)
    if not target_path.is_dir():
        _log.error(f'{target_path} is not a directory')
        return -1

    for candidate_pkg in target_path.iterdir():
        if not candidate_pkg.is_dir() or candidate_pkg.name == 'build':
            continue

        toolbox_doc_strings = _autodoc_toolbox(candidate_pkg)
        if len(toolbox_doc_strings) == 0:
            continue

        toolbox_name = candidate_pkg.stem
        ouput_path = Path(args.output)
        toolbox_output_path = ouput_path / toolbox_name
        toolbox_output_path.mkdir(exist_ok=True, parents=True)
        if not toolbox_output_path.is_dir():
            _log.error(f'Could not create output directory {toolbox_output_path}')
            return -1

        with open(ouput_path / f'{toolbox_name}.rst', 'wt') as f:
            f.write(f'''
{_sanitize(toolbox_name)} Toolbox
{'#'*len(toolbox_name+' Toolbox')}
This is a Matlab toolbox which provides the following packages/namespaces:

.. toctree::
   :maxdepth: 2
   :glob:

   {toolbox_name}/*
''')
        for pkg_name, pkg_doc_string in toolbox_doc_strings.items():
            if pkg_name is not None:
                api_filename = toolbox_output_path / (pkg_name + '.rst')
                if api_filename.is_file() and not args.force:
                    _log.error(f'Output file {api_filename} already exists, use --force to overwrite')
                    return - 1

                with open(api_filename, 'wt') as f:
                    f.write(pkg_doc_string)

    if not args.no_toc:
        raise NotImplementedError('TOC not implemented')

    return 0


def _autodoc_toolbox(path):
    """
    Autogenerates sphinx rst files for a matlab toolbox.

    :param      path:  The path to the toolbox
    """
    toolbox_doc_strings = collections.defaultdict(str)
    toolbox_name = path.stem

    pkg_string = '''
{pkg_sanitized} package
{equal_signs}

.. mat:automodule:: {pkg}
   :members:
   :undoc-members:
   :show-inheritance:

'''

    subpkg_string = '''
{subpkg_sanitized} {subtype}
{dash_signs}
.. mat:automodule:: {toolbox_name}.{subpkg}
   :members:
   :undoc-members:
   :show-inheritance:
'''

    n_sub_packages = 0
    for p in path.rglob('*'):
        if p.is_dir() and p.stem.startswith('+'):
            subtype = 'package'
        elif p.is_dir() and p.stem.startswith('@'):
            subtype = 'class'
        else:
            # TODO: handle scripts ?
            continue

        _log.info(p)

        subpkg_path = str(p.relative_to(path)).split('/')
        subpkg = '.'.join(subpkg_path)

        if subtype == 'package' and len(subpkg_path) == 1:
            _log.info(f'package = {subpkg_path}')
            parent_key = f'{toolbox_name}.{subpkg}'
            n_sub_packages += 1
            toolbox_doc_strings[parent_key] = pkg_string.format(
                pkg=parent_key,
                pkg_sanitized=_sanitize(subpkg),
                equal_signs='=' * len(_sanitize(subpkg)) + '=' * 7)
        else:
            _log.info(f'package = {subpkg_path}')
            parent_key = '.'.join(subpkg_path[:-1])
            parent_key = f'{toolbox_name}.{parent_key}'
            _log.info(f'parent_key = {parent_key} -> {toolbox_doc_strings[parent_key]}')

            toolbox_doc_strings[parent_key] += subpkg_string.format(
                toolbox_name=toolbox_name,
                subpkg=subpkg, subtype=subtype,
                subpkg_sanitized=_sanitize(subpkg),
                dash_signs='-' * (len(subpkg) + len(subtype) + 1)
            )

    if n_sub_packages > 0:
        return toolbox_doc_strings
    else:
        return {}


sys.exit(main(parse_args()))
