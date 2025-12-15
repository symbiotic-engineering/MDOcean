# How to build docs

Create/activate environment:
```
conda create -n sphinxenv3 -c conda-forge pip sphinxcontrib-matlabdomain myst-parser
pip install furo==2021.11.16
conda activate sphinxenv2
```
Build:
```
rm -rf docs/generated
sphinx-build -a -W --keep-going -b html docs docs/_build/html
```

# Docs debugging log

On 7/29/25 and 12/14/25, I was trying a number of ways to get nested toctrees for matlab.
1. Use the python apidoc extension directly. This does not support the matlab domain and I did not pursue it further.
2. Use the python autosummary extension directly. This also does not support the matlab domain so results in `py:obj ref not found` warnings, which I fixed by overriding the table to return none, but this had still the issue of containing only the high level modules and functions, not any submodules.
3. Custom matlab autosummary extension: see `_ext/matlab_autosummary.py` to detect the matlab submodules. This is what ultimately worked. The solution was to modify the context within the python autosummary (overload the `generate_autosummary_content` and `generate_autosummary_docs` functions), combined with a modified table and modifying the `_templates/autosummary/module.rst` template file).
4. matlab apidoc extension suggested by user zouhairm of matlabdomain: see `_ext/matlab_apidoc.py` (taken from [this link](https://github.com/sphinx-contrib/matlabdomain/issues/149#issuecomment-1435385185)). This does not actually produce table of contents because it uses notoc, so does not solve the problem.
5. matlab autosummary extension suggested by user changlichun of matlabdomain at [this PR](https://github.com/sphinx-contrib/matlabdomain/pull/211). I tried using this fork with Sphinx 8.2.3 and got an error about the `split_full_qualified_name` function not existing. I downgraded to sphinx 7.1.2 and it fixed it. Then I got error
    ```
    sphinx.errors.ExtensionError: Handler <function analyze at 0x7f31a56e2d40> for event 'builder-inited' threw an exception (exception: 'NoneType' object has no attribute 'safe_getmembers')
    ```
    I fixed that by using the python debugger to notice that it was hitting the following if statement: `if package.startswith("_") or package.startswith("."): return None`. I fixed this by using the absolute path to the module, which got the script to run. But, it has the problem I hoped it would solve: it doesn't recognize matlab modules (with or without +) to show them in the toc and expand them recursively.
6. Manually generate the rst files using a python script, without using any extension. See `_ext/matlab_walk.py`. This currently is only successful at producing .rst files for only the top-level modules, not submodules, although the code is intended to include submodules.