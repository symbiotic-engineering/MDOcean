from sphinx.ext.autosummary import Autosummary, generate as ag, autosummary_table
from pathlib import Path
from sphinxcontrib import mat_types
import sphinx
sphinx_version = int(sphinx.__version__.split(".")[0])

def generate_matlab_modulelist(app,config):
    srcdir = Path(app.srcdir)
    toolbox_root = Path(config.matlab_src_dir)
    if toolbox_root is None:
        raise RuntimeError("matlab_src_dir is not set in conf.py")
    rst_file = srcdir / "api_generate.rst"

    entries = []

    for path in toolbox_root.glob("*"):
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
    app.add_directive("autosummary", MatlabAutosummary, override=True)
    app.connect("config-inited", generate_matlab_modulelist)

    def patch_generate(app, config):
        user_context = config.autosummary_context

        original_generate = ag.generate_autosummary_docs
        original_generate_content = ag.generate_autosummary_content

        def patched_generate_autosummary_docs(sources, *args, **kwargs):
            results = []
            for source in sources:
            
                def patched_generate_autosummary_content(*args, **kwargs):
                    name = args[0]
                    ctx = user_context.get(name, {})
                    if sphinx_version == 8:
                        idx = 7
                    elif sphinx_version == 7:
                        idx = 8
                    else:
                        raise RuntimeError(f"Unsupported Sphinx version: {sphinx_version}")
                    
                    old_ctx = args[idx] # hardcode context (7th argument)
                    assert type(old_ctx) is dict, f"Expected dict for context, got {type(old_ctx)}"
                    new_args = args[:idx] + (ctx,) + args[idx+1:]
                    print(f"Generating autosummary content for {name} with context {ctx}")
                    return original_generate_content(*new_args, **kwargs)

                ag.generate_autosummary_content = patched_generate_autosummary_content

                files = original_generate([source], *args, **kwargs)
                if files is None:
                    files = []  # avoid NoneType
                results.extend(files)
            return results

        ag.generate_autosummary_docs = patched_generate_autosummary_docs

    app.connect("config-inited", patch_generate)

    return {
        'version': '0.1',
        'parallel_read_safe': False,
        'parallel_write_safe': False,
    }


class MatlabAutosummary(Autosummary):
    def get_table(self, *args, **kwargs):
        from docutils import nodes
        from sphinx.addnodes import pending_xref
        
        table = super().get_table(*args, **kwargs)
        
        # table is a list of nodes, traverse each one
        for table_node in table:
            # Traverse the entire node tree to find pending_xref and literal nodes
            for node in table_node.traverse():
                # Change pending_xref attributes from py to mat
                if isinstance(node, pending_xref):
                    if node.get('refdomain') == 'py':
                        node['refdomain'] = 'mat'
                    if node.get('reftype') == 'obj':
                        node['reftype'] = 'obj'  # Keep as obj
                    # Update py:class and py:module attributes
                    if node.get('py:class'):
                        node['mat:class'] = node.pop('py:class')
                    if node.get('py:module'):
                        node['mat:module'] = node.pop('py:module')
                
                # Change literal class attributes from "xref py py-obj" to "xref mat mat-obj"
                elif isinstance(node, nodes.literal):
                    classes = node.get('classes', [])
                    if 'py' in classes:
                        # Replace 'py' with 'mat' and 'py-obj' with 'mat-obj'
                        new_classes = []
                        for cls in classes:
                            if cls == 'py':
                                new_classes.append('mat')
                            elif cls == 'py-obj':
                                new_classes.append('mat-obj')
                            else:
                                new_classes.append(cls)
                        node['classes'] = new_classes
        return table

    def get_items(self, names):
        print("MatlabAutosummary get_items called with names:", names)
        env = self.env
        mat_data = env.domains.get('mat').data
        known_modules = set(mat_data.get('modules', {}).keys())
        known_classes = set(mat_data.get('classes', {}).keys())

        entities_table = mat_types.analyze(self.env.app)

        print("Known modules:", known_modules)
        print("Known classes:", known_classes)

        collected = set()

        def collect_recursive(mod_name):
            if mod_name.startswith('__'):
                return
            collected.add(mod_name)
            prefix = mod_name + '.'
            for mod in known_modules:
                if mod.startswith(prefix) and mod not in collected:
                    print(f"recursing on submodule: {mod}")
                    collect_recursive(mod)
            for cls in known_classes:
                if cls.startswith(prefix) and cls not in collected:
                    print(f"collecting class: {cls}")
                    collected.add(cls)

        for name in names:
            if not name.startswith('__'):
                if name in known_modules or name in known_classes:
                    collect_recursive(name)
                else:
                    collected.add(name)

        results = []
        for mod in sorted(collected):
            if not mod.split('.')[-1].startswith('__'):
            # (name, signature, summary, real_name)
                results.append((mod, '', '', mod))

        print("Collected modules:", results)
        return results
