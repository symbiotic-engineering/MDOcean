from sphinx.ext.autosummary import Autosummary, generate as ag

def setup(app):
    # app.add_directive("autosummary", MatlabAutosummary, override=True)
    # app.connect("autosummary-context", merge_contexts)

    # def patch_generate(app):
    #     user_context = app.config.autosummary_context

    #     original_generate = ag.generate.generate_autosummary_docs

    #     def patched_generate_autosummary_docs(sources, *args, **kwargs):
    #         for source in sources:
    #             name = Path(source).stem
    #             # Patch get_context just for this source file
    #             ag.generate.get_context = lambda *a, **k: {
    #                 **old_get_context(*a, **k),
    #                 **user_context.get(name, {})
    #             }
    #             original_generate([source], *args, **kwargs)

    #         # Restore original get_context after all sources processed
    #         ag.generate.get_context = old_get_context

    #     ag.generate.generate_autosummary_docs = patched_generate_autosummary_docs

    # app.connect("builder-inited", patch_generate)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }

def merge_contexts(app, name, obj, parent, modname, full_name, context):

    # name is the base name of the page (e.g. 'simulation')
    custom = user_context.get(name, {})

    # Merge your custom keys into the default context
    for key, val in custom.items():
        if key in context:
            # Extend without duplicates
            context[key] = list(dict.fromkeys(context[key] + val))
        else:
            context[key] = val

    return context


class MatlabAutosummary(Autosummary):
    def get_table(self, *args, **kwargs):
        # No table, to avoid py:obj ref warnings
        return []

    # def get_items(self, names):
    #     env = self.env
    #     mat_data = env.domains.get('mat').data
    #     known_modules = set(mat_data.get('modules', {}).keys())

    #     collected = set()

    #     def collect_recursive(mod_name):
    #         if mod_name.startswith('__'):
    #             return
    #         collected.add(mod_name)
    #         prefix = mod_name + '.'
    #         for mod in known_modules:
    #             if mod.startswith(prefix) and mod not in collected:
    #                 print(f"recursing on submodule: {mod}")
    #                 collect_recursive(mod)

    #     for name in names:
    #         if not name.startswith('__'):
    #             if name in known_modules:
    #                 collect_recursive(name)
    #             else:
    #                 collected.add(name)

    #     results = []
    #     for mod in sorted(collected):
    #         if not mod.split('.')[-1].startswith('__'):
    #         # (name, signature, summary, real_name)
    #             results.append((mod, '', '', mod))

    #     print("Collected modules:", results)
    #     return results
    
def patched_get_context(self, fullname, **kwargs):
    # Called to prepare template context for this module stub page

    ctx = super().get_context(fullname, **kwargs)
    env = self.env
    mat_data = env.domains['mat'].data

    known_modules = mat_data.get('modules', {})
    known_objects = mat_data.get('objects', {})

    # Submodules are modules prefixed by fullname + '.'
    submodules = sorted(
        mod for mod in known_modules
        if mod.startswith(fullname + '.') and mod != fullname and not mod.split('.')[-1].startswith('__')
    )

    # Members are objects whose modname == fullname
    members = sorted(
        objname for objname, (modname, kind, _, _) in known_objects.items()
        if modname == fullname
        and not objname.startswith('__')
    )

    print("submodules: ", submodules)
    print("members: ", members)

    ctx['modules'] = [
        m for m in submodules
        if not m.startswith('__') and not m.endswith('__')
    ]

    ctx['members'] = [
        m for m in members
        if not m.startswith('__') and not m.endswith('__')
    ]

    return ctx
