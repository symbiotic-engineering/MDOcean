from sphinx.ext.autosummary import Autosummary

class NoTableAutosummary(Autosummary):
    # prevent table generation which causes ref warnings with matlab domain
    def get_table(self, *args, **kwargs):
        return []

def setup(app):
    app.add_directive("autosummary", NoTableAutosummary, override=True)
    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }

