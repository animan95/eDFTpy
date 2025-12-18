
"""
Module to expose more detailed version info for the installed `dftpy`
"""
version = "0.0.1dev0+git20251204.355d2b9"
__version__ = version
full_version = version

git_revision = "355d2b90acd425a0170f51de982342cb1f2db22f"
release = 'dev' not in version and '+' not in version
short_version = version.split("+")[0]
