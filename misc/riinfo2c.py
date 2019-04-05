'''INCOMPLETE! This will read read the refractiveindex.info yaml files
and transforms the database into a C source.'''

import re
import os
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

# Right now, we can process only the 'tabulated nk' data
searchfor = '- type: tabulated nk'
searchfor = re.compile(searchfor)

ridatadir = "/u/46/necadam1/unix/repo/refractiveindex.info-database/database/data"

nktables = dict()

def find_files_by_pattern (pattern, dir):
  r = re.compile(pattern)
  for parent, dnames, fnames in os.walk(ridatadir):
    for fname in fnames:
        filename = os.path.join(parent, fname)
        if os.path.isfile(filename):
            with open(filename) as f:
                text = f.read()
                if r.search(text):
                    yield (



