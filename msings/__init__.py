from os.path import join, dirname

_data = join(dirname(__file__), 'data')

try:
    with open(join(_data, 'sha')) as s, open(join(_data, 'ver')) as v:
        sha = s.read().strip()
        ver = int(v.read())
except Exception, e:
    __version__ = ''
else:
    __version__ = '%04i.%s' % (ver, sha)

