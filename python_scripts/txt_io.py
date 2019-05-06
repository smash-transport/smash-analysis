from numpy import vstack, savetxt, genfromtxt
from collections import OrderedDict

def save_table(data, f, comments=None):
    """Make a text table from given dictionary and store it in a file.

    The keys will be used as column names. The values should be iterables of
    equal length.

    Optionally, comments can passed that will be prefixed with '#' and appended
    to the end of the file.
    """
    if type(f) is str:
        f = open(f, 'w')
        got_path = True
    else:
        got_path = False
    header = ' '.join(data.keys())
    array = vstack(data.values()).T
    f.write(header)
    f.write('\n')
    savetxt(f, array)
    if comments:
        f.write('# ')
        f.write(comments)
    if got_path:
        f.close()
    return array

def load_table(f):
    """Read a text table from a file and return it as a dictionary.

    The keys are the column names.
    """
    if type(f) is str:
        path = f
        f = open(f, 'r')
        got_path = True
    else:
        got_path = False
    header = '#'
    while header.startswith('#'):
        header = f.readline().rstrip('\n')
    if 'detailed_balance' in path:
        # Unicode strings in first column. Data type needs to be specified.
        array = genfromtxt(f, missing_values='-', dtype='unicode',
                           encoding='utf-8')
    else:
        array = genfromtxt(f, missing_values='-')
    if got_path:
        f.close()
    return OrderedDict(zip(header.split(), array.T))

def test_save_table():
    from StringIO import StringIO

    # save
    d = OrderedDict([('col1', range(5)), ('col2', range(5, 10))])
    f = StringIO()
    save_table(d, f)
    s = f.getvalue()
    f.close()

    # read
    f = StringIO(s)
    for i1, i2 in zip(load_table(f).items(), d.items()):
        assert i1[0] == i2[0]
        assert (i1[1] == i1[1]).all()

if __name__ == "__main__":
    test_save_table()
