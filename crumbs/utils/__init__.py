from sys import version_info

if version_info[0] < 3 or (version_info[0] == 3 and version_info[1] < 3):
    # until python 3.3 the standard file module has no support for
    # wrapping file object and required to open a new file
    # bz2file is a backport of the python 3.3 std library module
    try:
        from bz2file import BZ2File
    except ImportError:
        pass
else:
    from bz2 import BZ2File


def approx_equal(float_a, float_b, tolerance=0.01):
    'an aproach to compare two floating numbers'
    tol = tolerance
    return (abs(float_a - float_b) / max(abs(float_a), abs(float_b))) < tol
