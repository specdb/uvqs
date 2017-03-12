#!/usr/bin/env python
"""
Run a build of the DB
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import pdb

try:  # Python 3
    ustr = unicode
except NameError:
    ustr = str

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(
        description='Build the uvqs DB')
    parser.add_argument("-v", "--version", help="DB version to generate")
    parser.add_argument("-t", "--test", default=False, action='store_true', help="Test?")
    parser.add_argument("-c", "--clobber", default=False, action='store_true', help="Clobber existing file?")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args=None):
    """ Run
    Parameters
    ----------
    args

    Returns
    -------

    """
    from uvqs import build_db
    import h5py

    # Grab arguments
    pargs = parser(options=args)

    # Run
    if pargs.version is None:
        print("Building v01 of the uvqs DB")
        build_db.ver01(test=pargs.test, clobber=pargs.clobber)
    elif pargs.version == 'v01':
        print("Building v01 of the uvqs DB")
        build_db.ver01(test=pargs.test, clobber=pargs.clobber)
    elif pargs.version == 'v02':
        pdb.set_trace()
        print("Building v02 of the uvqs DB")
        #build_db.ver02(test=pargs.test, clobber=pargs.clobber)
    else:
        raise IOError("Bad version number")

if __name__ == '__main__':
    main()
