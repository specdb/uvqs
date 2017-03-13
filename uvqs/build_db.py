""" Module to build the hdf5 database file for IGMspec
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, warnings

import h5py
import json
import datetime
import pdb
from collections import OrderedDict

from astropy.table import Table, vstack, Column
from astropy import units as u

from linetools import utils as ltu

from specdb import defs
from specdb.build import utils as sdbbu

import uvqs
from uvqs import fuv

try:
    bstr = bytes
except NameError:  # For Python 2
    bstr = str

def ver01(test=False, skip_copy=False, publisher='J.X. Prochaska', clobber=False):
    """ Build version 2.X

    Reads previous datasets from v1.X

    Parameters
    ----------
    test : bool, optional
      Run test only
    skip_copy : bool, optional
      Skip copying the data from v01

    Returns
    -------
    """
    version = 'v01'
    # Read v01
    outfil = uvqs.__path__[0]+'/../DB/UVQS_DB_{:s}.hdf5'.format(version)

    # Main DB Table
    idkey = 'UVQS_ID'
    maindb, tkeys = sdbbu.start_maindb(idkey)

    # Clobber?
    if not chk_clobber(outfil, clobber=clobber):
        return
    # Begin
    hdf = h5py.File(outfil,'w')

    group_dict = {}
    new_groups = get_build_groups(version)
    meta_only = False
    # Loop over the new groups
    for gname in new_groups:
        print("Working on group: {:s}".format(gname))
        # Meta
        meta = new_groups[gname].grab_meta()
        # Survey flag
        flag_g = sdbbu.add_to_group_dict(gname, group_dict, skip_for_debug=True)
        # IDs
        maindb = sdbbu.add_ids(maindb, meta, flag_g, tkeys, idkey, mtch_toler=5*u.arcsec,
                               pair_sep=4.9*u.arcsec, first=(flag_g==1))#, close_pairs=(gname in pair_groups))
        # Spectra
        if not meta_only:
            new_groups[gname].hdf5_adddata(hdf, gname, meta)
            new_groups[gname].add_ssa(hdf, gname)

    # Check for duplicates -- There is 1 pair in SDSS (i.e. 2 duplicates)
    if not sdbbu.chk_for_duplicates(maindb, dup_lim=2):
        raise ValueError("Failed duplicates")

    # Check stacking
    if not sdbbu.chk_vstack(hdf):
        print("Meta data will not stack using specdb.utils.clean_vstack")
        print("Proceed to write at your own risk..")
        pdb.set_trace()

    # Finish
    zpri = defs.z_priority()
    sdbbu.write_hdf(hdf, str('uvqs'), maindb, zpri,
                    group_dict, version, Publisher=str(publisher))

    print("Wrote {:s} DB file".format(outfil))
    print("Update DB info in specdb.defs.dbase_info !!")


def chk_clobber(outfil, clobber=False):
    """ Simple clobber check

    outfil : str
    clobber : bool, optional
    """
    # Chk clobber
    if os.path.isfile(outfil):
        if clobber:
            warnings.warn("Overwriting previous DB file {:s}".format(outfil))
            return True
        else:
            warnings.warn("Not overwiting previous DB file.  Set clobber=True to do so")
            return False
    else:
        return True


def get_build_groups(version):
    """
    Parameters
    ----------
    version : str

    Returns
    -------
    build_groups : dict

    """

    groups = OrderedDict()
    if version == 'v01':
        groups['FUV'] = fuv
    elif version == 'v02':
        pdb.set_trace()
        groups = OrderedDict()
    else:
        raise IOError("Not ready for this version")
    # Return
    return groups
