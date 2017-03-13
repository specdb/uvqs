""" Module to ingest SDSS III (aka BOSS) data products
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import os, json
import pdb

from astropy.table import Table, Column, vstack
from astropy.time import Time
from astropy.io import fits

from linetools import utils as ltu
from linetools.spectra import io as lsio

from specdb.build.utils import chk_for_duplicates
from specdb.build.utils import chk_meta
from specdb.build.utils import init_data


def grab_meta():
    """ Grab meta Table

    Returns
    -------
    boss_meta : Table

    """

    #http://www.sdss.org/dr12/algorithms/boss-dr12-quasar-catalog/
    uvqs_fuv = Table.read(os.getenv('DROPBOX_DIR')+'/z1QSO/database/uvq_dr1_v1.7.fits')
    nuvqs = len(uvqs_fuv)
    # DATE-OBS -- Grab from header
    #t = Time(list(uvqs_fuv['MJD'].data), format='mjd', out_subfmt='date')  # Fixes to YYYY-MM-DD
    #boss_meta.add_column(Column(t.iso, name='DATE-OBS'))
    # Add columns -- From specfiles
    #boss_meta.add_column(Column(['BOSS']*nboss, name='INSTR'))
    #boss_meta.add_column(Column(['BOTH']*nboss, name='GRATING'))
    #boss_meta.add_column(Column([2100.]*nboss, name='R'))  # RESOLUTION
    #boss_meta.add_column(Column(['SDSS 2.5-M']*nboss, name='TELESCOPE'))
    # Redshift logic
    uvqs_fuv['zem_GROUP'] = uvqs_fuv['Z']
    uvqs_fuv['sig_zem'] = uvqs_fuv['Z_SIG']
    uvqs_fuv['flag_zem'] = [str('UVQS')]*nuvqs
    # Rename RA/DEC
    uvqs_fuv.rename_column('RA', 'RA_GROUP')
    uvqs_fuv.rename_column('DEC', 'DEC_GROUP')
    # STYPE
    uvqs_fuv['STYPE'] = [str('QSO')]*nuvqs
    # Check
    assert chk_meta(uvqs_fuv, chk_cat_only=True)
    # Return
    return uvqs_fuv


def hdf5_adddata(hdf, sname, meta, debug=False, chk_meta_only=False, **kwargs):
    """ Add BOSS data to the DB

    Parameters
    ----------
    hdf : hdf5 pointer
    IDs : ndarray
      int array of IGM_ID values in mainDB
    sname : str
      Survey name
    chk_meta_only : bool, optional
      Only check meta file;  will not write
    boss_hdf : str, optional


    Returns
    -------

    """
    from specdb import defs
    from astropy.time import Time
    # Init
    Rdicts = defs.get_res_dicts()
    dpath = os.getenv('DROPBOX_DIR') + 'z1QSO/database/1D_Spectra/'
    # Add Survey
    print("Adding {:s} survey to DB".format(sname))
    uvqs_grp = hdf.create_group(sname)

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 30000  # Just needs to be large enough
    data = init_data(max_npix, include_co=False)
    # Init
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    spec_set.resize((nspec,))
    wvminlist = []
    wvmaxlist = []
    speclist = []
    npixlist = []
    instrlist = []
    gratinglist = []
    Rlist = []
    timelist = []
    telelist = []
    # Loop
    maxpix = 0
    for jj,row in enumerate(meta):
        # Generate full file
        full_file = dpath+row['SPEC_FILE']
        #print("Reading {:s}".format(full_file))
        # Read
        if 'CAHA' in full_file:
            spec = lsio.readspec(full_file, masking='none')
        else:
            spec = lsio.readspec(full_file)
        # npix
        try:
            npix = spec.npix
        except ValueError:  # Bad CAHA data
            npix = spec.flux.size
        maxpix = max(npix,maxpix)
        # Parse name
        fname = full_file.split('/')[-1]
        # Fill
        for key in ['wave','flux','sig']:
            data[key] = 0.  # Important to init (for compression too)
        data['flux'][0][:npix] = spec.flux.value
        data['sig'][0][:npix] = spec.sig.value
        data['wave'][0][:npix] = spec.wavelength.value
        # Meta
        speclist.append(str(fname))
        wvminlist.append(np.min(data['wave'][0][:npix]))
        wvmaxlist.append(np.max(data['wave'][0][:npix]))
        npixlist.append(npix)
        if 'LCO' in full_file:
            instrlist.append('duPont-BCS')
            telelist.append('duPont')
            gratinglist.append('600/5000') # Not accurate for all data, I suspect
            Rlist.append(1200.)
            timelist.append(spec.header['UT-DATE']+'T'+spec.header['UT-TIME'])
        elif 'Lick' in full_file:
            instrlist.append('Kast')
            gratinglist.append('Both')
            Rlist.append(1000.)
            telelist.append('Lick-3m')
            timelist.append(spec.header['DATE-OBS'])
        elif 'CAHA' in full_file:
            instrlist.append('CAFOS')
            telelist.append('CAHA')
            gratinglist.append('??')
            Rlist.append(1000.)
            if 'Aug' in row['OBS_DATE']:
                timelist.append(row['OBS_DATE'][-4:]+'-08-01')
            else:
                pdb.set_trace() # TIME
        elif 'Keck' in full_file:
            instrlist.append('ESI')
            telelist.append('Keck-II')
            gratinglist.append('ECH')
            Rlist.append(Rdicts['ESI'][spec.header['SLMSKNAM']])
            timelist.append(spec.header['DATE-OBS']+'T'+spec.header['UT'])
        elif 'MMT' in full_file:
            instrlist.append('mmtbluechan')
            telelist.append('MMT')
            gratinglist.append('??')
            Rlist.append(1000.)
            timelist.append(spec.header['DATE-OBS']+'T'+spec.header['UT'])
        elif 'Magellan' in full_file:
            instrlist.append('MagE')
            telelist.append('Magellan')
            gratinglist.append('N/A')
            Rlist.append(5000.)
            timelist.append(spec.header['UT-DATE']+'T'+spec.header['UT-TIME'])
        else:
            pdb.set_trace()
        # Only way to set the dataset correctly
        spec_set[jj] = data

    #
    print("Max pix = {:d}".format(maxpix))
    # Add columns
    t = Time(timelist, out_subfmt='date')  # Fixes to YYYY-MM-DD
    meta.add_column(Column(t.iso, name='DATE-OBS'))
    #meta.add_column(Column(speclist, name='SPEC_FILE'))
    meta.add_column(Column(npixlist, name='NPIX'))
    meta.add_column(Column(wvminlist, name='WV_MIN'))
    meta.add_column(Column(wvmaxlist, name='WV_MAX'))
    meta.add_column(Column(Rlist, name='R'))
    meta.add_column(Column(gratinglist, name='DISPERSER'))
    meta.add_column(Column(telelist, name='TELESCOPE'))
    meta.add_column(Column(instrlist, name='INSTR'))
    meta.add_column(Column(np.arange(nspec,dtype=int),name='GROUP_ID'))
    meta.add_column(Column([2000.]*len(meta), name='EPOCH'))

    # Add HDLLS meta to hdf5
    if chk_meta(meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = meta
    else:
        pdb.set_trace()
        raise ValueError("meta file failed")
    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2016AJ....152...25M',
                 bib='uvqs'),
            ]
    jrefs = ltu.jsonify(refs)
    hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return


def add_ssa(hdf, dset):
    """  Add SSA info to meta dataset
    Parameters
    ----------
    hdf
    dset : str
    """
    from specdb.ssa import default_fields
    Title = '{:s}: UVQS FUV Quasars'.format(dset)
    ssa_dict = default_fields(Title, flux='flambda', fxcalib='RELATIVE')
    hdf[dset]['meta'].attrs['SSA'] = json.dumps(ltu.jsonify(ssa_dict))
