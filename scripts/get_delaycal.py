#! /usr/bin/env python

"""
Helper script to find the CBID of a valid delaycal (noise diode) 
observation assocated with a given observation. The CBID of the
valid delay calibration can be search for in the archive and 
its URL (with a token) used as input to the KATCalPipe.py script
with --polcal uding the option --delaycal_mvf
"""

import argparse
import katdal


def _get_noise_diode_scan(kd):
    # Get the scan we want from the delaycal observation.
    # A delaycal has two tracks and the second one should have CompScanLabel='corrected'
    kd.select(scans='track,stop', compscans='corrected')
    # Now there should only be the scans we want, find the longest one
    scan_select = -1
    max_scan = 0
    output = None 
    for scan, _, target in kd.scans():
        scan_length = len(kd.dumps)
        if scan_length > max_scan:
            max_scan, scan_select = scan_length, scan
    if scan_select >= 0:
        output = f'{kd.obs_params["capture_block_id"]:13s} {scan:<5d} {max_scan:<6d} {target.name}'
    return output


parser = argparse.ArgumentParser()
parser.add_argument("katdata", help="URL of observation for which to find the associated delaycal CBID.")

args = parser.parse_args()

base_katdata = katdal.open(args.katdata)
sb_cbids = base_katdata.source.telstate.get_range('sdp_capture_block_id', st=0)
delay_cbids = [cbid[0] for cbid in sb_cbids[::-1]
                       if 'calibrate_delays.py' in 
                       base_katdata.source.telstate.view(cbid[0])['obs_params']['script_name']]
# Check each candidate delaycal CBID for an appropriate noise-diode scan
out = []
for cbid in delay_cbids:
    dc_katdata = katdal.open(args.katdata, capture_block_id=cbid)
    dc_search = _get_noise_diode_scan(dc_katdata)
    if dc_search is not None:
        out += [dc_search]

if out == []:
    print(f'No valid delay calibration observations found.')
else:
    print(f'Valid delay calibration observations with noise diode:\n\n' \
          f'{"CBID":13s} {"Scan":5s} {"Dumps":6s} {"Target"}\n' \
          f'{"=" * 40}')
    for line in out:
        print(f'{line}')
