#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Create a json file from BIDS
Input: BIDS folder
Output: json file
"""

import os

import argparse
import bids
import json


def _build_args_parser():
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=__doc__)

    parser.add_argument(
            "bids",
            help="BIDS folder")

    parser.add_argument(
            "json",
            help="json output file")

    return parser


class readBIDS(object):
    def __init__(self, bids, json):
        self.bids = bids
        self.json = json
        self.data = []

    def run(self):
        self.ds = bids.BIDSLayout(self.bids)

        for nSub in self.ds.get_subjects():
            sessions = self.ds.get_sessions(subject=nSub)
            if sessions:
                for nSess in sessions:
                    dwis = self.ds.get(subject=nSub, session=nSess,
                                       datatype='dwi', extensions='nii.gz',
                                       suffix='dwi')
                    fmaps = self.ds.get(subject=nSub, session=nSess,
                                        datatype='fmap', extensions='nii.gz',
                                        suffix='epi')
                    t1s = self.ds.get(subject=nSub, session=nSess,
                                      datatype='anat', extensions='nii.gz',
                                      suffix='T1w')
                    
                    for nRun, dwi in enumerate(dwis):  # Possible runs
                        self.getData(nSub, dwi, fmaps, t1s, nSess, nRun)

            else:
                dwis = self.ds.get(subject=nSub, datatype='dwi',
                                   extensions='nii.gz', suffix='dwi')
                fmaps = self.ds.get(subject=nSub, datatype='fmap',
                                    extensions='epi.nii.gz', suffix='epi')
                t1s = self.ds.get(subject=nSub, datatype='anat',
                                  extensions='nii.gz', suffix='T1w')
                nSess = ''

                for nRun, dwi in enumerate(dwis):  # Possible runs
                    self.getData(nSub, dwi, fmaps, t1s, nSess, nRun)

            self.writeJson()

    def getData(self, nSub, dwi, fmaps, t1s, nSess, nRun):
        dwi_path = dwi.path
        bvec_path = self.ds.get_bvec(dwi_path)
        bval_path = self.ds.get_bval(dwi_path)

        dwi_PE = 'todo'
        dwi_revPE = -1
        conversion = {"i": "x", "j": "y", "k": "z"}
        if 'PhaseEncodingDirection' in dwi.get_metadata():
            dwi_PE = dwi.get_metadata()['PhaseEncodingDirection']
            dwi_PE = dwi_PE.replace(dwi_PE[0], conversion[dwi_PE[0]])
            if len(dwi_PE) == 1:
                dwi_revPE = dwi_PE+'-'
            else:
                dwi_revPE = dwi_PE[0]

        # Find b0 for topup, take the first one
        revb0_path = ''
        totalreadout = ''
        for nfmap in fmaps:
            if 'PhaseEncodingDirection' in nfmap.get_metadata() and\
               'IntendedFor' in nfmap.get_metadata() and\
               'TotalReadoutTime' in (dwi.get_metadata() and nfmap.get_metadata()):
                refDWI = nfmap.get_metadata()['IntendedFor']

                if os.path.basename(refDWI) == dwi.filename:
                    fmap_PE = nfmap.get_metadata()['PhaseEncodingDirection']
                    fmap_PE = fmap_PE.replace(fmap_PE[0], conversion[fmap_PE[0]])
                    dwi_RT = dwi.get_metadata()['TotalReadoutTime']
                    fmap_RT = nfmap.get_metadata()['TotalReadoutTime']
                    totalreadout = 'todo'
                    if fmap_PE == dwi_revPE:
                        revb0_path = nfmap.path
                        if dwi_RT == fmap_RT:
                            totalreadout = dwi_RT
                        break

        t1_path = 'todo'
        if len(t1s) == 1:
            t1_path = t1s[0].path
        elif 'run' in dwi.path:
            for t1 in t1s:
                if 'run-' + str(nRun + 1) in t1.path:
                    t1_path = t1.path

        self.data.append({'subject': nSub,
                          'session': nSess,
                          'run': nRun,
                          't1': t1_path,
                          'dwi': dwi_path,
                          'bvec': bvec_path,
                          'bval': bval_path,
                          'rev_b0': revb0_path,
                          'DWIPhaseEncodingDir': dwi_PE,
                          'TotalReadoutTime': totalreadout})

    def writeJson(self):
        with open(self.json, 'w') as outfile:
            json.dump(self.data[::-1],
                      outfile,
                      indent=4,
                      separators=(',', ': '),
                      sort_keys=True)
            # add trailing newline for POSIX compatibility
            outfile.write('\n')


def main():
    parser = _build_args_parser()
    args = parser.parse_args()
    BIDS_app = readBIDS(args.bids, args.json)
    BIDS_app.run()


if __name__ == '__main__':
    main()
