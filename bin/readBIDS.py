#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

import argparse
import bids
import json

def get_arguments():
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description="",
            epilog="""
            Create a json file from BIDS
            Input: BIDS folder
            Output: json file
            """)

    parser.add_argument(
            "-i", "--bids",
            required=True, nargs="+",
            help="BIDS folder",
            )

    parser.add_argument(
            "-o", "--json",
            required=True, nargs="+",
            help="json output file",
            )

    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        return args


class readBIDS(object):
    def __init__(self, bids, json):
        self.bids = bids[0]
        self.json = json[0]
        self.data = []
        #self.bidsValidated = True

    def run(self):
        self.ds = bids.BIDSLayout(self.bids)

        for nSub in self.ds.get_subjects():
            sessions = self.ds.get_sessions(subject=nSub)
            if sessions:
                dwis = []
                t1 = []
                fmaps = []
                for nSess in sessions:
                    dwis = self.ds.get(subject=nSub, session=nSess,
                                       datatype='dwi', extensions='nii.gz')
                    fmaps = self.ds.get(subject=nSub, session=nSess,
                                        datatype='fmap', extensions='nii.gz')
                    t1s = self.ds.get(subject=nSub, session=nSess,
                                     datatype='anat', extensions='nii.gz')

                    for nRun, dwi in enumerate(dwis):  # Possible runs
                        self.getData(nSub, dwi, fmaps, t1, nSess, nRun)
            else:
                dwis = self.ds.get(subject=nSub, datatype='dwi',
                                   extensions='nii.gz')
                fmaps = self.ds.get(subject=nSub, datatype='fmap',
                                    extensions='epi.nii.gz')
                t1s = self.ds.get(subject=nSub, datatype='anat',
                                 extensions='nii.gz')
                nSess = ''

                for nRun, dwi in enumerate(dwis):  # Possible runs
                    self.getData(nSub, dwi, fmaps, t1s, nSess, nRun)

            self.writeJson()

    def getData(self, nSub, dwi, fmaps, t1s, nSess, nRun):
        dwi_path = dwi.path
        bvec_path = self.ds.get_bvec(dwi_path)
        bval_path = self.ds.get_bval(dwi_path)

        try:
            dwi_PE = dwi.metadata['PhaseEncodingDirection']
            if len(dwi_PE) == 1:
                dwi_revPE = dwi_PE+'-'
            else:
                dwi_revPE = dwi_PE[0]
        except KeyError:
            # Crash violemment
            print('Violence')

        # Find b0 for topup, take the first one
        revb0_path = ''
        for nfmap in fmaps:
            if 'PhaseEncodingDirection' in nfmap.metadata and 'IntendedFor' in nfmap.metadata:
                fmap_PE = nfmap.metadata['PhaseEncodingDirection']
                refDWI = nfmap.metadata['IntendedFor']

                if fmap_PE == dwi_revPE:
                    if os.path.basename(refDWI) == dwi.filename:
                        revb0_path = nfmap.path
                        break

        totalreadout = 1
        try:
            dwi_RT = dwi.metadata['TotalReadoutTime']
            fmap_RT =   nfmap.metadata['TotalReadoutTime']
            if dwi_RT == fmap_RT:
                totalreadout = dwi_RT
        except:
            print('Pas grave')

        t1_path = 'todo'
        if len(t1s)==1:
            t1_path =  t1[0].path
        elif 'run' in dwi.path:
            for t1 in t1s:
                if 'run-'+str(nRun+1) in t1.path:
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
    """Let's go"""
    args = get_arguments()
    app = readBIDS(**vars(args))
    return app.run()


if __name__ == '__main__':
    sys.exit(main())
