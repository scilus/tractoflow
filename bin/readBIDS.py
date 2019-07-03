#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Create a json file from BIDS
Input: BIDS folder
Output: json file
"""

import os

import argparse
from bids import BIDSLayout
import json
import time


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
        self.layout = None
        self.associations = {}
        self.json = json
        self.data = []

    def run(self):
        self.layout = BIDSLayout(self.bids, index_metadata=False)
        subjects = self.layout.get_subjects()
        for nSub in subjects:
            dwis = self.layout.get(subject=nSub,
                          datatype='dwi', extension='nii.gz',
                          suffix='dwi')
            t1s = self.layout.get(subject=nSub,
                         datatype='anat', extension='nii.gz',
                         suffix='T1w')
            fmaps = self.layout.get(subject=nSub,
                            datatype='fmap', extension='nii.gz',
                            suffix='epi')
            bvals = self.layout.get(subject=nSub,
                            datatype='dwi', extension='bval',
                            suffix='dwi')
            bvecs = self.layout.get(subject=nSub,
                            datatype='dwi', extension='bvec',
                            suffix='dwi')
            self.get_dwi_associations(fmaps, bvals, bvecs)
            for nRun, dwi in enumerate(dwis):  # Possible runs
                self.getData(nSub, dwi, t1s, nRun)
        self.writeJson()

    def get_metadata(self, bf):
        filename = bf.path.replace(
                   '.' + bf.get_entities()['extension'], '')
        with open(filename + '.json', 'r') as handle:
            return json.load(handle)
    
    def get_dwi_associations(self, fmaps, bvals, bvecs):
        for bval in bvals:
            dwi_filename = os.path.basename(bval.path).replace('.bval', '.nii.gz')
            if dwi_filename not in self.associations.keys():
                self.associations[dwi_filename] = {"bval": bval.path}
            else:
                self.associations[dwi_filename]["bval"] = bval.path
        for bvec in bvecs:
            dwi_filename = os.path.basename(bvec.path).replace('.bvec', '.nii.gz')
            if dwi_filename not in self.associations.keys():
                self.associations[dwi_filename] = {"bvec": bvec.path}
            else:
                self.associations[dwi_filename]["bvec"] = bvec.path
        for fmap in fmaps:
            metadata = self.get_metadata(fmap)
            intended = [metadata.get('IntendedFor', [])]
            for target in intended:
                dwi_filename = os.path.basename(target)
                if dwi_filename not in self.associations.keys():
                    self.associations[dwi_filename]= {'fmap' : [fmap]}
                elif 'fmap' in self.associations[dwi_filename].keys():
                    self.associations[dwi_filename]['fmap'].append(fmap)
                else:
                    self.associations[dwi_filename]['fmap'] = [fmap]

    def getData(self, nSub, dwi, t1s, nRun):
            nSess=''
            if 'session' in dwi.get_entities().keys():
                nSess = dwi.get_entities()['session']
            dwi_path = dwi.path
            fmaps = []
            bval_path=''
            bvec_path=''
            if dwi.filename in self.associations.keys():
                if "bval" in self.associations[dwi.filename].keys():
                    bval_path = self.associations[dwi.filename]['bval']
                if "bvec" in self.associations[dwi.filename].keys():
                    bvec_path = self.associations[dwi.filename]['bvec']
                if "fmap" in self.associations[dwi.filename].keys():
                    fmaps = self.associations[dwi.filename]['fmap']

            dwi_PE = 'todo'
            dwi_revPE = -1
            conversion = {"i": "x", "j": "y", "k": "z"}
            dwi_metadata = self.get_metadata(dwi)
            if 'PhaseEncodingDirection' in dwi_metadata:
                dwi_PE = dwi_metadata['PhaseEncodingDirection']
                dwi_PE = dwi_PE.replace(dwi_PE[0], conversion[dwi_PE[0]])
                if len(dwi_PE) == 1:
                    dwi_revPE = dwi_PE+'-'
                else:
                    dwi_revPE = dwi_PE[0]

            # Find b0 for topup, take the first one
            revb0_path = ''
            totalreadout = ''
            for nfmap in fmaps:
                nfmap_metadata = self.get_metadata(nfmap)
                if 'PhaseEncodingDirection' in nfmap_metadata and\
                'TotalReadoutTime' in dwi_metadata and\
                'TotalReadoutTime' in nfmap_metadata:

                    fmap_PE = nfmap_metadata['PhaseEncodingDirection']
                    fmap_PE = fmap_PE.replace(fmap_PE[0], conversion[fmap_PE[0]])
                    dwi_RT = dwi_metadata['TotalReadoutTime']
                    fmap_RT = nfmap_metadata['TotalReadoutTime']
                    if fmap_PE == dwi_revPE and dwi_RT == fmap_RT:
                        revb0_path = nfmap.path
                        totalreadout = dwi_RT
                        break

            t1_path = 'todo'
            t1_nSess=[]
            for t1 in t1s:
                if 'session' in t1.get_entities().keys() and\
                    t1.get_entities()['session'] == nSess:
                    t1_nSess.append(t1)
                elif 'session' not in t1.get_entities().keys():
                    t1_nSess.append(t1)

            if len(t1_nSess) == 1:
                t1_path = t1_nSess[0].path
            elif 'run' in dwi.path:
                for t1 in t1_nSess:
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
            json.dump(self.data,
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
