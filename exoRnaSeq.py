#! /usr/bin/python3
import os
import sys
import csv
import json
import glob
import time
import logging
import argparse
import pickle
import subprocess
from collections import namedtuple
from collections import defaultdict
from collections import Counter
from difflib import SequenceMatcher
import vartest.variant as variant
import numpy as np
from sklearn.metrics import jaccard_similarity_score

class exoRnaProject:
    def __init__(self, inpath=None, outpath=None):
        if inpath == None:
            self.inpath = '/data/projects/exosome/macrophage_RNA_Seq'
        else:
            self.inpath = os.path.abspath(inpath)
        if outpath == None:
            self.outpath = '{0}/Outputs'.format(self.inpath)
        else:
            self.outpath = os.path.abspath(outpath)
        

class topHat:
    def __init__(self, outpath=None, tophat=None, gtf=None, btindex=None,
        name=None, threads=None):
        if name == None:
            self.name = 'exomsome'
        else:
            self.name = name
        if outpath == None:
            self.outpath = '{0}/topHat/{1}'.format(os.getcwd(), self.name)
            if not os.path.exists('{0}/topHat'.format(os.getcwd())):
                os.mkdir('{0}/topHat'.format(os.getcwd()))
        else:
            self.outpath = '{0}/{1}'.format(outpath, self.name)
        if tophat == None:
            self.tophat = '/projects/home/sravishankar9/tools/tuxedo_suite/tophat-2.0.13.Linux_x86_64/tophat2'
        else:
            self.tophat = tophat
        if gtf == None:
            self.gtf = '/data/db/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
        else:
            self.gtf = gtf
        if btindex == None:
            self.btindex = '/data/db/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome'
        else:
            self.btindex = btindex
        if threads == None:
            self.threads = 4
        else:
            self.threads = threads
        return

    def runTopHat(self, inpath, top_params=None):
        if top_params == None:
            top_params = ['--no-coverage-search', '--b2-very-sensitive',
                '--no-novel-juncs']
        tcmd = [self.tophat] + top_params + ['-p', self.threads, '-G',
            self.gtf, '-o', self.outpath, self.btindex] + inpath
        trun = subprocess.Popen(tcmd, stdout=PIPE, stderr=PIPE, shell=False)
        tlog = trun.communicate()
        if trun.returncode != 0:
            print('TopHat failed for sample: {0}'.format(self.name))
            sys.exit()
        else:
            return(self.outpath, trun.returncode)

class Cufflinks:
    def __init__(self, outpath=None, cuffpath=None, gtf=None, btindex=None,
        name=None, threads=None):
        if name == None:
            self.name = 'exomsome'
        else:
            self.name = name
        if outpath == None:
            self.outpath = '{0}/Cufflinks/{1}'.format(os.getcwd(), self.name)
            if not os.path.exists('{0}/Cufflinks'.format(os.getcwd())):
                os.mkdir('{0}/Cufflinks'.format(os.getcwd()))
            self.doutpath = '{0}/Cuffdiff'.format(self.outpath)
            self.moutpath = '{0}/Cuffmerge'.format(self.outpath)
        else:
            self.outpath = '{0}/{1}'.format(outpath, self.name)
            self.doutpath = '{0}/Cuffdiff'.format(self.outpath)
            self.moutpath = '{0}/Cuffmerge'.format(self.outpath)
        if cuffpath == None:
            self.cuffpath = '/projects/home/sravishankar9/tools/tuxedo_suite/cufflinks/'
            self.cufflinks = '{0}cufflinks'.format(self.cuffpath)
            self.cuffdiff = '{0}cuffdiff'.format(self.cuffpath)
            self.cuffmerge = '{0}cuffmerge'.format(self.cuffpath)
        else:
            self.cuffpath = cuffpath
            self.cufflinks = '{0}cufflinks'.format(self.cuffpath)
            self.cuffdiff = '{0}cuffdiff'.format(self.cuffpath)
            self.cuffmerge = '{0}cuffmerge'.format(self.cuffpath)
        if gtf == None:
            self.gtf = '/data/db/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
        else:
            self.gtf = gtf
        if btindex == None:
            self.btindex = '/data/db/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome'
            self.genome = '/data/db/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.fa'
        else:
            self.btindex = btindex
        if threads == None:
            self.threads = 4
        else:
            self.threads = threads
        return

        def runCufflinks(self, inpath, clinks_params=None):
            if clinks_params == None:
                clinks_params = []
            ccmd = [self.cufflinks, '-p', self.threads, '-G', self.gtf,
                '-o', self.outpath, inpath] + clinks_params
            crun = subprocess.Popen(ccmd, stdout=PIPE, stderr=PIPE, shell=False)
            clog = crun.communicate()
            if trun.returncode != 0:
                print('Cufflinks failed for sample: {0}'.format(self.name))
                sys.exit()
            else:
                return(self.outpath, crun.returncode)

        def runCuffMerge(self, inpath, cmerge_params=None):
            if cmerge_params == None:
                cmerge_params = []
            mcmd = [self.cuffmerge, '-g', self.gtf, '-s' self.genome, '-p',
                self.threads, '-o', self.moutpath, inpath]
            mrun = subprocess.Popen(mcmd, stdout=PIPE, stderr=PIPE, shell=False)
            mlog = mrun.communicate()
            if trun.returncode != 0:
                print('Cuffmerge failed for sample: {0}'.format(self.name))
                sys.exit()
            else:
                return(self.moutpath, mrun.returncode)

        def runCuffdiff(self, inpath, cdiff_params=None):
            if cdiff_params == None:
                cdiff_params = []
            dcmd = [self.cuffdiff, '-o', self.doutpath, '-b', self.genome,
                '-p', self.threads] + cdiff_params + inpath
            drun = subprocess.Popen(dcmd, stdout=PIPE, stderr=PIPE, shell=False)
            dlog = drun.communicate()
            if trun.returncode != 0:
                print('Cuffdiff failed for sample: {0}'.format(self.name))
                sys.exit()
            else:
                return(self.doutpath, drun.returncode)
