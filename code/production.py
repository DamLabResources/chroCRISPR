import pandas as pd
import numpy as np
from os.path import abspath, dirname, pardir, exists
import os
from tempfile import TemporaryDirectory, NamedTemporaryFile
from subprocess import check_call, check_output
import shlex
from Bio.Alphabet import generic_dna, RNAAlphabet, DNAAlphabet
import gzip
from Bio.Seq import Seq
import argparse, argcomplete
import sys
## U2OS used HG19 ,which does not have string 'chr' in fasta name
## HEK DNaseSeq data was mapped from NCBI SRA downloads, which has fasta name 'chr1, chr2 ...'

## GUIDE-Seq output used 'chr1, chr2 ...', while CIRCLE-Seq output used '1, 2, ...' for names of chromosomes.
## Of course every one has its own naming system....

class DNaseCrispr():
    def __init__(self, bona_file,  cell, window=100, bam_file=None, keep=False):
        self.bam_file = bam_file
        self.cell_read(cell)
        self.window = window
        self.keep = keep
        self.bona_file = bona_file
        self.bonaFileRead()


    def bonaFileRead(self):
        tag = True
        if "GUIDE" in self.bona_file:
            self.df = self.guidePro()
        elif "CIRCLE" in self.bona_file:
            self.df = self.circlePro()
        elif "S2" in self.bona_file:
            self.df = self.circleS2()
        else:
            tag =False
            self.dfPickCol = self.tctome()
            #raise NameError("Input file name has to have key word 'GUIDE' or 'CIRCLE'.")


        #if self.keep:
        #    self.dfPickCol = self.df.loc[:, ['chr', 'start', 'stop', 'name', 'read_count',
        #                                     'gRNA_seq', 'offsite_seq', 'Mismatches']]
            #self.dfPickCol = self.dfPickCol.rename(columns={'start': 'bStart',
            #                               'stop': 'bEnd'})

        if tag:
            self.df['bStart'] = self.df['Position'] - int(np.ceil(self.window / 2))
            self.df['bEnd'] = self.df['Position'] + int(np.ceil(self.window / 2))
            self.df['bStart'] = self.df.bStart.astype('int')
            self.df['bEnd'] = self.df.bEnd.astype('int')
            self.dfPickCol = self.df.loc[:, ['chr', 'bStart', 'bEnd', 'name', 
                                             'read_count', 'trc', 'cpm',
                                                 'start', 'stop', 'target_name',
                                                 'gRNA_seq', 'offsite_seq', 'Mismatches',
                                             ]]

        if self.bam_file:
            with NamedTemporaryFile('w', dir='/tmp', delete=False) as handle:
                self.dfPickCol.to_csv(handle.name, header=False, sep='\t', index=False)
                print('call DNase depth from bam file')
                sumCov = self.dnase_depth(self.dfPickCol, self.bam_file, handle.name)
            self.dfPickCol['sumCov'] = sumCov
            self.dfPickCol['normed_sumCov'] = self.dfPickCol['sumCov'] / self.bam_num_read(self.bam_file) / (self.dfPickCol['bEnd'] - self.dfPickCol['bStart'])
            if 'read_count' in list(self.dfPickCol.columns):
                self.dfPickCol['normed_read_count'] = self.dfPickCol['read_count'] / self.dfPickCol['read_count'].sum()


        # get cutting scores
        if self.keep:
            self.dfPickCol['MIT'] = self.cut_effi('MIT')
            self.dfPickCol['CFD'] = self.cut_effi('CFD')
            self.dfPickCol['Kinetic'] = self.cut_effi('Kinetic')

    def cell_read(self, cell):

        if cell in ["HEK", "U2OS", "HEK293T"]:
            self.cell = cell
        else:
            raise ValueError("Input cell lines available now are HEK293T and U2OS.")


    def guidePro(self):
        df = pd.read_csv(self.bona_file, header=0, index_col=False)
        # exclude non-true hit in the guideseq package output, which is the column 'Off-Target Sequence' with empty content.
        df = df[~df['Off-Target Sequence'].isnull()]
        df = df.rename(columns={'#BED Chromosome': 'chr', 'bi.sum.mi': 'read_count',
                                'BED Min.Position': 'start', 'BED Max.Position': 'stop',
                                'BED Name': 'name', 'Cells': 'target_cells', 'Position':'oriPos',
                                'Target Sequence': 'gRNA_seq', 'Targetsite': 'target_name',
                                'Off-Target Sequence': 'offsite_seq',
                                })
        if self.cell == 'U2OS':
            df['chr'] = df['chr'].str.replace('chr', '')
        
        groups = {'11449328': ['VEGFA_site2', 'VEGFA_site3', 'VEGFA_site1', 'EMX1'],
                  '15471209': ['FANCF', 'RNF2'],
                  '7287344': ['HEKgRNA4', 'HEKgRNA1', 'HEKgRNA3', 'HEKgRNA2'],
                  '100': ['PGP1_TAZ']}

        trc = [] # list of total read count from GUIDE-Seq NGS raw read count
        for idx, row in df.iterrows():
            for g in groups.keys():
                if row['target_name'] in groups[g]:
                    trc.append(int(g))

        df['trc'] = trc
        # Cleavage Per Million identified in GUIDE-Seq
        df['cpm'] = df['read_count'] /df['trc'] * 10e6  
        
        # select cell-specific hits
        df = df[df['target_cells'].str.contains(self.cell)]

        pos = []
        for idx, row in df.iterrows():
            if row.Strand == '+':
                pos.append(row['BED off-target end'] - 5)
            elif row.Strand == '-':
                pos.append(row['BED off-target start'] + 5)
        df['Position'] = pos

        return df



    def circlePro(self):
        df = pd.read_csv(self.bona_file, header=-1, sep='\t', index_col=False)
        circleCol = ['chr', 'start', 'stop', 'name', 'read_count', 'strand',
                     'iv', 'iv.chrom', 'iv.start', 'iv.end', 'window_sequence', 'sequence', 'distance', 'length',
                     'filename', 'target_name', 'target_cells', 'full_name', 'target_sequence', 'realigned_target',
                     'pval_pos', 'pval_nar', 'pval_one', 'control_pval_pos',
                     'control_pval_nar', 'control_pval_one']
        df.columns = circleCol
        df = df.rename(columns={'target_sequence': 'gRNA_seq',
                                'sequence': 'offsite_seq',
                                'distance': 'Mismatches'})
        pos = []
        for idx, row in df.iterrows():
            if row.strand == '+':
                pos.append(row['stop']-5)
            elif row.strand == '-':
                pos.append(row['start']+5)
        df['Position'] = pos
        if self.cell == 'HEK':
            df['chr'] = 'chr' + df['chr'].astype(str)

        groups = {'11449328': ['U2OS_combined_VEGFA_site_2', 'U2OS_combined_VEGFA_site_3', 'U2OS_exp2_VEGFA_site_1', 'U2OS_exp2_EMX1'],
                  '11557513': ['U2OS_exp2_FANCF', 'U2OS_exp2_RNF2'],
                  '7287344': ['HEK293_combined_Adli_site_4', 'HEK293_Adli_site_2', 'HEK293_Adli_site_3', 'HEK293_Adli_site1'],
                  '100': ['PGP1_TAZ']}
        aa=0
        trc = [] # list of total read count from GUIDE-Seq NGS raw read count
        for idx, row in df.iterrows():
            for g in groups.keys():
                if row['target_name'] in groups[g]:
                    trc.append(int(g))

        df['trc'] = trc
        # Cleavage Per Million identified in CIRCLE-Seq
        df['cpm'] = df['read_count'] /df['trc'] * 10e6
        
        # select cell-specific hits
        df = df[df['target_cells'].str.contains(self.cell)]
        
        return df

    def circleS2(self): #
        df = pd.read_excel(self.bona_file)
        df = df.rename(columns={'Chromosome': 'chr',
                                'Start': 'start',
                                'End': 'stop',
                                'Summary': 'name',
                                'Read': 'read_count',
                                'Strand': 'strand',
                                'Off-target Sequence': 'offsite_seq',
                                'Distance': 'Mismatches',
                                'Targetsite': 'target_name',
                                'TargetSequence': 'gRNA_seq',})
        pos = []
        for idx, row in df.iterrows():
            if row.strand == '+':
                pos.append(row['stop'] - 5)
            elif row.strand == '-':
                pos.append(row['start'] + 5)
        df['Position'] = pos
        if self.cell == 'HEK':
            df['chr'] = 'chr' + df['chr'].astype(str)

        # Take out hits with offsite_seq is empty
        df = df[~df['offsite_seq'].isnull()]
        # select cell-specific hits
        df = df[df['Cell'].str.contains(self.cell)]

        return df

    def tctome(self):
        if ".gz" in self.bona_file:
            handle = gzip.open(self.bona_file, 'tr')
        else:
            handle = open(self.bona_file, 'r')

        df = pd.read_csv(handle, header=-1, sep='\t')
        cols = df.columns
        df.columns = ['chr', 'bStart', 'bEnd', 'name', 'score', 'strand'] + list(np.arange(6,df.shape[1]))
        if self.cell == 'U2OS':
            df['chr'] = df['chr'].str.replace('chr', '')
        return df

    def dnase_depth(self, dfHits, bam_file, region_file):
        print(region_file)
        cmd = 'sambamba depth base -t 7 --min-coverage=0 -L %s %s' % (region_file, bam_file)
        print(cmd)
        sumCov = []
        with NamedTemporaryFile('w', dir='/tmp', delete=False) as handle:
            check_call(shlex.split(cmd), stdout=handle)
            #handle.write(check_output(shlex.split(cmd)))
            df = pd.read_csv(handle.name, header=0, sep='\t')
            print(df.shape)
            print('depth file is %s ' % handle.name)
            print('reading base depth file generated by sambamba')
            sumCov = []
            for idx, row in dfHits.iterrows():
                #print('row is %s %s %s' % (row['bStart'], row['bEnd'], row['chr']))
                su = df[(df['POS'] <= row['bEnd']) & (df['POS'] >= row['bStart'])]# & (df['REF'] == row['chr'])]
                #print(su.shape)
                #su = df[(df['POS'] <= row['bEnd']) & (df['POS'] >= row['bStart']) & (df['REF'] == row['chr'])].COV.sum()
                #sumCov.append(su)
                if su.shape[0] == 0:
                    sumCov.append(0)
                else:
                    sumCov.append(su.COV.sum())
            print(len(sumCov))
        return sumCov

    def bam_num_read(self, bam_file):
        cmd = 'sambamba view -c %s' % bam_file
        num = int(check_output(shlex.split(cmd)).decode("utf-8"))
        print('read number in bam file is %s' % num)
        return (num/1000000.0)


    def shuffle(self):
        cols = ['chr', 'bStart', 'bEnd', 'name', 'read_count', 'sumCov']
        dfTotalShuf = pd.Dataframe(columns=cols)
        with NamedTemporaryFile('w', dir='/tmp', delete=False) as f:
            for i in range(1):
                self.dfPickCol.to_csv(f.name, sep='\t', header=False, index=False)
                if self.cell == 'HEK':
                    region_file = '../raw_data/hg19.chromSize.noChrY.chr.sorted.bed'
                elif self.cell == 'U2OS':
                    region_file = '../raw_data/hg19.chromSize.noChrY.num.sorted.bed'
                cmd = 'bedtools shuffle -seed %s -i %s -g %s' % (i, f.name, region_file)

                print(cmd)
                check_call(shlex.split(cmd), stdout = f)
                dfShuf = pd.read_csv(f.name, header=-1, sep='\t')
                dfShuf.columns = cols
                dfShuf['bStart'] = dfShuf.bStart.astype('int')
                dfShuf['bEnd'] = dfShuf.bEnd.astype('int')
                print(dfShuf.head())
                sumCov = self.dnaseDepth(dfShuf, self.bam_file, f.name)
                dfShuf['sumCov'] = sumCov
                dfShuf.to_csv(f.name, sep='\t', header=False, index=False)
                dfTotalShuf = pd.concat(dfTotalShuf, dfShuf, ignore_index=True)

        return dfTotalShuf

    def cut_effi(self, method):
        # calculate MIT and CFD using crispr_tree
        import sys
        rpath = dirname(dirname(dirname(dirname(abspath(__file__)))))
        crispr_path = rpath + '/submodules/crseek'
        sys.path.append(crispr_path)
        from crseek import estimators, preprocessing, utils
        # Evaluate the built-in models on the same data.
        ests = {'MIT': estimators.MITEstimator.build_pipeline(),
                'CFD': estimators.CFDEstimator.build_pipeline(strict=False),
                'Kinetic': estimators.KineticEstimator.build_pipeline()}
        est = ests[method]
        seqs = []
        for idx, row in self.dfPickCol.iterrows():
            on = row['gRNA_seq']
            off = row['offsite_seq']
            #print(on ,off)
            off = off.replace('-', 'N')
            on = on[:20]
            if len(off)<23:
                off = 'N'*(23-len(off))+off
                #seqs.append([on[:20], 'N'*(23-len(off))+off])
            elif len(off) > 23:
                off = off[-23:]
                #seqs.append([on[:20], off[-23:]])
            on = on.replace('T', 'U')

            seqs.append([on, off])

        cutting_data = pd.DataFrame(seqs)
        cutting_data.columns = ['spacer', 'target']
        cutting_data['target'] = list(utils.smrt_seq_convert('Seq',
                                                             cutting_data['target'].values,
                                                             alphabet=generic_dna))

        spacer_dna = utils.smrt_seq_convert('Seq',
                                            cutting_data['spacer'].values,
                                            alphabet=generic_dna)
        cutting_data['spacer'] = [s.transcribe() for s in spacer_dna]
        data = cutting_data[['spacer', 'target']].values
        scores = est.predict_proba(data)


        return scores

if __name__ == '__main__':
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('ERROR: %s\n' % message)
            self.print_help()
            sys.exit(2)
    
    parser = MyParser()
    parser.add_argument('--regions', '-L', dest='L',
                        required=True,
                        help="List of chromosomal regions from [GUIDE-Seq, CIRCLE-Seq, Human gene transcripts]")
    parser.add_argument('--celltype', '-c', dest='c',
                        required=True,
                        help="Cell type of DNase-Seq. Options: ['HEK293T', 'U2OS']")
    parser.add_argument('--bam', '-b', dest='b',
                        required=True,
                        help="Sorted DNase-Seq bam file generated from preprocessing.py")
    parser.add_argument('--window', '-w', dest='w',
                        help="Window size of DNase-Seq flanking cleavage site (CRISPR) or transcription start site (transcriptome). Default: 100bp for CRISPR")
    
    parser.add_argument('--output', '-o', dest='o',
                        help="Output file with DNase-Seq RPM, CPM, and predicted cleavage efficiency.")
    argcomplete.autocomplete(parser)
    args = parser.parse_args()
    #print(args.c)
    #print(type(args.c))
    if not exists(args.L):
        parser.error("Region file does not exist.\n")
    elif not exists(args.b):
        parser.error("DNase-Seq bam file does not exist.\n")
    else:
        
    #hitFile = '../raw_data/GUIDEseq_allgRNAs_identified.csv'
    #bamFile = '../raw_data/HEK293T/HEK.se50.DNaseSeq.sorted.bam'
        wrap = DNaseCrispr(bona_file=args.L, cell=args.c, bam_file=args.b, keep=False)
    
    if args.o:
        if not exists(dirname(args.o)):
            print(dirname(args.o))
            os.makedirs(dirname(args.o))
        wrap.dfPickCol.to_csv(args.o, index=False)
    #wrap.dnaseDepth(wrap.bam_file, wrap.dfPickCol)
    #wrap.shuffle()

