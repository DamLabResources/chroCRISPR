from subprocess import check_call
import ruffus
import tasks
import os
import sys
import argparse
# Read mapping preprocessing for DNase-seq on U2OS and HEK293T
# Using GRCh37 (hg19) as ref since GUIDE-seq and CIRCLE-seq
# are both annotated by GRCh37

def mappipe(ifold, ref_file, minlen=20, rclip=0):
    
    ifold = os.path.join(ifold, '')
    ifile = '*.fastq.gz'#  '*.fastq.gz'
    #ref_file = '/data/index/HG19.fasta'
    trim_regex = r'(.*)\/(SRR.+).fastq.gz$'
    pipeline = ruffus.Pipeline('FastqDNaseSeq')
    trim_task = pipeline.collate(tasks.trimmer,
                                 name = 'TrimGalore',
                                 input = ifold+ifile, 
                                 filter = ruffus.regex(trim_regex),
                                 output = r'\1/\2_trimmed.fq.gz',
                                 # extras[0]: minimum length,
                                 # [1]:right end clip size
                                 extras= [[minlen, rclip]])
    trfile = '*_trimmed.fq.gz'
    aln_regex = r'(.*)\/(.*).fq.gz$'
    align_task = pipeline.collate(tasks.bwa_aln,
                                name = 'bwa_aln',
                                input =  ifold+ trfile,
                                filter = ruffus.regex(aln_regex),
                                output = r'\1/\2.sai',
                                extras= [ref_file])
    align_task.follows('TrimGalore')

    ## sai to sam file using bwa samse
    sai_file = '*.sai'
    samse_regex = r'(.*)\/(.*).sai$'
    samse_task = pipeline.collate(tasks.bwa_samse,
                                  name = 'bwa_samse',
                                  input = ifold + sai_file,
                                  filter = ruffus.regex(samse_regex),
                                  output = r'\1/\2.sam',
                                  # extras[0]: fastq required for samse,
                                  # [1]: ref indexed fasta,
                                  # [2]: max multiple mapped reads [Default=3]
                                  extras = [[r'\1/\2.fq.gz', ref_file, 10]])
    samse_task.follows('bwa_aln')

    ## sam to bam using sambamba view
    sam_file = '*.sam'
    tobam_regex = r'(.*)\/(.*).sam$'
    tobam_task = pipeline.collate(tasks.sam_to_bam,
                                  name = 'sam_bam',
                                  input = ifold + sam_file,
                                  filter = ruffus.regex(tobam_regex),
                                  output = r'\1/\2.bam')
    tobam_task.follows('bwa_samse')

    ## sorting bam with sambamba sort
    bam_file = '*trimmed.bam'
    sort_bam_regex = r'(.*)\/(.*).bam$'
    sort_bam_task  = pipeline.collate(tasks.sort_bam,
                                      name = 'sorting_bam',
                                      input = ifold + bam_file,
                                      filter = ruffus.regex(sort_bam_regex),
                                      output = r'\1/\2.sorted.bam')
    sort_bam_task.follows('sam_bam')

    ## bam to bed using bam2bed
    sorted_bam_file = '*trimmed.sorted.bam'
    sorted_bam_regex = r'(.*)\/(.*).sorted.bam$'
    sorted_bam_task = pipeline.collate(tasks.bam2bed,
                                       name = 'bam2bed',
                                       input = ifold + sorted_bam_file,
                                       filter = ruffus.regex(sorted_bam_regex),
                                       output = r'\1/\2.sorted.bed')
    sorted_bam_task.follows('sorting_bam')

    full_pipe = ruffus.Pipeline('Full pipeline',
                                input = ['bam2bed'])

    full_pipe.run()

def quick(ifold):
    
    # sorting bam file
    pipeline = ruffus.Pipeline('BamDNaseSeq')
    bam_file = '*.bam'
    sort_bam_regex = r'(.*)\/(.*).bam$'
    sort_bam_task = pipeline.collate(tasks.sort_bam,
                                     name='sorting_bam',
                                     input= os.path.join(ifold, bam_file),
                                     filter=ruffus.regex(sort_bam_regex),
                                     output=r'\1/\2.sorted.bam')
    ## bam to bed using bam2bed
    sorted_bam_file = '*.sorted.bam'
    sorted_bam_regex = r'(.*)\/(.*).sorted.bam$'
    sorted_bam_task = pipeline.collate(tasks.bam2bed,
                                       name='bam2bed',
                                       input=os.path.join(ifold,sorted_bam_file),
                                       filter=ruffus.regex(sorted_bam_regex),
                                       output=r'\1/\2.sorted.bed')
    sorted_bam_task.follows('sorting_bam')

    full_pipe = ruffus.Pipeline('Full pipeline',
                                input=['bam2bed'])

    full_pipe.run()

if __name__ == '__main__':
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('ERROR: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser()
    
    parser.add_argument('--infolder', '-i', dest='i',
                        help='Input file directory. Fetching all fastq/bam files under given directory.', required=True)
    parser.add_argument('--pipeline', '-p', dest='p',
                        choices=['bam', 'fastq'], required=True,
                        help="Pipeline selection ['bam', 'fastq']")
    parser.add_argument('--reference', '-r', dest='r',
                        help='Reference genome, e.g. HG19')
    
    args = parser.parse_args()
    
    
    fullpath = os.path.dirname(os.path.abspath(args.i))
    if args.p == 'bam':
        print('Pipeline bam...')
        print('Reference genome is not required.')
        if not os.path.exists(fullpath):
            print('Folder does not exist. Creating directory..')
            os.makedirs(fullpath)
        else:
            print('Folder exists...')
        quick(fullpath)
            
    elif args.p == 'fastq':
        print('pipeline mappipe')
        if not args.r:
            parser.error('The --reference argument is required for [fastq] pipeline.\n')
        elif not os.path.exists(args.r):
            parser.error('The --reference file does not exist\n')
        else:
            refpath = os.path.abspath(args.r)
        if not os.path.exists(fullpath):
            print('Folder does not exist. Creating directory..\n%s\n'%fullpath)
            os.makedirs(fullpath)
        else:
            print('Folder exists...')
        mappipe(fullpath, refpath)
    else:
        print('help message')
