from subprocess import check_call
import shlex
import os
import gzip
import shutil
from tempfile import TemporaryDirectory
import csv

def trimmer(fastq_path, out_path, para_list):
    if para_list[1]>0:
        rclipCmd = '--three_prime_clip_R1 %(rclip)i'
    else:
        rclipCmd =''
    cmd1 = 'trim_galore --fastqc --illumina --length %(len)i --gzip '
    cmd2 = '-o %(direc)s %(fastq)s'
    cmd = cmd1 + rclipCmd + cmd2
    direc = out_path.rsplit(os.sep,1)[0]
    tdict = {'fastq': fastq_path[0],
             'direc': direc,
             'len': para_list[0],
             'rclip': para_list[1]}
    print(cmd % tdict)
    check_call(shlex.split(cmd%tdict))

def bwa_aln(fastq_path, sai_path, ref_path):
    cmd = 'bwa aln -l %(seedlen)i -t ' \
          '%(threads)i -f %(out_path)s ' \
          '%(ref)s %(fasta_path)s'

    tdict = {'seedlen': 17,
             'threads' : 7,
             'out_path': sai_path,
             'ref': ref_path,
             'fasta_path': fastq_path[0]}
    print(cmd % tdict)
    check_call(shlex.split(cmd % tdict))

def bwa_samse(sai_path, sam_path, para_list):
    cmd = 'bwa samse -n %(mmm)i -f %(sam)s %(ref)s %(sai)s %(fastq)s'
    try:
        fastq_path = para_list[0]
        ref_path = para_list[1]
    except IndexError:
        print('[fastq_path] and [ref_path] are required')

    try:
        mul_map_max = para_list[2]
    except IndexError:
        mul_map_max =3
        print('Only reads with multiple mapping less than %i will be listed ' % mul_map_max)
        pass

    tdict = {'ref': ref_path,
             'sai': sai_path[0],
             'fastq': fastq_path,
             'sam': sam_path,
             'mmm': mul_map_max}

    check_call(shlex.split(cmd %tdict))

def sam_to_bam(sam_file, bam_file):
    cmd = 'sambamba view -S -f bam -l 9 -t 7 -o %s %s' % (bam_file,sam_file[0])
    check_call(shlex.split(cmd))

def sort_bam(bam_file, sorted_file):
    cmd = 'sambamba sort -t 7 -o %s %s' %(sorted_file, bam_file[0])
    check_call(shlex.split(cmd))

def bam2bed(bam_file, bed_file):
    cmd = 'bedtools bamtobed -i %s' %(bam_file)
    with open(bed_file, 'w') as handle:
        check_call(shlex.split(cmd), stdout=handle)

def test_bwa_aln():
    fastq = '/data/DNaseSeq/U2OS/alignment/SRR991/SRR991.trim.1000.fq'
    ref = '/data/index/HG19.fasta'
    sai = fastq.replace('.fq', '.sai')
    bwa_aln(fastq, ref, sai)

def beddepth(region_file, bam_file, depth_file):
    cmd = 'sambamba depth base -t 7 --min-coverage=0 -L %s %s' % (region_file, bam_file)
    with open(depth_file, 'w') as handle:
        check_call(shlex.split(cmd), stdout=handle)

