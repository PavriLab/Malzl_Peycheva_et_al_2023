import argparse as ap
import pysam as ps
import subprocess as sp
import ntpath, os
import logging

def get_tags(read):
    tags = {}
    for tag in ['AS', 'XM', 'XO', 'XG', 'NM']:
        try:
            tags[tag] = read.get_tag(tag)

        except KeyError:
            tags[tag] = 100000

    return tags


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser('this is a python implementation of the HARP algorithm used in doi:10.1101/gr.232561.117')
parser.add_argument('-f', '--fastq', nargs ='+',
                    help = 'one or more fastqs to process')
parser.add_argument('-bti1', '--bowtie2index1', required = True,
                    help = 'bowtie2 index for genome1')
parser.add_argument('-bti2', '--bowtie2index2', required = True,
                    help = 'bowtie2 index for genome2')
parser.add_argument('-mapq', '--minMapQual', default = 20, type = int,
                    help = 'minimum mapping quality for an alignment')
parser.add_argument('-t', '--threads', default = 1, type = int,
                    help = 'number of threads to use for read mapping with bowtie2')
parser.add_argument('-o', '--outdir', default = '.',
                    help = 'directory to write the outputfiles to')
args = parser.parse_args()

for fq in args.fastq:
    fqname = ntpath.basename(fq)
    # mapping reads to reference genomes
    samfiles = {}
    for btindex in [args.bowtie2index1, args.bowtie2index2]:
        logging.info('mapping %s against %s using %i threads' % (fq, ntpath.basename(btindex), args.threads))

        sam = os.path.join(args.outdir, '_'.join([fqname.split('.')[0], ntpath.basename(btindex)]) + '.sam')
        samfiles[ntpath.basename(btindex)] = sam
        cmdstring = 'bowtie2 -x {0} -p {1} --end-to-end --very-sensitive --reorder -U {2} -S {3}'.format(btindex,
                                                                                                         args.threads,
                                                                                                         fq,
                                                                                                         sam)
        logging.info(cmdstring)
        align = sp.Popen(cmdstring, shell = True)
        align.wait()

    # generating outputfile names
    genome1 = open(samfiles[ntpath.basename(args.bowtie2index1)].split('.')[0] + '.fq', 'w')
    genome2 = open(samfiles[ntpath.basename(args.bowtie2index2)].split('.')[0] + '.fq', 'w')
    unmapped = open('_'.join([fq.split('.')[0], 'unmapped']) + '.fq', 'w')
    ambiguous = open('_'.join([fq.split('.')[0], 'ambiguous']) + '.fq', 'w')


    # parsing samfiles
    align1 = ps.AlignmentFile(samfiles[ntpath.basename(args.bowtie2index1)], 'r')
    align2 = ps.AlignmentFile(samfiles[ntpath.basename(args.bowtie2index2)], 'r')

    readcount = {k: 0 for k in ['g1', 'g2', 'umap', 'ambi']}
    mapqt = args.minMapQual
    while True:
        try:
            r1, r2 = align1.__next__(), align2.__next__()

        except StopIteration:
            break

        r1name, r2name = r1.query_name, r2.query_name
        r1tags, r2tags = get_tags(r1), get_tags(r2)

        #logging.info('{0} {1} {2}'.format(str(r1tags), str(r1.is_unmapped), r1.mapq))
        #logging.info('{0} {1} {2}'.format(str(r2tags), str(r2.is_unmapped), r2.mapq))
        # if the readnames do not match we raise an error
        if r1name != r2name:
            raise Exception('samfiles seem to be reordered. make sure they are in the same order as the input fastq')

        # if the read does not map to either genome we write to unmapped
        elif (r1.is_unmapped and r2.is_unmapped):
            readcount['umap'] += 1
            unmapped.write('@{0}\n{1}\n+{0}\n{2}\n'.format(r1name, r1.seq, r1.qqual))

        # if r is unmapped on genome 2 but mapped to genome 1 with decent quality
        # we write to genome1
        elif (r2.is_unmapped and r1.mapq >= mapqt):
            readcount['g1'] += 1
            genome1.write('@{0}\n{1}\n+{0}\n{2}\n'.format(r1name, r1.seq, r1.qqual))

        # same but with genome 2
        elif (r1.is_unmapped and r2.mapq >= mapqt):
            readcount['g2'] += 1
            genome2.write('@{0}\n{1}\n+{0}\n{2}\n'.format(r1name, r1.seq, r1.qqual))

        # if the read maps to both genomes we assign it to the genome with
        # the best match
        elif (r1tags['NM'] < r2tags['NM'] and r1tags['XM'] < r2tags['XM'] and
              r1.mapq >= mapqt and r2.mapq >= mapqt):
            readcount['g1'] += 1
            genome1.write('@{0}\n{1}\n+{0}\n{2}\n'.format(r1name, r1.seq, r1.qqual))

        elif (r2tags['NM'] < r1tags['NM'] and r2tags['XM'] < r1tags['XM'] and
              r1.mapq >= mapqt and r2.mapq >= mapqt):
            readcount['g2'] += 1
            genome2.write('@{0}\n{1}\n+{0}\n{2}\n'.format(r1name, r1.seq, r1.qqual))

        # if the read maps equally good to both genomes it is ambiguous
        elif (r1tags['NM'] == r2tags['NM'] and r1tags['XM'] == r2tags['XM']):
            readcount['ambi'] += 1
            ambiguous.write('@{0}\n{1}\n+{0}\n{2}\n'.format(r1name, r1.seq, r1.qqual))

        # if it does not fulfil any of the above it is unmapped
        else:
            readcount['umap'] += 1
            unmapped.write('@{0}\n{1}\n+{0}\n{2}\n'.format(r1name, r1.seq, r1.qqual))

    logging.info('summary statistics')
    for k, c in readcount.items():
        logging.info('{0}:\t{1}'.format(k, c))