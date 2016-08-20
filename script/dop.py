#!/usr//bin/env python

import argparse
import subprocess
import os.path
import plot

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description= "Draw dot plot by calculating k-mer match between two DNA sequences")

parser.add_argument('-d', '--outdir', default='.',
                    help='output directory')

parser.add_argument('-k', '--kmer_size', type=int, default=20,
                    help='k-mer size for matching')

parser.add_argument('-s', '--scale', action='store_true',
                    help='Scale output image size (to look into details of long sequence comparison)')

parser.add_argument('seqs', nargs='+',
                    help='FASTA files of sequences to be compared. If only one file is specified, \
                          the sequence is compared to itself')

args = parser.parse_args()

if len(args.seqs) == 1:
    q_seq_file = t_seq_file = args.seqs[0]
    q_seq_filename = t_seq_filename = os.path.basename(q_seq_file)
    out_filename = q_seq_filename + '.' + str(args.kmer_size) + '.png'
elif len(args.seqs) == 2:
    q_seq_file, t_seq_file = args.seqs
    q_seq_filename = os.path.basename(q_seq_file)
    t_seq_filename = os.path.basename(t_seq_file)
    out_filename = q_seq_filename + '__' + t_seq_filename + '.' + str(args.kmer_size) + '.png'
else:
    sys.exit('Error: Specify one or two FASTA files')

if not os.path.exists(args.outdir):
    os.mkdir(args.outdir)

out_file = os.path.join(args.outdir, out_filename)

script_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(os.path.dirname(script_dir), 'src')

cmd = '%s/comptool align %s %s -s -k %d' % (src_dir, t_seq_file, q_seq_file, args.kmer_size)
subprocess.check_call(cmd, shell=True)

f_match_file = 'alignments-forward-startpos_%s_%s.tsv'  % (t_seq_filename, q_seq_filename)
b_match_file = 'alignments-backward-startpos_%s_%s.tsv' % (t_seq_filename, q_seq_filename)
if not os.path.exists(f_match_file) or not os.path.exists(b_match_file):
    sys.exit('Cannot find match files. Exit')
cmd = 'tail -n +1 %s >> %s' % (b_match_file, f_match_file)
subprocess.check_call(cmd, shell=True)

def get_seq_len(seq_file):
    with open(seq_file, 'r') as f:
        return sum([len(l.strip('\n')) for l in f.readlines()[1:]])

q_seq_len = get_seq_len(q_seq_file)
t_seq_len = get_seq_len(t_seq_file)
ratio = float(q_seq_len) / t_seq_len

plot.plot(q_seq_filename, t_seq_filename, q_seq_len, t_seq_len, args.kmer_size, f_match_file, out_file, args.scale)

os.remove(f_match_file)
os.remove(b_match_file)
