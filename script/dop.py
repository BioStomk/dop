#!/usr//bin/env python

import sys
import os.path
import argparse
import subprocess
import plot

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description= "Draw dot plot by calculating k-mer match between two DNA sequences")

parser.add_argument('seq_files', nargs='+',
                    help='FASTA files of sequences to be compared. If only one file is specified, \
                          the sequence is compared to itself')

parser.add_argument('-d', '--outdir', default='.',
                    help='output directory')

parser.add_argument('-k', '--kmer_size', type=int, default=20,
                    help='k-mer size for matching')

parser.add_argument('-s', '--scale', action='store_true',
                    help='Scale output image size (to look into details of long sequence comparison)')

parser.add_argument('-c', '--color', action='store_true',
                    help='Draw dots with two colors to distinguish forward and reverse matches')

parser.add_argument('--no_labels', action='store_true',
                    help='Do not draw labels')

parser.add_argument('--no_tics', action='store_true',
                    help='Do not draw tics')

parser.add_argument('-p', '--plain', action='store_true',
                    help='Alias for --no_labels --no_tics')

parser.add_argument('-m', '--store_matches', action='store_true',
                    help='Store match files')

parser.add_argument('-f', '--force', dest='recalculate_match', action='store_true',
                    help='Force to calculate matches even if the match files exist')


args = parser.parse_args()

if len(args.seq_files) == 1:
    q_seq_file = t_seq_file = args.seq_files[0]
    q_seq_filename = t_seq_filename = os.path.basename(q_seq_file)
    out_filename = q_seq_filename + '.' + str(args.kmer_size) + '.png'
elif len(args.seq_files) == 2:
    q_seq_file, t_seq_file = args.seq_files
    q_seq_filename = os.path.basename(q_seq_file)
    t_seq_filename = os.path.basename(t_seq_file)
    out_filename = q_seq_filename + '__' + t_seq_filename + '.' + str(args.kmer_size) + '.png'
else:
    sys.exit('Error: Specify one or two FASTA files')

if args.plain:
    args.no_labels = True
    args.no_tics   = True

if not os.path.exists(args.outdir):
    os.mkdir(args.outdir)

out_file = os.path.join(args.outdir, out_filename)

script_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(os.path.dirname(script_dir), 'src')

f_match_file = '%s__%s.match.%d.forward' % (t_seq_filename, q_seq_filename, args.kmer_size)
b_match_file = '%s__%s.match.%d.reverse' % (t_seq_filename, q_seq_filename, args.kmer_size)

if not (os.path.exists(f_match_file) and os.path.exists(b_match_file)) or args.recalculate_match:
    cmd = '%s/comptool search %s %s -k %d' % (src_dir, t_seq_file, q_seq_file, args.kmer_size)
    subprocess.check_call(cmd, shell=True)

def get_seq_len(seq_file):
    with open(seq_file, 'r') as f:
        return sum([len(l.strip('\n')) for l in f.readlines()[1:]])

q_seq_len = get_seq_len(q_seq_file)
t_seq_len = get_seq_len(t_seq_file)
ratio = float(q_seq_len) / t_seq_len

plot.plot(q_seq_filename, t_seq_filename, q_seq_len, t_seq_len, args.kmer_size, f_match_file, b_match_file, out_file, not args.no_tics, not args.no_labels, args.scale, args.color)

if not args.store_matches:
    os.remove(f_match_file)
    os.remove(b_match_file)
