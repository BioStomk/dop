import os
import subprocess
import math

def plot(x_name, y_name, x_len, y_len, kmer_size, f_file, b_file, out_file, scale=False, color=False):
    cmd = 'gnuplot -e "'
    cmd += 'set x2range [0 : %d];' % x_len
    cmd += 'set yrange  [%d : 0];' % y_len
    cmd += 'set x2tics %d;' % 10000
    cmd += 'set ytics  %d;' % 10000
    cmd += 'set mx2tics 5;'
    cmd += 'set mytics  5;'
    cmd += 'set tics out;'
    cmd += 'set tics scale 1.5, 0.5;'
    cmd += 'set x2tics offset 0, -0.5;'
    cmd += 'set ytics  offset 0.8, 0;'
    cmd += "set x2label '%s';" % x_name.replace('_', '\_')
    cmd += "set ylabel  '%s';" % y_name.replace('_', '\_')
    cmd += 'set x2label offset 0, -0.8;'
    cmd += 'set ylabel  offset 2.0, 0;'
    cmd += "set label 'k-mer size = %d' at screen 0.01, 0.98 ;" % kmer_size
    cmd += 'unset key;'
    cmd += 'unset xtics;'
    if scale:
        x_size = x_len / 100
        y_size = y_len / 100
    else:
        if x_len > y_len:
            x_size = 50 * math.log(x_len, 2)
            y_size = x_size * y_len / x_len * 1.2
        else:
            y_size = 50 * math.log(y_len, 2)
            x_size = y_size * x_len / y_len * 1.2
    cmd += 'set terminal png size %d, %d;' % (x_size, y_size)
    cmd += "set output '%s';" % out_file
    if color:
        cmd += "plot '%s' w p ps 0.15 lc rgb 'red', '%s' w p ps 0.15 lc rgb 'dark-green';" % (f_file, b_file)
    else:
        if not os.path.exists(f_file) or not os.path.exists(b_file):
            sys.exit('Cannot find match files. Exit')
        cat_cmd = 'tail -n +1 %s >> %s' % (b_file, f_file)
        subprocess.check_call(cat_cmd, shell=True)
        cmd += "plot '%s' w p ps 0.15 lc rgb 'black';" % f_file
    cmd += '"'
    subprocess.check_call(cmd, shell=True)
