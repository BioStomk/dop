import subprocess
import math

def plot(x_name, y_name, x_len, y_len, data, out_file, scale=False):
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
    cmd += 'unset key;'
    cmd += 'unset xtics;'
    if scale:
        if x_len > y_len:
            x_size = 50 * math.log(x_len, 2)
            y_size = x_size * y_len / x_len * 1.2
        else:
            y_size = 50 * math.log(y_len, 2)
            x_size = y_size * x_len / y_len * 1.2
    else:
        x_size = x_len / 100
        y_size = y_len / 100
    print x_size, y_size
    cmd += 'set terminal png size %d, %d;' % (x_size, y_size)
    cmd += "set output '%s';" % out_file
    cmd += "plot '%s' w p ps 0.15 lc rgb 'black'" % data
    cmd += '"'
    subprocess.check_call(cmd, shell=True)
