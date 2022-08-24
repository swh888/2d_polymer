import math
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Calculate block averages.")

parser.add_argument("-i", "--inputfile", default="input.in", metavar='', help="input file name")
#parser.add_argument("-o", "--outputfile", default="output.out", metavar='', help="output file name")
#parser.add_argument("-cx", "--column", default=0, type=int, metavar='', help="input column, col 0 is time")
args = parser.parse_args()

inputfile=args.inputfile
#outputfile=args.outputfile
#cx=args.column
f=open(inputfile)
g=open("s-data.txt","w")
h=open("sd-total.txt","w")

data = np.loadtxt(inputfile, skiprows=1)
data = np.log(data)

blksize = []
l = len(data)
for j in range(1, math.floor(np.log2(l))):
    blksize.append(int(math.pow(j, 2)))

total_avg = np.mean(data)
total_var = np.sum(np.square(data - total_avg))/(l - 1)

for b in blksize:
    nblk = math.floor(len(data)/b)
    blk_var = 0.0
    for n in range(0,nblk):
        blkavg = 0.0
        min_ind = n * b
        max_ind = (n + 1) * b
        for ind in range(min_ind, max_ind):
            blkavg += data[ind]/b
        blk_var += (blkavg - total_avg) * (blkavg - total_avg) / (nblk - 1)
    g.write("{:f} {:f}\n".format(1/b, b * blk_var / total_var))
g.close()

# fit the linear region to calculate corrected uncertainty
s_data = np.loadtxt('s-data.txt')

s_data = s_data[np.where(s_data[:,0] < 0.05)[0]]

xavg = np.mean(s_data[:,0])
yavg = np.mean(s_data[:,1])
ssxx = np.sum(np.square(s_data[:,0] - xavg))
ssyy = np.sum(np.square(s_data[:,1] - yavg))
ssxy = np.sum((s_data[:,0] - xavg) * (s_data[:,1] - yavg))

s_slope = ssxy/ssxx
s_intercept = yavg - s_slope * xavg

# use the intercept to calculate the standard error of the average
point_error = 1.96 * np.sqrt(total_var * s_intercept / len(data))

h.write(f"{point_error}")
h.close()



