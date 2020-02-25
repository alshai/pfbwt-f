import sys
import argparse
from scipy import sparse
import numpy as np

def read_marker_file(f, nbytes=8, endian="little"):
    # get size of file
    f.seek(0,2)
    fsize = f.tell()
    if (fsize % nbytes) or (fsize % 3):
        sys.stderr.write("invalid file ({} bytes)!\n".format(fsize))
        exit(1)
    else:
        print("{} records".format(fsize / nbytes / 3))
    f.seek(0,0)
    marr = np.ndarray((int(fsize / nbytes / 3), 3), dtype=np.int64)
    bytestr = f.read(nbytes)
    i=0
    while bytestr:
        x = int.from_bytes(bytestr, endian)
        if (i % 3 < 2):
            marr[int(i/3), i % 3] = x
        else:
            marr[int(i/3), i % 3] = 0
        bytestr = f.read(nbytes)
        i += 1
    return sparse.dok_matrix(sparse.coo_matrix((marr[:,1], (marr[:,0], marr[:,2]))))


def generate_marker_array(args):
    nbytes = args.bytes
    endian = args.endian
    with open(args.marker_info, "rb") as f:
        mar = read_marker_file(f, nbytes, endian)
    of = open(args.o + ".ma", "wb")
    if args.sa == "-":
        f = sys.stdin.buffer
    else:
        f = open(args.sa, "rb")
    bytestr = f.read(nbytes)
    while bytestr:
        s = int.from_bytes(bytestr, endian)
        if s >= mar.shape[0]:
            of.write(int(0).to_bytes(nbytes, endian))
        else:
            marker = mar[s,0]
            of.write(int(marker).to_bytes(nbytes, endian))
        bytestr = f.read(nbytes)
    f.close()
    of.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("sa")
    parser.add_argument("marker_info")
    parser.add_argument("-o", default="haplotypes")
    parser.add_argument("--bytes", type=int, default=8)
    parser.add_argument("--endian", default="little")
    args = parser.parse_args()

    generate_marker_array(args)
