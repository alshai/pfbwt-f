import sys
import struct

ALE_MASK = 0xF000000000000000
SEQ_MASK = 0x0FFFF00000000000
POS_MASK = 0x00000FFFFFFFFFFF
SEQ_SHIFT = 46
ALE_SHIFT = 60

if __name__ == "__main__":
    buffer = open(sys.argv[1], 'rb').read()
    vals = struct.iter_unpack('Q', buffer) 
    line = []
    i = 0
    marker = (None,None,None)
    rn = [None, None]
    for val in vals:
        val = val[0]
        if val != 2**64-1:
            if i > 1:
                pos = val & POS_MASK
                seq = (val & SEQ_MASK) >> SEQ_SHIFT
                ale = (val & ALE_MASK) >> ALE_SHIFT
                marker = (seq, pos, ale)
            else:
                rn[i] = val
                i += 1
        else:
            for j in range(rn[0], rn[1]+1):
                sys.stdout.write("{} {} {} {}\n".format( j, marker[0], marker[1], marker[2]))
            i = 0