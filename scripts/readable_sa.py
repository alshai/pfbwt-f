import sys
import struct

if __name__ == "__main__":
    buffer = open(sys.argv[1], 'rb').read()
    vals = struct.iter_unpack('Q', buffer) 
    line = []
    print("\n".join(list(map(str, [v[0] for v in vals]))))