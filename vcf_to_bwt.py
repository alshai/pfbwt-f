import sys
import argparse
import subprocess as sp
import multiprocessing as mp

class VcfToParseArgs:

    def __init__(self, i, sample, h, other):
        self.i = i
        self.sample = sample
        # self.o = other
        self.fasta = other.fasta
        self.vcf = other.vcf
        self.h = h
        self.prefix = "{}".format(other.o)
        self.full_prefix = '{}.{}.{}'.format(self.prefix, self.sample, self.h)


# todo: write to log also
def vcf_to_parse(args):
    print(args.h)
    vcf_scan_cmd = ["./vcf_scan", "--stdout", 
            '-f', args.fasta, 
            "-c", '21', 
            "-H", args.h, 
            "-S", args.sample,
            "-o", args.prefix,
            args.vcf]
    vcf_scan_proc = sp.Popen(vcf_scan_cmd, stdout=sp.PIPE)
    pfbwt_cmd = ['./pfbwt-f64', '--parse-only', '--print-docs', '-s', 
                 '-o', args.full_prefix]
    pfbwt_proc = sp.run(pfbwt_cmd, stdin=vcf_scan_proc.stdout, check=True)
    vcf_scan_proc.wait()

def vcf_to_bwt(args):
    if not args.samples:
        sys.stderr.write("no support yet for specification of zero samples\n")
        exit(1)
    samples = open(args.samples).read().strip().split('\n')

    thread_args = [VcfToParseArgs(i, s, h, args) for i, s in enumerate(samples) for h in ['0', '1']]
    parse_pool = mp.Pool(processes=args.threads)
    parse_pool.map(vcf_to_parse, thread_args)
    parse_pool.close()
    parse_pool.join()

    merge_pfp_cmd = ['./merge_pfp', '--docs', '-o', args.o, '-t', str(args.threads)] + [a.full_prefix for a in thread_args]
    print(' '.join(merge_pfp_cmd))
    sp.run(merge_pfp_cmd, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta")
    parser.add_argument("vcf")
    parser.add_argument("--samples", "-S", help="file containing new-line delimited desired samples in vcf. defaults to all")
    parser.add_argument("--threads", "-t", type=int, default=int(1), help="number of threads (default: 1)")
    parser.add_argument("-o", default="out", help="output prefix")
    args = parser.parse_args()

    vcf_to_bwt(args)
