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
        self.save_fasta = other.save_fasta

def vcf_to_parse(args):
    vcf_scan_cmd = ["./vcf_scan", "--stdout",
                    '-f', args.fasta,
                    "-c", '21',
                    "-H", args.h,
                    "-S", args.sample,
                    "-o", args.prefix] + args.vcf
    log = open(args.full_prefix + ".log", "w")
    if args.save_fasta:
        fa_fname = args.full_prefix + ".fa"
        with open(fa_fname, "w") as fa_fp:
            vcf_scan_proc = sp.run(vcf_scan_cmd, stdout=fa_fp, stderr=log)
        pfbwt_cmd = ['./pfbwt-f64', '--parse-only', '--print-docs', '-s', '-o', args.full_prefix, fa_fname]
        pfbwt_proc = sp.run(pfbwt_cmd, check=True, stdout=log, stderr=sp.PIPE)
    else:
        vcf_scan_proc = sp.Popen(vcf_scan_cmd, stdout=sp.PIPE, stderr=log)
        pfbwt_cmd = ['./pfbwt-f64', '--parse-only', '--print-docs', '-s', '-o', args.full_prefix]
        pfbwt_proc = sp.run(pfbwt_cmd, stdin=vcf_scan_proc.stdout, stdout=log, stderr=sp.PIPE, check=True)
        vcf_scan_proc.wait()
    log.close()


def vcf_to_parse_ref(args):
    full_prefix = args.prefix + ".ref"
    vcf_scan_cmd = ["./vcf_scan", "--stdout", "-r",
                    '-f', args.fasta,
                    "-c", '21',
                    "-o", args.prefix] + args.vcf
    log = open(full_prefix + ".log", "w")
    if args.save_fasta:
        fa_fname = full_prefix + ".fa"
        with open(fa_fname, "w") as fa_fp:
            vcf_scan_proc = sp.run(vcf_scan_cmd, stdout=fa_fp, stderr=log)
        pfbwt_cmd = ['./pfbwt-f64', '--parse-only', '--print-docs', '-s', '-o', full_prefix, fa_fname]
        pfbwt_proc = sp.run(pfbwt_cmd, check=True, stdout=log, stderr=sp.PIPE)
    else:
        vcf_scan_proc = sp.Popen(vcf_scan_cmd, stdout=sp.PIPE, stderr=log)
        pfbwt_cmd = ['./pfbwt-f64', '--parse-only', '--print-docs', '-s', '-o', full_prefix]
        pfbwt_proc = sp.run(pfbwt_cmd, stdin=vcf_scan_proc.stdout, stdout=log, stderr=sp.PIPE, check=True)
        vcf_scan_proc.wait()
    log.close()


def vcf_to_bwt(args):
    if not args.samples:
        sys.stderr.write("no support yet for specification of zero samples\n")
        exit(1)
    samples = open(args.samples).read().strip().split('\n')
    thread_args = [VcfToParseArgs(i, s, h, args) for i, s in enumerate(samples) for h in ['0', '1']]
    all_prefixes = [args.o + ".ref"] + [a.full_prefix for a in thread_args]
    # parsing
    print("generating parses")
    parse_pool = mp.Pool(processes=args.threads)
    ref_parse = parse_pool.apply_async(vcf_to_parse_ref, (VcfToParseArgs(0, "", 0, args),))
    parse_pool.imap_unordered(vcf_to_parse, thread_args)
    ref_parse.get()
    parse_pool.close()
    parse_pool.join()
    log = open(args.o + ".log", 'w')
    if args.no_merge: # note: triggers args.save_fasta
        args.save_fasta = True
        with open(args.o + ".all.fa", "w") as fa_fp:
            sp.run(['cat'] + [p + ".fa" for p in all_prefixes], check=True, stdout=fa_fp)
        all_pfbwt_cmd = ['./pfbwt-f64', '-s', '--print-docs', '-o', args.o + ".all", args.o + ".all.fa"]
        print(' '.join(all_pfbwt_cmd))
        sp.run(all_pfbwt_cmd, check=True, stderr=log)
    else: # merging and BWT construction
        merge_pfp_cmd = ['./merge_pfp', '-s', '--parse-bwt', '--docs', '-o', args.o, '-t', str(args.threads)] + all_prefixes
        print(' '.join(merge_pfp_cmd))
        sp.run(merge_pfp_cmd, stdout=log, stderr=sp.PIPE, check=True)
        # TODO: cleanup of unneeded files here
        pfbwt_cmd = ['./pfbwt-f64', '-s', '--pfbwt-only', '--print-docs', '-o', args.o]
        print(' '.join(pfbwt_cmd))
        sp.run(pfbwt_cmd, stdout=log, stderr=sp.PIPE, check=True)
        log.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help='reference fasta file')
    parser.add_argument("vcf", nargs='+', help='vcf files containing haplotype panel')
    parser.add_argument("--samples", "-S", help="file containing new-line delimited desired samples in vcf. defaults to all")
    parser.add_argument("--threads", "-t", type=int, default=int(1), help="number of threads (default: 1)")
    parser.add_argument("--save_fasta", "-f", action='store_true', help="store fasta sequences generated from VCFs")
    parser.add_argument("-o", default="out", help="output prefix")
    parser.add_argument("--no_merge", action='store_true', help="generate a BWT from a non-merged parse of text collection")
    args = parser.parse_args()

    vcf_to_bwt(args)
