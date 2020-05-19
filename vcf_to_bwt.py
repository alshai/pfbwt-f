import sys
import argparse
import subprocess as sp
import multiprocessing as mp
import os

PARSE_EXTS  = ['bwlast', 'bwsai', 'dict', 'docs', 'fa', 'ilist', 'log', 'mai',
               'n', 'occ', 'parse']

def clean_parse_files(prefix):
    for ext in PARSE_EXTS:
        try:
            os.remove("{}.{}".format(prefix, ext))
        except FileNotFoundError:
            pass

class VcfToXArgs:

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


def vcf_to_fa(args):
    vcf_scan_cmd = ["./vcf_scan", "--stdout",
                    '-f', args.fasta,
                    "-c", '21',
                    "-H", args.h,
                    "-S", args.sample,
                    "-o", args.prefix] + args.vcf
    log = open(args.full_prefix + ".log", "w")
    fa_fname = args.full_prefix + ".fa"
    with open(fa_fname, "w") as fa_fp:
        vcf_scan_proc = sp.run(vcf_scan_cmd, stdout=fa_fp, stderr=log)
    log.close()


def vcf_to_fa_ref(args):
    full_prefix = args.prefix + ".ref"
    vcf_scan_cmd = ["./vcf_scan", "--stdout", "-r",
                    '-f', args.fasta,
                    "-c", '21',
                    "-o", args.prefix] + args.vcf
    log = open(full_prefix + ".log", "w")
    fa_fname = full_prefix + ".fa"
    with open(fa_fname, "w") as fa_fp:
        vcf_scan_proc = sp.run(vcf_scan_cmd, stdout=fa_fp, stderr=log)
    log.close()


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


def vcf_to_bwt_no_merge(args):
    samples = open(args.samples).read().strip().split('\n')
    thread_args = [VcfToXArgs(i, s, h, args) for i, s in enumerate(samples) for h in ['0', '1']]
    all_prefixes = [args.o + ".ref"] + [a.full_prefix for a in thread_args]
    # parsing
    print("generating fastas")
    fasta_pool = mp.Pool(processes=args.threads)
    ref_worker = fasta_pool.apply_async(vcf_to_fa_ref, (VcfToXArgs(0, "", 0, args),))
    fasta_pool.imap_unordered(vcf_to_fa, thread_args)
    ref_worker.get()
    fasta_pool.close()
    fasta_pool.join()
    print("running pfbwt")
    log = open(args.o + ".log", 'w')
    cat_cmd = ['cat'] + [p + ".fa" for p in all_prefixes]
    cat_proc = sp.Popen(cat_cmd, stdout=sp.PIPE)
    all_pfbwt_cmd = ['./pfbwt-f64', '-s', '--print-docs', '-o', args.o]
    print(' '.join(cat_cmd) + " | " + ' '.join(all_pfbwt_cmd))
    sp.run(all_pfbwt_cmd, check=True, stdin=cat_proc.stdout, stderr=log)
    cat_proc.wait()


def vcf_to_bwt_w_merge(args):
    samples = open(args.samples).read().strip().split('\n')
    thread_args = [VcfToXArgs(i, s, h, args) for i, s in enumerate(samples) for h in ['0', '1']]
    all_prefixes = [args.o + ".ref"] + [a.full_prefix for a in thread_args]
    # parsing
    print("generating parses")
    parse_pool = mp.Pool(processes=args.threads)
    ref_proc = parse_pool.apply_async(vcf_to_parse_ref, (VcfToXArgs(0, "", 0, args),))
    parse_pool.imap_unordered(vcf_to_parse, thread_args)
    ref_proc.get()
    parse_pool.close()
    parse_pool.join()
    log = open(args.o + ".log", 'w')
    # merge
    merge_pfp_cmd = ['./merge_pfp', '-s', '--parse-bwt', '--docs', '-o', args.o, '-t', str(args.threads)] + all_prefixes
    print(' '.join(merge_pfp_cmd))
    sp.run(merge_pfp_cmd, stdout=log, stderr=sp.PIPE, check=True)
    if args.clean:
        for prefix in all_prefixes:
            clean_parse_files(prefix)
    # construct BWT
    pfbwt_cmd = ['./pfbwt-f64', '-s', '--pfbwt-only', '--print-docs', '-o', args.o]
    print(' '.join(pfbwt_cmd))
    sp.run(pfbwt_cmd, stdout=log, stderr=sp.PIPE, check=True)
    log.close()


def vcf_to_bwt(args):
    if not args.samples:
        sys.stderr.write("no support yet for specification of zero samples\n")
        exit(1)
    if args.no_merge:
        vcf_to_bwt_no_merge(args)
    else:
        vcf_to_bwt_w_merge(args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help='reference fasta file')
    parser.add_argument("vcf", nargs='+', help='vcf files containing haplotype panel')
    parser.add_argument("--samples", "-S", help="file containing new-line delimited desired samples in vcf. defaults to all")
    parser.add_argument("--threads", "-t", type=int, default=int(1), help="number of threads (default: 1)")
    parser.add_argument("--save_fasta", "-f", action='store_true', help="store fasta sequences generated from VCFs")
    parser.add_argument("-o", default="out", help="output prefix")
    parser.add_argument("--no_merge", action='store_true', help="generate a BWT from a non-merged parse of text collection")
    parser.add_argument("--clean", action='store_true', help="cleanup intermediate files as we go")
    args = parser.parse_args()

    vcf_to_bwt(args)
