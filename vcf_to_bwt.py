import sys
import argparse
import subprocess as sp
import multiprocessing as mp
import os
import logging

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
        self.ref_only = False
        self.m = other.marker_array

class VcfScanCmd:
    def __init__(self, args):
        self.H = args.h
        self.f = args.fasta
        self.S = args.sample
        self.o = args.prefix
        self.c = '21' # TODO: make this an option
        self.stdout = True
        self.ref = False
        self.m = False
        self.vcfs = args.vcf

    def set_ref(self):
        self.ref = True

    def set_m(self):
        self.m = True

    def get_cmd(self):
        cmd = ["./vcf_scan", '-f', self.f, "-c", self.c,  "-o", self.o]
        if self.stdout:
            cmd += ['--stdout']
        if self.ref:
            cmd +=  ['-r']
        else:
            cmd += ["-H", self.H, "-S", self.S]
        if self.m:
            cmd += ['-m']
        cmd += self.vcfs
        return cmd

def vcf_to_fa(args, ref=False):
    cmd_builder = VcfScanCmd(args)
    full_prefix = ""
    if ref:
        cmd_builder.set_ref()
        full_prefix = args.prefix + ".ref"
    else:
        full_prefix = args.full_prefix
    if args.m:
        cmd_builder.set_m()
    log = open(full_prefix + ".log", "w")
    fa_fname = full_prefix + ".fa"
    with open(fa_fname, "w") as fa_fp:
        vcf_scan_proc = sp.run(cmd_builder.get_cmd(), stdout=fa_fp, stderr=log)
    log.close()


def vcf_to_parse(args, ref=False):
    cmd_builder = VcfScanCmd(args)
    full_prefix = ""
    if ref:
        cmd_builder.set_ref()
        full_prefix = args.prefix + ".ref"
    else:
        full_prefix = args.full_prefix
    if args.m:
        cmd_builder.set_m()
    vcf_scan_cmd = cmd_builder.get_cmd()
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

class vcf_to_fa_builder(ref=False):
    def __init__(self, ref=False):
        self.ref = ref
    def __call__(self, args):
        return vcf_to_fa(args, self.ref)

class vcf_to_parse_builder:
    def __init__(self, ref=False):
        self.ref = ref
    def __call__(self, args):
        return vcf_to_parse(args, self.ref)

def vcf_to_bwt_no_merge(args):
    samples = open(args.samples).read().strip().split('\n')
    thread_args = [VcfToXArgs(i, s, h, args) for i, s in enumerate(samples) for h in ['0', '1']]
    all_prefixes = [args.o + ".ref"] + [a.full_prefix for a in thread_args]
    # parsing
    print("generating fastas")
    fasta_pool = mp.Pool(processes=args.threads)
    ref_worker = fasta_pool.apply_async(vcf_to_fa_builder(True), (VcfToXArgs(0, "", 0, args),))
    fasta_pool.imap_unordered(vcf_to_fa_builder(False), thread_args)
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
    # set up logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    log_fp = open(args.o + ".log", 'w')
    log_h = logging.StreamHandler(log_fp)
    log_h.setFormatter(logging.Formatter('[%(asctime)s] %(message)s') )
    stderr_h = logging.StreamHandler(sys.stderr)
    stderr_h.setFormatter(logging.Formatter('[%(asctime)s] %(message)s') )
    logger.addHandler(log_h)
    logger.addHandler(stderr_h)

    # parsing
    logger.info("generating parses from VCF")
    parse_pool = mp.Pool(processes=args.threads)
    ref_proc = parse_pool.apply_async(vcf_to_parse_builder(True), (VcfToXArgs(0, "", 0, args),))
    parse_pool.imap_unordered(vcf_to_parse_builder(False), thread_args)
    ref_proc.get()
    parse_pool.close()
    parse_pool.join()
    logger.info("done generating parses from VCF")
    # merge
    logger.info("merging parses")
    merge_pfp_cmd = ['./merge_pfp', '-s', '--parse-bwt', '--docs', '-o', args.o, '-t', str(args.threads)] + all_prefixes
    logger.info(" ".join(merge_pfp_cmd))
    sp.run(merge_pfp_cmd, stdout=log_fp, stderr=sp.PIPE, check=True)
    logger.info("done merging parses")
    if args.clean:
        logger.info("cleaning files")
        for prefix in all_prefixes:
            clean_parse_files(prefix)
    # construct BWT
    logger.info("constructing BWT")
    pfbwt_cmd = ['./pfbwt-f64', '-s', '--pfbwt-only', '--print-docs', '-o', args.o]
    logger.info(" ".join(pfbwt_cmd))
    sp.run(pfbwt_cmd, stdout=log_fp, stderr=sp.PIPE, check=True)
    logger.info("done constructing BWT")

    log_h.close()
    log_fp.close()


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
    parser.add_argument("--marker_array", "-m", action='store_true', help="cleanup intermediate files as we go")
    args = parser.parse_args()

    vcf_to_bwt(args)
