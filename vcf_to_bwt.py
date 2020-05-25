import sys
import argparse
import subprocess as sp
import multiprocessing as mp
import os
import logging

PARSE_EXTS  = ['bwlast', 'bwsai', 'dict', 'docs', 'fa', 'ilist', 'log', 'mai', 'n', 'occ', 'parse']

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
        self.mmap = other.mmap


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
        if args.mmap:
            pfbwt_cmd = pfbwt_cmd[:1] + ['-m'] + pfbwt_cmd[1:]
        pfbwt_proc = sp.run(pfbwt_cmd, check=True, stdout=log, stderr=sp.PIPE)
    else:
        vcf_scan_proc = sp.Popen(vcf_scan_cmd, stdout=sp.PIPE, stderr=log)
        pfbwt_cmd = ['./pfbwt-f64', '--parse-only', '--print-docs', '-s', '-o', full_prefix]
        if args.mmap:
            pfbwt_cmd = pfbwt_cmd[:1] + ['-m'] + pfbwt_cmd[1:]
        pfbwt_proc = sp.run(pfbwt_cmd, stdin=vcf_scan_proc.stdout, stdout=log, stderr=sp.PIPE, check=True)
        vcf_scan_proc.wait()
    log.close()


class vcf_to_fa_builder:
    def __init__(self, ref=False):
        self.ref = ref
    def __call__(self, args):
        return vcf_to_fa(args, self.ref)


class vcf_to_parse_builder:
    def __init__(self, ref=False):
        self.ref = ref
    def __call__(self, args):
        return vcf_to_parse(args, self.ref)



def merge_marker_indexes(args, thread_args, logger, log_fp):
    logger.info("merging marker indexes")
    length = int(open(args.o + ".ref.n").read().strip())
    cmd = ['./merge_marker_indexes', str(length), args.o + ".mai"] + [args.o + ".ref.mai"] + [a.full_prefix + ".mai" for a in thread_args]
    logger.info(" ".join(cmd))
    sp.run(cmd, check=True, stdout=log_fp, stderr=sp.PIPE)
    logger.info("merged marker indexes")


def vcf_to_bwt(args):
    if not args.samples:
        sys.stderr.write("no support yet for specification of zero samples\n")
        exit(1)
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
    if args.no_merge:
        logger.info("generating fasta files from VCF")
        fasta_pool = mp.Pool(processes=args.threads)
        ref_proc = fasta_pool.apply_async(vcf_to_fa_builder(True), (VcfToXArgs(0, "", 0, args),))
        fasta_pool.imap_unordered(vcf_to_fa_builder(False), thread_args)
        ref_proc.get()
        fasta_pool.close()
        fasta_pool.join()
        logger.info("done generating fastas from VCF")
        # merging
        if args.marker_array:
            merge_marker_indexes(args, thread_args, logger, log_fp)
        logger.info("generating parse (directly from fasta files)")
        parse_cmd = ['./pfbwt-f64', '--parse-only', '-s', '--print-docs', '-o', args.o]
        cat_cmd = ['cat'] + [args.o + ".ref.fa"] + [a.full_prefix + ".fa" for a in thread_args]
        cat_proc = sp.Popen(cat_cmd, stdout=sp.PIPE, stderr=log_fp)
        logger.info(" ".join(cat_cmd) + " | " + " ".join(parse_cmd))
        parse_proc = sp.run(parse_cmd, stdin=cat_proc.stdout, stderr=sp.STDOUT, stdout=log_fp)
        cat_proc.wait()
        logger.info("generated parse")
    else:
        logger.info("generating parses from VCF")
        parse_pool = mp.Pool(processes=args.threads)
        ref_proc = parse_pool.apply_async(vcf_to_parse_builder(True), (VcfToXArgs(0, "", 0, args),))
        parse_pool.imap_unordered(vcf_to_parse_builder(False), thread_args)
        ref_proc.get()
        parse_pool.close()
        parse_pool.join()
        logger.info("done generating parses from VCF")
        # merging
        if args.marker_array:
            merge_marker_indexes(args, thread_args, logger, log_fp)
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
    pfbwt_cmd = ['./pfbwt-f64', '--stdout', 'sa', '-s', '--pfbwt-only', '--print-docs', '-o', args.o]
    if args.mmap:
        pfbwt_cmd = pfbwt_cmd[:1] + ['-m'] + pfbwt_cmd[1:]
    logger.info(" ".join(pfbwt_cmd))
    pfbwt_proc = sp.Popen(pfbwt_cmd, stdout=sp.PIPE, stderr=log_fp)
    if args.marker_array:
        logger.info("also constructing marker array from marker index")
        # marker_array_cmd = ['./marker_index_to_array', args.o + ".mai", '-', args.o + ".ma"]
        marker_array_cmd = ['./marker_index_to_array', "-o", args.o + ".ma", args.o + ".mai", '-']
        if args.mmap:
            marker_array_cmd = marker_array_cmd[:1] + ['-m'] + marker_array_cmd[1:]
        logger.info(" ".join(marker_array_cmd))
        marker_array_proc = sp.run(marker_array_cmd, stdin=pfbwt_proc.stdout, stderr=sp.STDOUT, stdout=log_fp, check=True)
        pfbwt_proc.wait()
    else:
        sa_fp = open(args.o + ".sa", "wb")
        sa_save_proc = sp.run(['cat'], stdin=pfbwt_proc.stdout, stdout=sa_fp, stderr=log_fp, check=True)
        pfbwt_proc.wait()
        sa_fp.close()
    logger.info("done constructing BWT")
    if args.clean and not args.keep_parse:
        logger.info("cleaning final parse files")
        clean_parse_files(args.o)
    log_h.close()
    log_fp.close()


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
    parser.add_argument("--keep_parse", action='store_true', help="keeps the final parse (for when --clean is used)")
    parser.add_argument("--mmap", '-M', action='store_true', help="tell pfbwt-f64 to use mmap (use this for very large files)")
    args = parser.parse_args()

    vcf_to_bwt(args)
