#!/bin/env python3

import sys
import argparse
import subprocess as sp
import multiprocessing as mp
import os
import logging
import pysam

PARSE_EXTS  = ['bwlast', 'bwsai', 'dict', 'docs', 'fa', 'ilist', 'log', 'mps', 'n', 'occ', 'parse']
PFBWTF_EXE = './pfbwt-f64' if os.path.exists('./pfbwt-f64') else 'pfbwt-f64'
MERGE_PFP_EXE = './merge_pfp' if os.path.exists('./merge_pfp') else 'merge_pfp'
VCF_SCAN_EXE = './vcf_scan' if os.path.exists('./vcf_scan') else 'vcf_scan'
MERGE_MPS_EXE = './merge_mps' if os.path.exists('./merge_mps') else 'merge_mps'
MPS_TO_MA_EXE = './mps_to_ma' if os.path.exists('./mps_to_ma') else 'mps_to_ma'

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
        self.m = other.ma
        self.mmap = other.mmap
        self.wsize = other.wsize
        self.ma_wsize = other.ma_wsize


class VcfScanCmd:
    def __init__(self, args):
        self.H = args.h
        self.f = args.fasta
        self.S = args.sample
        self.o = args.prefix
        # self.c = args.contigs # TODO: make this an option
        self.stdout = True
        self.ref = False
        self.m = False
        self.vcfs = args.vcf
        self.wsize = args.wsize
        self.ma_wsize = args.ma_wsize

    def set_ref(self):
        self.ref = True

    def set_m(self):
        self.m = True

    def get_cmd(self):
        cmd = [VCF_SCAN_EXE, '-f', self.f, "-o", self.o]
        if self.stdout:
            cmd += ['--stdout']
        if self.ref:
            cmd +=  ['-r']
        else:
            cmd += ["-H", self.H, "-S", self.S]
        if self.m:
            cmd += ['-m']
        cmd += ['-w', str(self.wsize)]
        cmd += ['-x', str(self.ma_wsize)]
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
    cmd = cmd_builder.get_cmd()
    log.write("{}\n".format(' '.join(cmd)))
    fa_fname = full_prefix + ".fa"
    with open(fa_fname, "w") as fa_fp:
        vcf_scan_proc = sp.run(cmd, stdout=fa_fp, stderr=log, check=True)
        vcf_scan_proc.check_returncode()
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
        log.write("{}\n".format((' '.join(vcf_scan_cmd))))
        with open(fa_fname, "w") as fa_fp:
            vcf_scan_proc = sp.run(vcf_scan_cmd, stdout=fa_fp, stderr=log, check=True)
        pfbwt_cmd = [PFBWTF_EXE, '--non-acgt-to-a', '--parse-only', '--print-docs', '-s', '-o', full_prefix, fa_fname]
        if args.mmap:
            pfbwt_cmd = pfbwt_cmd[:1] + ['-m'] + pfbwt_cmd[1:]
        log.write("{}\n".format((' '.join(pfbwt_cmd))))
        pfbwt_proc = sp.run(pfbwt_cmd, check=True, stdout=log, stderr=sp.PIPE)
    else:
        pfbwt_cmd = [PFBWTF_EXE, '--non-acgt-to-a', '--parse-only', '--print-docs', '-s', '-o', full_prefix]
        if args.mmap:
            pfbwt_cmd = pfbwt_cmd[:1] + ['-m'] + pfbwt_cmd[1:]
        log.write("{}\n".format((' '.join(vcf_scan_cmd) + "|" + ' '.join(pfbwt_cmd))))
        vcf_scan_proc = sp.Popen(vcf_scan_cmd, stdout=sp.PIPE, stderr=log)
        pfbwt_proc = sp.run(pfbwt_cmd, stdin=vcf_scan_proc.stdout, stdout=log, stderr=sp.PIPE, check=True)
        vcf_scan_proc.wait()
        if vcf_scan_proc.returncode != 0:
            sys.stderr.write("vcf scan failed w/ error: {}\n".format(vcf_scan_proc.returncode))
            raise Exception
        if pfbwt_proc.returncode != 0:
            sys.stderr.write("pfbwt parsing failed w/ error: {}\n".format(pfbwt_proc.returncode))
            raise Exception
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


def merge_mps(args, thread_args, logger, log_fp):
    logger.info("merging marker positions")
    cmd = [MERGE_MPS_EXE, args.o + ".mps", args.o + ".ref"] + [a.full_prefix for a in thread_args]
    logger.info(" ".join(cmd))
    sp.run(cmd, check=True, stdout=log_fp, stderr=sp.PIPE)
    logger.info("merged marker positions")


class PfbwtCmd:

    def __init__(self, args):
        self.wsize = str(args.wsize)
        self.mod = str(args.mod)
        self.sa = args.sa
        self.rssa = args.rssa
        self.ma = args.ma
        self.mmap = args.mmap
        self.o = args.o

    def get_cmd(self):
        cmd = [PFBWTF_EXE, '--pfbwt-only', '--print-docs', '-o', self.o, '-w', self.wsize, '-m', self.mod]
        if self.sa or self.ma:
            cmd += ['--stdout', 'sa', '-s']
        if self.mmap:
            cmd += ['-m']
        if self.rssa:
            cmd += ['-r']
        return cmd


def vcf_to_bwt(args):
    sp.run(["samtools", "faidx", args.fasta], check=True)
    if not args.samples:
        sys.stderr.write("no sample file specified, defaulting to samples first VCF {}\n".format(args.vcf[0]))
        samples = list(pysam.VariantFile(args.vcf[0], "r").header.samples)
    else:
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
        if args.ma:
            merge_mps(args, thread_args, logger, log_fp)
        logger.info("generating parse (directly from fasta files)")
        parse_cmd = [PFBWTF_EXE, '--parse-only', '-s', '--print-docs', '-o', args.o]
        cat_cmd = ['cat'] + [args.o + ".ref.fa"] + [a.full_prefix + ".fa" for a in thread_args]
        cat_proc = sp.Popen(cat_cmd, stdout=sp.PIPE, stderr=log_fp)
        logger.info(" ".join(cat_cmd) + " | " + " ".join(parse_cmd))
        parse_proc = sp.run(parse_cmd, check=True, stdin=cat_proc.stdout, stderr=sp.STDOUT, stdout=log_fp)
        cat_proc.wait()
        if cat_proc.returncode != 0:
            raise Exception
        if parse_proc.returncode != 0:
            sys.stderr.write("parsing failed\n")
            raise Exception
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
        if args.ma:
            merge_mps(args, thread_args, logger, log_fp)
        logger.info("merging parses")
        merge_pfp_cmd = [MERGE_PFP_EXE, '-w', str(args.wsize), '-s', '--parse-bwt', '--docs', '-o', args.o, '-t', str(args.threads)] + all_prefixes
        logger.info(" ".join(merge_pfp_cmd))
        try:
            sp.run(merge_pfp_cmd, stdout=log_fp, stderr=sp.PIPE, check=True)
        except sp.CalledProcessError as e:
            logger.info("MERGE_PFP ERROR: {} {}".format(e.returncode, e.stderr.decode()))
            exit(1)
        logger.info("done merging parses")
        if args.clean:
            logger.info("cleaning files")
            for prefix in all_prefixes:
                clean_parse_files(prefix)
    # construct BWT
    logger.info("constructing BWT")
    pfbwt_cmd_builder = PfbwtCmd(args)
    pfbwt_cmd = pfbwt_cmd_builder.get_cmd()
    pfbwt_proc = sp.Popen(pfbwt_cmd, stdout=sp.PIPE, stderr=log_fp)
    if args.sa and not args.ma:
        # cat output to <output>.sa, don't call mps_to_ma
        sa_fp = open(args.o + ".sa", "wb")
        logger.info(" ".join(pfbwt_cmd) + " | cat")
        cat_sa_proc = sp.run(['cat'], stdin=pfbwt_proc.stdout, stdout=sa_fp, stderr=log_fp, check=True)
        cat_sa_proc.check_returncode()
        pfbwt_proc.wait()
        sa_fp.close()
    elif args.ma and not args.sa:
        logger.info("constructing marker array along with BWT. SA will not be saved")
        marker_array_cmd = [MPS_TO_MA_EXE, "-o", args.o + ".ma", args.o + ".mps", '-']
        if args.mmap:
            marker_array_cmd = marker_array_cmd[:1] + ['-m'] + marker_array_cmd[1:]
        logger.info(" ".join(pfbwt_cmd) + " | " + " ".join(marker_array_cmd))
        marker_array_proc = sp.run(marker_array_cmd, stdin=pfbwt_proc.stdout, stderr=sp.STDOUT, stdout=log_fp, check=True)
        marker_array_proc.check_returncode()
        pfbwt_proc.wait()
    elif args.sa and args.ma:
        logger.info("constructing marker array along with BWT. SA will also be saved")
        tee_proc = sp.Popen(['tee', args.o + ".sa"], stdin=pfbwt_proc.stdout, stderr=log_fp, stdout=sp.PIPE)
        marker_array_cmd = [MPS_TO_MA_EXE, "-o", args.o + ".ma", args.o + ".mps", '-']
        if args.mmap:
            marker_array_cmd = marker_array_cmd[:1] + ['-m'] + marker_array_cmd[1:]
        logger.info(" ".join(pfbwt_cmd) + " | tee " + args.o + ".sa" + " | " + " ".join(marker_array_cmd))
        marker_array_proc = sp.run(marker_array_cmd, stdin=tee_proc.stdout, stderr=sp.STDOUT, stdout=log_fp, check=True)
        tee_proc.wait()
        pfbwt_proc.wait()
    else:
        logger.info(" ".join(pfbwt_cmd))
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
    parser.add_argument("--samples", "-S",
                        help="file containing new-line delimited desired samples in vcf. defaults to all")
    parser.add_argument("--threads", "-t", type=int, default=int(1), help="number of threads (default: 1)")
    parser.add_argument("--save_fasta", "-f", action='store_true', help="store fasta sequences generated from VCFs")
    parser.add_argument("-o", default="out", help="output prefix")
    parser.add_argument("--no_merge", action='store_true',
                        help="generate a BWT from a non-merged parse of text collection")
    parser.add_argument("--clean", action='store_true', help="cleanup intermediate files as we go")
    parser.add_argument("--ma", "-m", action='store_true', help="build marker array")
    parser.add_argument("--keep_parse", action='store_true', help="keeps the final parse (for when --clean is used)")
    parser.add_argument("-s", "--sa", action='store_true', help="save SA to <output>.sa")
    parser.add_argument("-r", "--rssa", action='store_true',
                        help="save run-length sampled SA to <output>.ssa (run-starts) and <output>.esa (run-ends)")
    parser.add_argument("--mmap", '-M', action='store_true',
                        help="tell pfbwt-f64 to use mmap (use this for very large files)")
    parser.add_argument("--ma_wsize", default=10, type=int, help="window size to use for marker array")
    parser.add_argument("--wsize", default=10, type=int, help="window size for prefix-free parsing")
    parser.add_argument("--mod", default=10, type=int, help="mod for prefix-free parsing")
    args = parser.parse_args()

    vcf_to_bwt(args)
