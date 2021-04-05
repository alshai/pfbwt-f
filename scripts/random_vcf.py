import sys
import subprocess
import argparse
import os
import random

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def vcf_write_header(vcf_fp, fasta_fn, nsamples):
    contigs = []
    if not os.path.exists(fasta_fn + ".fai"):
        sys.stderr.write("no fai file found. scanning file to get seq names and lengths...\n")
        for name, seq, qual in readfq(open(fasta_fn, "r")):
            contigs.append((name, len(seq)))
    else:
        sys.stderr.write("reading {}.fai\n".format(fasta_fn))
        for line in open(fasta_fn + ".fai", "r"):
            fields = line.strip().split()
            contigs.append((fields[0], int(fields[1])))
    vcf_fp.write('##fileformat=VCFv4.3\n')
    vcf_fp.write('##source=random_vcf.py\n')
    vcf_fp.write('##reference={}\n'.format(fasta_fn))
    for name, length in contigs:
        vcf_fp.write('##contig=<ID={},length={}>\n'.format(name, length))
    vcf_fp.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf_fp.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')
    for i in range(nsamples):
        vcf_fp.write('\tsample{}'.format(i))
    vcf_fp.write('\n')


def vcf_from_fasta(args):
    if args.out:
        out = args.out
    else:
        out = args.fasta
    vcf = open(out + ".vcf", "w")
    vcf_write_header(vcf, args.fasta, args.nsamples)
    ACGT = set(["A", "C", "G", "T"])
    for name, seq, qual in readfq(open(args.fasta, "r")):
        i = 0
        nvars = int(args.var_density * len(seq))
        sys.stderr.write("will generate {} variants for seq of length {}\n".format(nvars, len(seq)))
        pos_v = sorted(random.sample(range(len(seq)), nvars))
        for pos in pos_v:
            # TODO: support indels
            # TODO: support >1 alt alleles
            ref = seq[pos].upper()
            alt = random.choice(list(ACGT - set(ref)))
            vcf.write('{}\t{}\tvar{}\t{}\t{}\t.\t.\t.\tGT'.format(name, pos+1, i, ref, alt))
            # generate genotypes here
            for j in range(args.nsamples * args.ploidy):
                g = 1 if random.random() < args.allele_freq else 0
                if j % 2:
                    vcf.write('{}'.format(g))
                else:
                    vcf.write('\t{}|'.format(g))
            vcf.write('\n')
            i += 1
    vcf.close()
    subprocess.run(["bcftools", "view",  "-O", "z", out + ".vcf"], stdout=open(out+".vcf.gz", "w"))
    subprocess.run(["bcftools", "index", out + ".vcf.gz"])

"""
usage: python random_vcf.py <fasta> <nsamples>
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta")
    parser.add_argument("--nsamples", "-n", type=int, default=1)
    parser.add_argument("--ploidy", "-p", type=int, default=2)
    parser.add_argument("--out", "-o")
    parser.add_argument("--var_density", type=float, default=0.001, help="fraction of desired variants w.r.t reference length")
    parser.add_argument("--allele_freq", type=float, default=0.5, help="desired alt allele frequency")
    args = parser.parse_args()

    if args.ploidy % 2:
        sys.stderr.write("currently odd ploidy is not supported!\n")
        exit(1)

    vcf_from_fasta(args)
