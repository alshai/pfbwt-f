import sys
import argparse
import pysam
import subprocess
from pysuffixarray.core import SuffixArray

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta")
    parser.add_argument("vcf")
    parser.add_argument("-w", default=10, type=int, help="number of As to interleave between sequences")
    parser.add_argument("-o", default="out", help="output prefix")
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf, "r")
    contigs = [f.split()[0].strip() for f in open(args.fasta + ".fai", "r").readlines()]
    out_seq = open(args.o + ".txt", "w")
    out_sa = open(args.o + ".sa", "w")
    out_bwt = open(args.o + ".bwt", "wb")
    out_markers = open(args.o + ".markers", "w")
    to_append = "A"*args.w

    big_seq = ""
    headers = []
    markers = {}
    # add reference sequence first
    for contig in contigs:
        pos = len(big_seq)
        for rec in vcf.fetch(contig):
            gt = 0
            markers[pos + rec.start]  = (rec.rid, rec.start, gt, rec.alleles[gt])
            pos += len(rec.alleles[gt]) - 1
        # append contig as is to big_seq
        p1 = subprocess.Popen(["samtools", "faidx", args.fasta, contig], stdout=subprocess.PIPE, text=True)
        seq, p_err = p1.communicate()
        fasta_lines = seq.strip().split('\n')
        headers.append(fasta_lines[0])
        big_seq += "".join(fasta_lines[1:]) + to_append

    for sample in vcf.header.samples:
        for h in [0, 1]:
            for contig in contigs:
                # first gather the markers
                pos = len(big_seq)
                for rec in vcf.fetch(contig):
                    gt = rec.samples[sample]["GT"][h]
                    markers[pos + rec.start]  = (rec.rid, rec.start, gt, rec.alleles[gt])
                    pos += len(rec.alleles[gt]) - 1
                # then generate the string
                p1 = subprocess.Popen(["samtools", "faidx", args.fasta, contig], stdout=subprocess.PIPE)
                p2 = subprocess.Popen(["bcftools", "consensus", "-p", "{}.{}.".format(sample, h), "-f", "-", "-H", str(h+1), "-s", sample, args.vcf], stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
                seq, p_err = p2.communicate()
                fasta_lines = seq.strip().split("\n")
                headers.append(fasta_lines[0])
                big_seq += "".join(fasta_lines[1:]) + to_append

    for p, (c, rp, g, a) in markers.items():
        # print(p, c, rp, g, big_seq[p], a)
        assert(big_seq[p:p+len(a)] == a)
    sa = SuffixArray(big_seq).suffix_array()
    for i, s in enumerate(sa):
        if s in markers:
            c, rp, g, a = markers[s]
            out_markers.write("{} {} {} {}\n".format(i, c, rp, g))
    for s in sa:
        out_bwt.write(big_seq[s-1].encode() if s != 0 else bytes.fromhex('00'))
    out_seq.write("{}\n".format(big_seq))
    out_sa.write("\n".join(str(i) for i in sa))
    vcf.close()

    out_seq.close()
    out_bwt.close()
    out_sa.close()
    out_markers.close()