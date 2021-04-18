import sys
import argparse
import pysam
import subprocess
import os
from pysuffixarray.core import SuffixArray

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta")
    parser.add_argument("vcf")
    parser.add_argument("-w", default=10, type=int, help="number of As to interleave between sequences")
    parser.add_argument("-o", default="out", help="output prefix")
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf, "r")
    if not os.path.exists(args.fasta + ".fai"):
        p = subprocess.run(["samtools", "faidx", args.fasta], check=True)
        p.check_returncode()
    contigs = [f.split()[0].strip() for f in open(args.fasta + ".fai", "r").readlines()]
    out_seq = open(args.o + ".txt", "w")
    out_sa = open(args.o + ".sa", "w")
    out_bwt = open(args.o + ".bwt", "wb")
    out_markers = open(args.o + ".markers", "w")
    to_append = "A"*args.w

    big_seq = ""
    headers = []
    markers = {}
    pos = 0
    # add reference sequence first
    for contig in contigs:
        for rec in vcf.fetch(contig):
            gt = 0
            ale_string = "{}->{}".format(rec.alleles[0], rec.alleles[1])
            if len(rec.alleles[0]) == 1 and len(rec.alleles[1]) == 1:
                markers[pos + rec.start]  = (rec.rid, rec.start, gt, ale_string)
            elif len(rec.alleles[0]) != len(rec.alleles[1]): # indels
                for i in range(rec.start, rec.stop + 1): # record marker at breakpoints and in between if del
                    markers[pos + i]  = (rec.rid, rec.start, gt, ale_string)
            else:
                sys.stderr.write("skipping marker at {}:{} - not SNP or indel\n".format(rec.contig, rec.start+1))
        pos += vcf.header.contigs[contig].length + len(to_append)
        # append contig as is to big_seq
        p1 = subprocess.Popen(["samtools", "faidx", args.fasta, contig], stdout=subprocess.PIPE, text=True)
        seq, p_err = p1.communicate()
        fasta_lines = seq.strip().split('\n')
        headers.append(fasta_lines[0])
        big_seq += "".join(fasta_lines[1:]) + to_append

    sys.stdout.write("{}\n".format(pos))
    for sample in vcf.header.samples:
        for h in [0, 1]:
            for contig in contigs:
                # first gather the markers
                pos = len(big_seq)
                bias = 0
                for rec in vcf.fetch(contig):
                    gt = rec.samples[sample]["GT"][h]
                    # remember that rlen and alen includes base before indel
                    rlen = len(rec.alleles[0])
                    alen = len(rec.alleles[1])
                    ale_string = "{}->{}".format(rec.alleles[0], rec.alleles[1])
                    if rlen == 1 and alen == 1:
                        markers[pos + bias + rec.start]  = (rec.rid, rec.start, gt, ale_string)
                    elif rlen != alen and gt == 0: # indel
                        for i in range(rec.start, rec.stop + 1):
                            # if del, includes everything between breakpoints 
                            markers[pos + bias + i]  = (rec.rid, rec.start, gt, ale_string)
                    elif rlen > alen and gt > 0: # deletion
                        markers[pos + bias + rec.start] = (rec.rid, rec.start, gt, ale_string)
                        markers[pos + bias + rec.start + 1] = (rec.rid, rec.start, gt, ale_string)
                        bias = bias - (rlen - 1)
                    elif rlen < alen and gt > 0: # insertion
                        for i in range(0, alen + 1):
                            markers[pos + bias + rec.start + i] = (rec.rid, rec.start, gt, ale_string)
                        bias += alen - 1
                    else:
                        sys.stderr.write("skipping marker at {}:{}\n".format(rec.contig, rec.start+1))
                pos += vcf.header.contigs[contig].length + len(to_append) + bias
                # then generate the string
                p1 = subprocess.Popen(["samtools", "faidx", args.fasta, contig], stdout=subprocess.PIPE)
                p2 = subprocess.Popen(["bcftools", "consensus", "-p", "{}.{}.".format(sample, h), "-f", "-", "-H", str(h+1), "-s", sample, args.vcf], stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
                seq, p_err = p2.communicate()
                fasta_lines = seq.strip().split("\n")
                headers.append(fasta_lines[0])
                big_seq += "".join(fasta_lines[1:]) + to_append

    for pos, (rid, start, gt, allele) in markers.items():
        sys.stdout.write("{}: {} {} {} {}\n".format(pos, rid, start, gt, allele))
    sys.stderr.write("building suffix array\n")
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
