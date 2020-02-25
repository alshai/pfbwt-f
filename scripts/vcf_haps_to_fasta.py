from cyvcf2 import VCF
import argparse
from Bio import SeqIO
from Bio.Seq import MutableSeq
import sys
def die(msg):
    sys.stderr.write(msg + "\n")
    exit(1)

def get_samples_from_bcf(fname):
    bcf = VCF(fname, mode='r')
    # TODO: check if vcf is sorted and indexed
    samples = list(bcf.samples)
    bcf.close()
    return samples

def write_marker(fp, hpos, rpos, allele, bytes=8, endian="little"):
    fp.write(hpos.to_bytes(bytes, "little"))
    fp.write(rpos.to_bytes(bytes, "little"))
    fp.write(allele.to_bytes(bytes, "little"))

# assumes diploid for now
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta")
    parser.add_argument("vcf")
    parser.add_argument("--stdout", action='store_true', help="output fasta to stdout")
    parser.add_argument("--samples", help="file name specifying which samples to print")
    parser.add_argument("--nsamples", type=int, default=None, help="print first n samples from this VCF")
    parser.add_argument("-o", help="output prefix", default="haplotypes")
    parser.add_argument("--endian", default="little")
    parser.add_argument("--bytes", type=int, default=8)
    # add optional arguments here
    args = parser.parse_args()
    if args.samples:
        samples = open(args.samples).read().strip().split("\n")
    elif args.nsamples:
        samples = get_samples_from_bcf(args.vcf)[:args.nsamples]
    else:
        samples = get_samples_from_bcf(args.vcf)
    if args.fasta == "-":
        fasta_in = sys.stdin
    else:
        fasta_in = open(args.fasta)
    pos = 0
    if args.stdout:
        fasta_fp = sys.stdout
    else:
        fasta_fp = open(args.o + ".fa", "w")
    marker_fp = open(args.o + ".markers", "wb")
    # TODO: check if vcf is sorted and indexed
    for sample in samples:
        for i in range(2): # process two haplotypes
            for fasta_record in SeqIO.parse(fasta_in, "fasta"):
                prec_end = 0
                prec_start = 0
                seq = fasta_record.seq
                contig = fasta_record.name
                fasta_fp.write(">{}.{}.{}\n".format(sample, i+1, contig))
                bcf = VCF(args.vcf)
                bcf.set_samples([sample])
                recs = bcf.__call__(contig)
                for rec in recs:
                    gts = rec.genotypes[0]
                    if (len(gts) != 3):
                        sys.stderr.write("error: number of genotypes for {} at marker {} is not 2!\n".format(sample, rec.ID))
                        exit(1)
                    if gts[i]:
                        if rec.start in range(prec_start, prec_end):
                            sys.stderr.write("warning: skipping {} at {}:{} because it overlaps a previous marker\n".format(rec.ID, rec.CHROM, rec.start))
                            continue
                        pos += rec.start - prec_end # ffwd current position
                        if len(rec.REF) >= len(rec.ALT[gts[i]-1]): # del and snp
                            write_marker(marker_fp, pos, rec.start, gts[i], args.bytes, args.endian)
                            pos += 1
                        elif len(rec.REF) < len(rec.ALT[gts[i]-1]):
                            # insertion
                            ins_size = len(rec.ALT[gts[i]-1]) - len(rec.REF) + 1
                            for j in range(ins_size):
                                write_marker(marker_fp, pos, rec.start, gts[i], args.bytes, args.endian)
                                pos += 1
                        fasta_fp.write(str(seq[prec_end:rec.start]) + str(rec.ALT[gts[i]-1]))
                        prec_start = rec.start
                        # prec_end skips over to position just after record
                        prec_end = rec.start + len(rec.REF)
                bcf.close()
                fasta_fp.write(str(seq[prec_end:]) + "\n")
                pos += len(seq) - prec_end
            fasta_in.seek(0)
    marker_fp.close()
    fasta_fp.close()
    fasta_in.close()
