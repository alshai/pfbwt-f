from cyvcf2 import VCF
import argparse
from Bio import SeqIO
from Bio.Seq import MutableSeq
import sys
import subprocess
from scipy import sparse
import numpy as np

""" pfbwt_marker_array.py

Author: Taher Mun
Date: Feb 27th, 2020

usage: python pfbwt_marker_array.py <fasta> <VCF>

Given a VCF file describing a set of haplotypes with respect to a reference
genome, generate the BWT (and optinally the SA) of the concatenation of these
haplotype sequences using Prefix-Free Parsing algorithm, and also generate
and store the marker array.

The Marker Array (MA) is an auxillary structure to the Suffix Array, much like
the Document Array.  It stores, for each variant contained in each of the
haplotypes, its corresponding position in the reference sequence. 
For example, if SA[i] lies at a SNP in some haplotype, MA[i] contains the
position in the reference sequence of that SNP. 

Technically, the marker array is the same length as the SA, but since these
data structures can end up being very large, we only store marker values at
known variant positions, so it is a sparse array.

As mentioned above, if the variant at SA[i] is a SNP, then MA[i] contains the
position of that SNP in the reference sequence.

If SA[i] lies just before a deletion, then MA[i] contains the position in the
reference sequence that lies just before the deletion.

If SA[i] lies just before or inside an insertion, MA[i] contains the position
in the reference sequence just before the insertion.

The BWT tends to contain runs of identical characters which correspond to
suffixes with the similar right-contexts. Due to this unique property, we would
expect a run in the BWT to represent, for the most part, corresponding posiions
within group of similar haplotypes. So if a group of haplotypes contain the
same variant at a certain position, the variant positions for each of these
haplotypes would probably be grouped in the same run of the BWT. By this logic,
the marker array would probably be highly repetitive at regions corresponding
to BWT runs, making it highly compressible.
"""

"""
credit: https://github.com/lh3/readfq/
"""
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


def die(msg):
    sys.stderr.write(msg + "\n")
    exit(1)


def get_samples_from_bcf(fname):
    bcf = VCF(fname, mode='r')
    samples = list(bcf.samples)
    bcf.close()
    return samples


def write_marker(fp, suffix_pos, ref_pos, allele, bytes=8, endian="little"):
    fp.write(suffix_pos.to_bytes(bytes, "little"))
    fp.write(ref_pos.to_bytes(bytes, "little"))


class VCFToFastaMarkersArgs:
    """
    input: return value of ArgumentParser.parse_args
    """
    def __init__(self, args):
        self.fasta = args.fasta
        self.vcf = args.vcf
        self.samples = args.samples
        self.nsamples = args.nsamples
        self.stdout = args.stdout
        self.bytes = args.bytes
        self.endian = args.endian
        if args.out:
            self.out = args.out
        else:
            self.out = args.fasta
        self.save_fasta = args.save_fasta


"""
Generator for converting VCF to sequence
input: VCFToFastaMarkersArgs object
results: yields, for each marker in VCF, the DNA segment that starts from just
         after the end of the last marker up to the end of the current marker
         As a side effect:
         optionally writes marker positions w/i the total concatenated strings
         along with the positions wrt the reference to args.out + ".markers",
         optionally writes fasta file representing whole haplotype sequences in
         VCF to args.out + ".from_vcf.fa"
"""
def vcf_to_fasta_markers(args):
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
    if args.save_fasta:
        fasta_out = open(args.out + ".from_vcf.fa", "w")
    marker_out = open(args.out + ".markers", "wb")
    # TODO: check if vcf is sorted and indexed
    for sample in samples:
        for i in range(2): # process two haplotypes
            for name, seq, qual in readfq(fasta_in):
            # for fasta_record in SeqIO.parse(fasta_in, "fasta"):
                prec_end = 0
                prec_start = 0
                # seq = fasta_record.seq
                # contig = fasta_record.name
                # to_write = ">{}.{}.{}\n".format(sample, i+1, contig)
                to_write = ">{}.{}.{}\n".format(sample, i+1, name)
                if args.save_fasta:
                    fasta_out.write(to_write)
                yield to_write.encode()
                bcf = VCF(args.vcf)
                bcf.set_samples([sample])
                # recs = bcf.__call__(contig)
                recs = bcf.__call__(name)
                for rec in recs:
                    gts = rec.genotypes[0]
                    if (len(gts) != 3):
                        die("error: number of genotypes for {} at marker {} is not 2!".format(sample, rec.ID))
                    if gts[i]:
                        if rec.start in range(prec_start, prec_end):
                            sys.stderr.write("warning: skipping {} at {}:{} because it overlaps a previous marker\n".format(rec.ID, rec.CHROM, rec.start))
                            continue
                        pos += rec.start - prec_end # ffwd current position
                        if len(rec.REF) >= len(rec.ALT[gts[i]-1]): # del and snp
                            write_marker(marker_out, pos, rec.start, gts[i], args.bytes, args.endian)
                            pos += 1
                        elif len(rec.REF) < len(rec.ALT[gts[i]-1]):
                            # insertion
                            ins_size = len(rec.ALT[gts[i]-1]) - len(rec.REF) + 1
                            for j in range(ins_size):
                                write_marker(marker_out, pos, rec.start, gts[i], args.bytes, args.endian)
                                pos += 1
                        # to_write = str(seq[prec_end:rec.start]) + str(rec.ALT[gts[i]-1])
                        to_write = seq[prec_end:rec.start] + rec.ALT[gts[i]-1]
                        if args.save_fasta:
                            fasta_out.write(to_write)
                        yield to_write.encode()
                        prec_start = rec.start
                        # prec_end skips over to position just after record
                        prec_end = rec.start + len(rec.REF)
                bcf.close()
                pos += len(seq) - prec_end
                to_write = seq[prec_end:] + "\n"
                if args.save_fasta:
                    fasta_out.write(to_write)
                yield to_write.encode()
            fasta_in.seek(0)
    if args.save_fasta:
        fasta_out.close()
    marker_out.close()
    fasta_in.close()


class MarkerArrayArgs:
    def __init__(self, args, sa_fp, err_fp):
        self.bytes = args.bytes
        self.endian = args.endian
        if args.out:
            self.out = args.out
        else:
            self.out = args.fasta
        self.sa_fp = sa_fp
        self.err_fp = err_fp
        self.save_sa = args.save_sa


# TODO add support for 2+ alt alleles
def markers_to_sparse(f, nbytes=8, endian="little"):
    if nbytes == 8:
        dtype='uint64'
    elif nbytes == 4:
        dtype='uint32'
    f.seek(0,2)
    fsize = f.tell()
    if (fsize % nbytes) or (fsize % 2):
        sys.stderr.write("invalid file ({} bytes)!\n".format(fsize))
        exit(1)
    f.seek(0,0)
    marr = np.ndarray((int(fsize / nbytes / 2), 2), dtype=dtype)
    bytestr = f.read(nbytes)
    i=0
    while bytestr:
        x = int.from_bytes(bytestr, endian)
        marr[int(i/2), i % 2] = x
        bytestr = f.read(nbytes)
        i += 1
    return sparse.dok_matrix(sparse.coo_matrix((marr[:,1], (marr[:,0], np.zeros((marr.shape[0]), dtype=dtype)))))


"""
rearrange markers by permutation set by SA
input: MarkerArrayArgs
results: writes marker array to args.out + ".ma"
         optionally writes suffix array to args.out + ".sa"
         ".ma" file formatted as follows:
         <Array Position> <Marker Value>
"""
def marker_array(args):
    i = int(0)
    nbytes = args.bytes
    endian = args.endian
    sa_in = args.sa_fp
    err_out = args.err_fp
    # read markers and save to dok_matrix
    with open(args.out + ".markers", "rb") as fp:
        marr = markers_to_sparse(fp, nbytes, endian)
    if args.save_sa:
        sa_out = open(args.out + ".sa", "wb")
    ma_out = open(args.out + ".ma", "wb")
    bytestr = sa_in.read(nbytes)
    while bytestr:
        s = int.from_bytes(bytestr, endian)
        if args.save_sa:
            sa_out.write(bytestr)
        if s < marr.shape[0]:
            marker = marr[s,0]
            if marker:
                ma_out.write(i.to_bytes(nbytes, endian))
                ma_out.write(marker.tobytes())
        bytestr = sa_in.read(nbytes)
        i += 1
    ma_out.close()
    if args.save_sa:
        sa_out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="fasta file containing reference sequence")
    parser.add_argument("vcf", help="VCF file containing haplotypes. Must be sorted and indexed")
    parser.add_argument("--stdout", action='store_true', help="output fasta to stdout")
    parser.add_argument("--samples", help="file name specifying which samples to print")
    parser.add_argument("--nsamples", type=int, default=None, help="print first n samples from this VCF")
    parser.add_argument("--out", "-o", help="output prefix (<fasta> by default)")
    parser.add_argument("--endian", default="little", choices=["big", "little"])
    parser.add_argument("--bytes", type=int, default=8, choices=[4,8])
    parser.add_argument("--save_fasta", action='store_true', default=False)
    parser.add_argument("--save_sa", action='store_true', default=False)
    args = parser.parse_args()

    # sanity check arguments
    if args.bytes == 8:
        exe = "./pfbwt-f64"
    elif args.bytes == 4:
        exe == "./pfbwt-f"
    else:
        die("--bytes can only be 4 or 8!") # redundant
    if args.out:
        out = args.out
    else:
        out = args.fasta

    # feed haplotypes into parsinng step
    parse_cmd = [exe, "--parse-only", "-o", args.fasta, "-s"]
    err1 = open("{}.parse.err".format(out), "w")
    with subprocess.Popen(parse_cmd, stdin=subprocess.PIPE, stderr=err1) as p:
        for s in vcf_to_fasta_markers(VCFToFastaMarkersArgs(args)):
            p.stdin.write(s)
        p.stdin.close()
        if p.wait():
            err1.write("error in pfbwtf --parse-only!\n")
        err1.close()

    # build the BWT and generate marker array
    pfbwt_cmd = [exe, "--pfbwt-only", "-o", args.fasta, "-s", "--stdout", "sa"]
    err2 = open("{}.pfbwt.err".format(out), "w")
    with subprocess.Popen(pfbwt_cmd, stderr=err2, stdout=subprocess.PIPE) as p:
        marker_array(MarkerArrayArgs(args, p.stdout, err2))
        err2.close()
