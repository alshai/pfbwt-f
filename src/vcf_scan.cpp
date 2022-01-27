#include <deque>
#include <vector>
#include <string>
#include <getopt.h>
#include "vcf_scanner.hpp"
#include "marker_array.hpp"

struct Args {
    int w = 10;
    int ma_w = 1;
    std::string fa_fname;
    std::vector<std::string> vcf_fnames;
    std::vector<std::string> contigs;
    std::vector<std::string> samples;
    std::string ref_fasta;
    std::string out = "out";
    int nthreads = 1;
    int verb = 0;
    int haplotype = 0;
    int ref_only = 0;
    int to_stdout = 0;
    int mai = 0;
};

VCFScannerArgs ArgsToVCFScannerArgs(Args args) {
    VCFScannerArgs v;
    v.vcf_fnames = args.vcf_fnames;
    v.contigs = args.contigs;
    v.verb = args.verb;
    v.ref_fasta = args.ref_fasta;
    return v;
}

Args parse_args(int argc, char** argv) {
    Args args;
    std::vector<std::string> vcf_fnames;
    std::vector<std::string> contigs;
    std::vector<std::string> samples;
    int c;
    char* pch;
    static struct option lopts[] = {
        {"window-size", required_argument, NULL, 'w'},
        {"marker-window-size", required_argument, NULL, 'x'},
        {"output", required_argument, NULL, 'o'},
        {"threads", required_argument, NULL, 't'},
        {"contigs", required_argument, NULL, 'c'},
        {"samples", required_argument, NULL, 'S'},
        // {"sample", required_argument, NULL, 's'},
        {"verbose", no_argument, NULL, 'v'},
        {"fasta", required_argument, NULL, 'f'},
        {"ref-only", no_argument, NULL, 'r'},
        {"stdout", no_argument, &args.to_stdout, 1},
        {"marker-index", no_argument, NULL, 'm'}
    };

    while ((c = getopt_long( argc, argv, "rf:w:x:o:t:c:S:vH:m", lopts, NULL) ) != -1) {
        switch(c) {
            case 'f':
                args.ref_fasta = optarg; break;
            // case 's':
            //     args.samples = {optarg}; break;
            case 'w':
                args.w = atoi(optarg); break;
            case 'x':
                args.ma_w = atoi(optarg); break;
            case 'o':
                args.out = optarg; break;
            case 't':
                args.nthreads = atoi(optarg); break;
            case 'c':
                pch = strtok(optarg, ",");
                while (pch != NULL) {
                    contigs.push_back(pch);
                    pch = strtok(NULL, ",");
                }
                break;
            case 'S':
                pch = strtok(optarg, ",");
                while (pch != NULL) {
                    samples.push_back(pch);
                    pch = strtok(NULL, ",");
                }
                break;
            case 'v':
                args.verb = 1; break;
            case 'H':
                args.haplotype = atoi(optarg);
                break;
            case 'r':
                args.ref_only = 1; break;
            case 'm':
                args.mai = 1; break;
            case '?':
                fprintf(stderr,  "Unknown option.\n");
                exit(1);
                break;
            case ':':
                fprintf(stderr,  "no argument specified for option\n");
                exit(1);
        }
    }
    for (int i = optind; i < argc; ++i) {
        vcf_fnames.push_back(argv[optind++]);
    }
    args.vcf_fnames = vcf_fnames;
    args.contigs = contigs;
    args.samples = samples;
    if (args.out == "") {
        fprintf(stderr, "no output prefix specified. Defaulting to 'out'\n");
    }
    return args;
}

void update_sequence(char* seq, int rlen, bcf1_t* rec, int ppos, int32_t gt, FILE* ofp, FILE* log) {
    (void) log;
    if (rec != NULL && gt > -1)  { // normal case
        fwrite(seq+ppos, sizeof(char), rec->pos - ppos, ofp); // write ref sequence
        char* allele = rec->d.allele[gt];
        fwrite(allele, sizeof(char), strlen(allele), ofp); // add allele
    } else { // end case
        fwrite(seq+ppos, sizeof(char), rlen - ppos, ofp); // write ref sequence
        fputc('\n', ofp);
    }
}

// produces fasta and marker array for one haplotypes of a single simple
void scan_vcf_sample(Args args, std::string sample) {
    VCFScannerArgs vargs(ArgsToVCFScannerArgs(args));
    vargs.sample = sample;
    FILE *log = NULL, *ma_fp = NULL, *fa_fp = NULL;
    int i = args.haplotype;
    std::string fa_fname = args.out + "." + sample + "." + std::to_string(i) + ".fa";
    std::string fa_header = args.out + "." + sample + "." + std::to_string(i);
    std::string ma_fname = args.out + "." + sample + "." + std::to_string(i) + ".mps";
    std::string log_fname = args.out + "." + sample + "." + std::to_string(i) + ".log";
    if (args.ref_only) {
        ma_fname = args.out + ".ref.mps";
        log_fname = args.out + ".ref.log";
        fa_fname = args.out + ".ref.fa";
    }
    log = fopen(log_fname.data(), "w");
    if (args.mai) ma_fp = fopen(ma_fname.data(), "wb");
    if (args.to_stdout) {
        fa_fp = stdout;
    } else {
        fa_fp = fopen(fa_fname.data(), "w");
    }

    MarkerPositionsWriter mi_writer(([&]() {
            if (args.mai) return MarkerPositionsWriter(args.ma_w, ma_fp, args.verb ? NULL : log);
            else return MarkerPositionsWriter();
        })()
    );
    int ppos = 0;
    int ppos_after = 0;
    int64_t seq_start = 0;
    int64_t bias = 0;
    std::string pseq("");
    int prange[2] = {0,0};
    int range[2] = {0,0};

    auto out_fn = [&](bcf_hdr_t* hdr, bcf1_t* rec, BCFGenotype& gtv, std::vector<size_t>& posv, char* ref_seq, int32_t ref_len, int rid) {
        (void) posv;
        const char* ref_name = bcf_hdr_id2name(hdr, rid);
        if (strcmp(ref_name, pseq.data())) {
            std::string to_write(ref_name);
            if (!args.ref_only) {
                to_write = sample + "." + std::to_string(args.haplotype) + "." + to_write;
            }
            fprintf(fa_fp, ">%s\n", to_write.data());
            pseq.assign(ref_name);
        }

        if (rec != NULL) {
            range[0] = rec->pos;
            range[1] = rec->pos + strlen(rec->d.allele[0]);
            if (prange[1] <= range[0]) {
                int gt =  args.ref_only ? 0 : gtv[i];
                gt = (gt == -1) ? 0 : gt; // hack for when gt is malformed. ie: 0 1|1 0|0 etc
                if (args.mai) {
                    // if indel, [ra]len includes base before indel
                    size_t rlen = strlen(rec->d.allele[0]);
                    size_t alen = strlen(rec->d.allele[1]);
                    if (rlen == 1 && alen == 1) {
                        mi_writer.update(static_cast<size_t>(seq_start + bias + rec->pos), rec->pos, gt, rid);
                    } else if (rlen != alen and gt == 0) {
                        for (size_t i = 0; i <= rlen; ++i) {
                            mi_writer.update(static_cast<size_t>(seq_start + bias + rec->pos + i), rec->pos, gt, rid);
                        }
                        // shouldn't have to change bias for gt==0
                    } else if (rlen > alen && gt > 0) { // deletion
                        mi_writer.update(static_cast<size_t>(seq_start+bias+rec->pos), rec->pos, gt, rid);
                        mi_writer.update(static_cast<size_t>(seq_start+bias+rec->pos+1), rec->pos, gt, rid);
                        bias = bias - (rlen - 1);
                        // fprintf(stderr, "bias after del: %ld\n", bias);
                    } else if (rlen < alen && gt > 0) { // insertion
                        for (size_t i = 0; i < alen + 1; ++i) {
                            mi_writer.update(static_cast<size_t>(seq_start+bias+rec->pos+i), rec->pos, gt, rid);
                        }
                        bias = bias + alen - 1;
                        // fprintf(stderr, "bias after ins: %ld\n", bias);
                    }
                }
                if (args.ref_only) assert(gt==0);
                update_sequence(ref_seq, ref_len, rec, ppos_after, gt, fa_fp, log);
                ppos = rec->pos;
                ppos_after = ppos + strlen(rec->d.allele[0]);
                prange[0] = range[0];
                prange[1] = range[1];
            } else {
                fprintf(stderr, "Warning: skipping overlapping variant %d:%ld (%s). ", rec->rid, rec->pos, rec->d.id);
                fprintf(stderr, "To allow overlapping variants, use multiallelic lines: ie. $ bcftools norm -m+any <vcf> (dels after snps not supported)\n");
            }
        } else { // end of contig
            if (args.mai) {
                mi_writer.finish_sequence();
                //fprintf(stderr, "bias: %ld\n", bias);
            }
            update_sequence(ref_seq, ref_len, NULL, ppos_after, -1, fa_fp, log);
            ppos = 0;
            ppos_after = 0;
            seq_start += ref_len + args.w + bias;
            bias = 0;
            prange[0] = 0;
            prange[1] = 0;
            range[0] = 0;
            range[1] = 0;
        }
    };
    VCFScanner v(vargs);
    v.vcf_for_each(out_fn);
    fprintf(stderr, "done loop\n");
    if (args.mai) fclose(ma_fp);
    if (!args.ref_only) fclose(fa_fp);
}

int main(int argc, char** argv) {
    Args args(parse_args(argc, argv));
    if (args.ref_only) scan_vcf_sample(args, "");
    else for (auto s: args.samples) {
        scan_vcf_sample(args,s);
    }
    return 0;
}
