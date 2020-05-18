#include <deque>
#include <vector>
#include <string>
#include <getopt.h>
#include "marker_array/vcf_scanner.hpp"
#include "marker_array/marker_index.hpp"

struct Args {
    int w = 10;
    std::string fa_fname;
    std::vector<std::string> vcf_fnames = {"/home/taher/r-index2/pfbwt-f/data/chr21.snps.vcf.gz"};
    std::vector<std::string> contigs = {"21"};
    std::vector<std::string> samples = {"HG00096"};
    std::string ref_fasta = "/home/taher/pfbwt-f/data/chr21.fa";
    std::string out = "out";
    int nthreads = 1;
    int verb = 0;
    int haplotype = 0;
    int ref_only = 0;
    int to_stdout = 0;
};

VCFScannerArgs ArgsToVCFScannerArgs(Args args) {
    VCFScannerArgs v;
    v.vcf_fnames = args.vcf_fnames;
    v.contigs = args.contigs;
    v.verb = args.verb;
    v.wsize = args.w;
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
        {"output", required_argument, NULL, 'o'},
        {"threads", required_argument, NULL, 't'},
        {"contigs", required_argument, NULL, 'c'},
        {"samples", required_argument, NULL, 'S'},
        // {"sample", required_argument, NULL, 's'},
        {"verbose", no_argument, NULL, 'v'},
        {"fasta", required_argument, NULL, 'f'},
        {"ref-only", no_argument, NULL, 'r'},
        {"stdout", no_argument, &args.to_stdout, 1}
    };

    while ((c = getopt_long( argc, argv, "rf:w:o:t:c:S:vH:", lopts, NULL) ) != -1) {
        switch(c) {
            case 'f':
                args.ref_fasta = optarg; break;
            // case 's':
            //     args.samples = {optarg}; break;
            case 'w':
                args.w = atoi(optarg); break;
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
    int ppos = 0;
    int i = args.haplotype;
    std::string fa_fname = args.out + "." + sample + "." + std::to_string(i) + ".fa";
    std::string fa_header = args.out + "." + sample + "." + std::to_string(i);
    std::string ma_fname = args.out + "." + sample + "." + std::to_string(i) + ".mai";
    std::string log_fname = args.out + "." + sample + "." + std::to_string(i) + ".log";
    if (args.ref_only) {
        ma_fname = args.out + ".ref.mai";
        log_fname = args.out + ".ref.log";
        fa_fname = args.out + ".ref.fa";
    }
    log = fopen(log_fname.data(), "w");
    ma_fp = fopen(ma_fname.data(), "wb");
    if (args.to_stdout) {
        fa_fp = stdout;
    } else {
        fa_fp = fopen(fa_fname.data(), "w");
    }
    MarkerIndexWriter mi_writer(args.w, ma_fp, args.verb ? NULL : log);
    int ppos_after = 0;
    std::string pseq("");
    // auto out_fn = [&](bcf1_t* rec, BCFGenotype& gtv, std::vector<size_t>& posv, char* ref_seq, int ref_len) {
    auto out_fn = [&](bcf1_t* rec, BCFGenotype& gtv, std::vector<size_t>& posv, char* ref_seq, int32_t ref_len, const char* ref_name) {
        if (strcmp(ref_name, pseq.data())) {
            std::string to_write(ref_name);
            if (!args.ref_only) {
                to_write = sample + "." + std::to_string(args.haplotype) + "." + to_write;
            }
            fprintf(fa_fp, ">%s\n", to_write.data());
            pseq.assign(ref_name);
        }

        if (rec != NULL) {
            if (rec->pos != ppos) {
                int pos = args.ref_only ? rec->pos : posv[i];
                int gt =  args.ref_only ? 0        : gtv[i];
                mi_writer.update(rec->pos, gt, pos);
                update_sequence(ref_seq, ref_len, rec, ppos_after, gt, fa_fp, log);
                ppos = rec->pos;
                ppos_after = ppos + strlen(rec->d.allele[0]);
            } else fprintf(stderr, "warning: overlapping variants at %d. skipping... \n", rec->pos);
        } else {
            mi_writer.update(ref_len, -1, ref_len);
            update_sequence(ref_seq, ref_len, NULL, ppos_after, -1, fa_fp, log);
        }
    };
    VCFScanner v(vargs);
    v.vcf_for_each(out_fn);
    fprintf(stderr, "done loop\n");
    fclose(ma_fp);
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
