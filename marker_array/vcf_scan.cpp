#include <deque>
#include <vector>
#include <string>
#include <getopt.h>
#include "marker_array/vcf_scanner.hpp"

struct Args {
    int w = 10;
    std::string fa_fname;
    std::vector<std::string> vcf_fnames = {"/home/taher/r-index2/pfbwt-f/data/chr21.snps.trunc.vcf.gz"};
    std::vector<std::string> contigs = {"21"};
    std::vector<std::string> samples = {"HG00096"};
    std::string ref_fasta = "/home/taher/pfbwt-f/data/chr21.fa";
    std::string out = "out";
    int nthreads = 2;
    int verb = 0;
    int haplotype = 0;
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
        {"verbose", required_argument, NULL, 'v'}
    };

    while ((c = getopt_long( argc, argv, "w:o:t:c:S:vH:", lopts, NULL) ) != -1) {
        switch(c) {
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

struct Marker {
    Marker() {}
    Marker(size_t p, size_t r, size_t a) : textpos(p), refpos(r), allele(a) {}
    size_t textpos;
    size_t refpos;
    size_t allele;
};

class MarkerWindow : public std::deque<Marker> {

    using MDeq = std::deque<Marker>;
    public:

    using iterator = MDeq::iterator;

    MarkerWindow() {}
    MarkerWindow(size_t w) : w_(w) {}

    // returns end() if last element is still within window
    // otherwise, return iterator to last element within begin()->textpos+w
    iterator end_of_window() {
        size_t spos = MDeq::begin()->textpos;
        for (auto it = MDeq::begin(); it != MDeq::end() ; ++it) {
            if (it->textpos >= spos + w_) {
                return it;
            }
        }
        return MDeq::end();
    }

    void push_back(size_t p, size_t r, size_t a) {
        std::deque<Marker>::push_back(Marker(p,r,a));
    }

    size_t get_pos() const { return pos_; }
    void set_pos(size_t p) { pos_ = p; }
    int get_w() const { return w_; }
    void set_w(int w) { w_ = w; }

    private:

    size_t pos_ = 0;
    int w_ = 10;
};

constexpr size_t zero = 0;
constexpr size_t negone = -1;

void update_marker_index(size_t rpos, int32_t gt, size_t mpos, MarkerWindow& win, FILE* ofp, FILE* log) {
    if (gt > -1) win.push_back(mpos, rpos, gt);
    size_t pos = win.get_pos();
    // skip pos ahead if needed
    if (pos < win.front().textpos - win.get_w() + 1) {
        pos = win.front().textpos - win.get_w() + 1;
        win.set_pos(pos);
    }
    // figure out if window is completed
    if (gt > -1) {
        typename MarkerWindow::iterator win_end;
        while ((win_end = win.end_of_window()) != win.end()) {
            // process window every position upto and including pos
            for (; pos <= win.front().textpos; ++pos) {
                win.set_pos(pos);
                fwrite(&pos, sizeof(size_t), 1, ofp);
                fprintf(log, "%lu: ", pos);
                for (auto it = win.begin(); it != win_end; ++it) {
                    if (pos <= it->textpos && it->textpos < pos + win.get_w()) {
                        fprintf(log, "(%lu %lu %lu) ", it->textpos, it->refpos, it->allele);
                        fwrite(&(it->refpos), sizeof(size_t), 1, ofp);
                        fwrite(&(it->allele), sizeof(size_t), 1, ofp);
                    } else break;
                }
                fwrite(&negone, sizeof(size_t), 1, ofp);
                fprintf(log, "\n");
            }
            win.pop_front();
        }
    } else { // last rec in file
        for (; pos <= win.front().textpos; ++pos) {
            win.set_pos(pos);
            fwrite(&pos, sizeof(size_t), 1, ofp);
            fprintf(log, "%lu: ", pos);
            for (auto it = win.begin(); it != win.end(); ++it) {
                if (pos <= it->textpos && it->textpos < pos + win.get_w()) {
                    fprintf(log, "(%lu %lu %lu) ", it->textpos, it->refpos, it->allele);
                    fwrite(&(it->refpos), sizeof(size_t), 1, ofp);
                    fwrite(&(it->allele), sizeof(size_t), 1, ofp);
                } else break;
            }
            fwrite(&negone, sizeof(size_t), 1, ofp);
            fprintf(log, "\n");
        }
    }
}

void update_sequence(char* seq, int rlen, bcf1_t* rec, int ppos, int32_t gt, FILE* ofp, FILE* log) {
    (void) log;
    if (rec != NULL && gt > -1)  {
        fwrite(seq+ppos, sizeof(char), rec->pos - ppos, ofp); // write ref sequence
        char* allele = rec->d.allele[gt];
        fwrite(allele, sizeof(char), strlen(allele), ofp); // add allele
    } else {
        fwrite(seq+ppos, sizeof(char), rlen - ppos, ofp); // write ref sequence
    }
}

// produces fasta and marker array for two haplotypes of a single simple
void scan_vcf_sample(Args args, std::string sample) {
    VCFScannerArgs vargs(ArgsToVCFScannerArgs(args));
    vargs.sample = sample;
    FILE *log, *ma_fp, *fa_fp;
    MarkerWindow winv;
    int ppos = 0;
    int i = args.haplotype;
    std::string fa_fname = args.out + "." + sample + "." + std::to_string(i) + ".fa";
    std::string fa_header = args.out + "." + sample + "." + std::to_string(i);
    std::string ma_fname = args.out + "." + sample + "." + std::to_string(i) + ".ma2";
    std::string log_fname = args.out + "." + sample + "." + std::to_string(i) + ".log";
    log = fopen(log_fname.data(), "w");
    ma_fp = fopen(ma_fname.data(), "wb");
    fa_fp = fopen(fa_fname.data(), "w");
    fprintf(fa_fp, ">%s\n", fa_header.data());
    winv.set_w(args.w);
    int ppos_after = 0;
    auto out_fn = [&](bcf1_t* rec, BCFGenotype& gtv, std::vector<size_t>& posv, char* ref_seq, int ref_len) {
        if (rec != NULL && rec->pos != ppos) { // default cause
            update_marker_index(rec->pos, gtv[i], posv[i], winv, ma_fp, log);
            update_sequence(ref_seq, ref_len, rec, ppos_after, gtv.size() ? gtv[i] : -1, fa_fp, log);
            ppos = rec->pos;
            ppos_after = ppos + strlen(rec->d.allele[0]);
        } else if (rec == NULL) { // last case
            update_marker_index(ref_len, -1, posv[i], winv, ma_fp, log);
            update_sequence(ref_seq, ref_len, NULL, ppos_after,  -1, fa_fp, log);
        } else {
            fprintf(stderr, "warning: overlapping variants at %d. skipping... \n", rec->pos); 
        }
    };
    VCFScanner v(vargs);
    v.vcf_for_each(out_fn);
    fprintf(stderr, "done loop\n");
    fclose(ma_fp);
    fclose(fa_fp);
}

int main(int argc, char** argv) {
    (void) argc;
    (void) argv;
    Args args;
    for (auto s: args.samples) {
        scan_vcf_sample(args,s);
    }
    return 0;
}
