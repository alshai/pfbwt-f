#ifndef VCF_SCAN_HPP
#define VCF_SCAN_HPP

#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <vector>
#include <string>
extern "C" {
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/faidx.h>
}

class BCFGenotype {

    public:

    BCFGenotype() {};

    BCFGenotype(bcf_hdr_t* hdr, bcf1_t* rec) {
        bcf_unpack(rec, BCF_UN_FMT);
        n_ = bcf_get_genotypes(hdr, rec, &gt_arr_, &c_);
    }

    BCFGenotype(bcf_hdr_t* hdr, bcf1_t* rec, int32_t* g) : gt_arr_(g) {
        extern_buff = true;
        bcf_unpack(rec, BCF_UN_FMT);
        n_ = bcf_get_genotypes(hdr, rec, &gt_arr_, &c_);
    }

    size_t operator[](size_t i) {
        if (n_ < 0) return 0;
        return bcf_gt_allele(gt_arr_[i]);
    }

    size_t size() const { return n_; }
    size_t capacity() const {return c_; }

    ~BCFGenotype() {
        if (!extern_buff) free(gt_arr_);
    }

    private:

    int32_t* gt_arr_ = NULL;
    int32_t n_ = 0;
    int32_t c_ = 0;
    bool extern_buff = false;
};

struct VCFScannerArgs {
    std::vector<std::string> vcf_fnames;
    std::vector<std::string> contigs;
    std::string samples_fname;
    std::string sample;
    std::string ref_fasta;
    int wsize = 1;
    int verb = 0;
};

class VCFScanner {

    public:

    VCFScanner(VCFScannerArgs args)
        : vcf_fnames_(args.vcf_fnames)
        , samples_fname_(args.samples_fname)
        , contig_names_(args.contigs)
        , sample_(args.sample)
    {
        if (!vcf_fnames_.size()) {
            fprintf(stderr, "error: VCF files are required\n");
            exit(1);
        }
        // synced_reader helps us iterate by contig
        synced_readers_ = bcf_sr_init();
        bcf_sr_set_opt(synced_readers_, BCF_SR_REQUIRE_IDX);
        bcf_sr_set_opt(synced_readers_, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_SNPS, BCF_SR_PAIR_SNP_REF);
        for (auto vcf: vcf_fnames_)  {
            fprintf(stderr, "adding reader: %s\n", vcf.data());
            bcf_sr_add_reader(synced_readers_, vcf.data());
        }
        int err, err2 = 0;
        if (samples_fname_ != "") { // default to using the first name, if provided
            err = bcf_sr_set_samples(synced_readers_, samples_fname_.data(), 1);
            for (int i = 0; i < synced_readers_->nreaders; ++i) {
                bcf_hdr_t* hdr = synced_readers_->readers[i].header;
                err2 += bcf_hdr_set_samples(hdr, samples_fname_.data(), 1);
            }
        } else if (sample_ != "") {
            err = bcf_sr_set_samples(synced_readers_, sample_.data(), 0);
            for (int i = 0; i < synced_readers_->nreaders; ++i) {
                bcf_hdr_t* hdr = synced_readers_->readers[i].header;
                err2 += bcf_hdr_set_samples(hdr, sample_.data(), 0);
            }
        } else {
            fprintf(stderr, "warning - dropping all samples. outputting ref alleles only\n");
            err = bcf_sr_set_samples(synced_readers_, "-", 0);
            for (int i = 0; i < synced_readers_->nreaders; ++i) {
                bcf_hdr_t* hdr = synced_readers_->readers[i].header;
                err2 += bcf_hdr_set_samples(hdr, NULL, 0);
            }
        }
        if (!err) {
            fprintf(stderr, "failed to set samples in bcf_srs_t\n");
            exit(1);
        }
        if (err2) {
            fprintf(stderr, "failed to set samples in bcf_hdr_t\n");
            exit(1);
        }
        fprintf(stderr, "there are %d files ", synced_readers_->nreaders);
        fprintf(stderr, "subsetted to %d samples\n", synced_readers_->n_smpl);
        int nseqs = 0;
        const char** seqnames = bcf_hdr_seqnames(synced_readers_->readers[0].header, &nseqs);
        fprintf(stderr, "VCF describes %d contigs\n", nseqs);
        if (!contig_names_.size()) {
            fprintf(stderr, "using all contigs in VCF(s)\n");
            for (int i = 0; i < nseqs; ++i) {
                contig_names_.push_back(seqnames[i]);
            }
        } // TODO check for contigs that don't exist in VCF
        if (args.ref_fasta != "") {
            ref_faidx = fai_load(args.ref_fasta.data());  
        }
    }

    ~VCFScanner() {
        if (synced_readers_) bcf_sr_destroy(synced_readers_);
    }

    template<typename F>
    void vcf_for_each(F out_fn) {
        for (auto contig: contig_names_) {
            fprintf(stderr, "scanning contig \"%s\"\n", contig.data());
            vcf_contig_for_each(contig, out_fn);
        }
    }

    template<typename F>
    void vcf_contig_for_each(std::string contig, F out_fn) {
        // prep
        bcf_hdr_t* hdr = synced_readers_->readers[0].header;
        int id = bcf_hdr_name2id(hdr, contig.data());
        // size_t length = hdr->id[BCF_DT_CTG][id].val->info[0];
        bcf1_t* rec = bcf_init();
        bcf_sr_seek(synced_readers_, contig.data(), 0);
        int32_t* gt_buf = NULL;
        int32_t ppos = 0;
        std::vector<size_t> deltav(nhaplotypes());
        std::vector<size_t> posv(nhaplotypes());
        char* ref_seq = NULL;
        int ref_length = 0;
        if (ref_faidx != NULL) {
            ref_seq = fai_fetch(ref_faidx, bcf_hdr_id2name(hdr, id), &ref_length);
        }
        // iterate through records
        while (bcf_sr_next_line(synced_readers_)) {
            bcf1_t* rec = bcf_sr_get_line(synced_readers_, 0);
            if (rec->rid != id) break;
            int ref_len = strlen(rec->d.allele[0]);
            // if all samples droppped, gt_buf has size <= 0
            BCFGenotype gts(hdr, rec, gt_buf); // get genotypes for this line
            for (size_t i = 0; i < posv.size(); ++i ) { // update positions for this line
                posv[i] = posv[i] + (rec->pos - ppos) + deltav[i];
                int alt_len = strlen(rec->d.allele[gts[i]]);
                deltav[i] = alt_len - ref_len; // negative if deletion
            }
            out_fn(rec, gts, posv, ref_seq, ref_length);
            ppos = rec->pos;
        }
        BCFGenotype gts;
        out_fn(NULL, gts, posv, ref_seq, ref_length);
        free(gt_buf);
        bcf_destroy(rec);
        free(ref_seq);
    }

    void set_vcfs(std::vector<std::string> vs) {
        vcf_fnames_ = vs;
    }

    size_t nsamples() const {
        return static_cast<size_t>(synced_readers_->n_smpl);
    }

    // TODO: check ploidy then return. for now, assumes diploid
    size_t nhaplotypes() const {
        return static_cast<size_t>(synced_readers_->n_smpl * 2);
    }

    private:

    // parameters
    std::vector<std::string> vcf_fnames_;
    std::string samples_fname_;

    // vcf file handling stuff
    bcf_srs_t* synced_readers_;

    // useful info
    std::vector<std::string> contig_names_;
    std::vector<std::string> contigs_;
    std::string sample_;
    faidx_t* ref_faidx = NULL;
};

#endif
