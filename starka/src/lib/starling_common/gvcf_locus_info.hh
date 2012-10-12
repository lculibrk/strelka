// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file
///
/// \author Chris Saunders
///

#ifndef __GVCF_LOCUS_INFO_HH
#define __GVCF_LOCUS_INFO_HH


#include "blt_common/position_snp_call_pprob_digt.hh"
#include "starling_common/align_path.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"
#include "starling_common/starling_shared.hh"

#include <bitset>
#include <iosfwd>


namespace VCF_FILTERS {

    enum index_t {
        HighDepth,
        LowGQX,
        HighSB,
        HighHpol,
        HighBaseFilt,
        HighRepeatCount,
        IndelConflict,
        SiteConflict,
        SIZE
    };

    inline 
    const char*
    get_label(const unsigned idx) {
        switch(idx) {
        case HighDepth: return "HighDepth";
        case LowGQX: return "LowGQX";
        case HighSB: return "HighSB";
        case HighHpol: return "HighHpol";
        case HighBaseFilt: return "HighBaseFilt";
        case HighRepeatCount: return "HighRepeatCount";
        case IndelConflict: return "IndelConflict";
        case SiteConflict: return "SiteConflict";
        default:
            assert(0);
            return NULL;
        }
    }
}



struct shared_modifiers {
    
    shared_modifiers() { clear(); }

    void
    set_filter(const VCF_FILTERS::index_t i) {
        filters.set(i);
    }

    void
    write_filters(std::ostream& os) const;

    void
    clear() {
        filters.reset();
    }

    int gqx;
    std::bitset<VCF_FILTERS::SIZE> filters;
};


struct indel_modifiers : public shared_modifiers {
    indel_modifiers() { clear(); }

    void
    clear() {
        shared_modifiers::clear();
        is_overlap=false;
    }

    ALIGNPATH::path_t cigar;

    bool is_overlap;
};

struct site_modifiers : public shared_modifiers {

    site_modifiers() { clear(); }

    void
    clear() {
        shared_modifiers::clear();
        is_block=false;
    }

    bool
    is_gqx() const {
        return ((!is_unknown) && is_used_covered);
    }

    bool is_unknown;
    bool is_covered;
    bool is_used_covered;
    bool is_block;

    unsigned max_gt;
};



struct indel_info {

    void
    init(const pos_t init_pos,
         const indel_key& init_ik,
         const starling_diploid_indel_core& init_dindel,
         const starling_indel_report_info& init_iri,
         const starling_indel_sample_report_info& init_isri)
    { 
        pos=(init_pos);
        ik=(init_ik);
        dindel=(init_dindel);
        iri=(init_iri);
        isri=(init_isri);
        imod.clear();
    }

    const char*
    get_gt() {
        if(imod.is_overlap) {
            return "1/2";
        }
        return STAR_DIINDEL::get_gt_label(dindel.max_gt);
    }

    pos_t pos;
    indel_key ik;
    starling_diploid_indel_core dindel;
    starling_indel_report_info iri;
    starling_indel_sample_report_info isri;

    indel_modifiers imod;
};



struct site_info {

    site_info()
        : pos(0)
        , ref('N')
        , n_used_calls(0)
        , n_unused_calls(0)
        , hpol(0)
    {}

    void
    init(const pos_t init_pos,
         const char init_ref,
         const snp_pos_info& good_pi,
         const bool used_allele_count_min_qscore) {

        pos=(init_pos);
        ref=(init_ref);
        good_pi.get_known_counts(known_counts,used_allele_count_min_qscore);
        smod.clear();
    }

    const char*
    get_gt() const {
        if(smod.is_unknown || (!smod.is_used_covered)) {
            return "./.";
        } else {
            return DIGT::get_vcf_gt(smod.max_gt,dgt.ref_gt);
        }
    }

    pos_t pos;
    char ref;
    unsigned n_used_calls;
    unsigned n_unused_calls;
    boost::array<unsigned,N_BASE> known_counts;
    diploid_genotype dgt;
    unsigned hpol;

    site_modifiers smod;
};



#endif
