//---------------------------------------------------------
// Copyright 2023 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
extern crate clap;
use clap::{Arg, App, SubCommand};
use rust_htslib::{bam, bam::Read };
use bio::stats::{Prob, PHREDProb};
use bio::alphabets::dna;

mod pileup_stats;
use crate::pileup_stats::PileupStats;
use crate::pileup_stats::base2index;

fn main() {
    let matches = App::new("dmel-mutations")
        .version("0.1")
        .author("Jared Simpson <jared.simpson@oicr.on.ca>")
        .subcommand(SubCommand::with_name("pileup-frequency")
                .about("calculate the frequency of observed bases per position")
                .arg(Arg::with_name("genome")
                    .short("g")
                    .long("genome")
                    .required(true)
                    .takes_value(true)
                    .help("the reference genome"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .get_matches();

    if let Some(matches) = matches.subcommand_matches("pileup-frequency") {
        calculate_pileup_frequency(matches.value_of("input-bam").unwrap(), matches.value_of("genome").unwrap());
    }
}

pub fn get_chromosome_sequence(reference_genome: &str,
                               bam_header: &rust_htslib::bam::HeaderView,
                               tid: u32) -> String {

    let faidx = rust_htslib::faidx::Reader::from_path(reference_genome).expect("Could not read reference genome:");
    let chromosome_length = bam_header.target_len(tid).unwrap() as usize;
    let chromosome_name = std::str::from_utf8(bam_header.tid2name(tid)).unwrap();

    let mut chromosome_sequence = faidx.fetch_seq_string(&chromosome_name, 0, chromosome_length).unwrap();
    chromosome_sequence.make_ascii_uppercase();
    return chromosome_sequence;
}

fn calculate_pileup_frequency(input_bam: &str, reference_genome: &str) {
    let mut bam = bam::IndexedReader::from_path(input_bam).expect("Could not read input bam file:");
    let header_view = bam.header().clone();
    let max_depth = 1000000;

    let chromosome_name = "chrM";
    let tid = header_view.tid(chromosome_name.as_bytes()).unwrap();
    let chromosome_bytes = get_chromosome_sequence(reference_genome, &header_view, tid).into_bytes();
    let output_strand = [ /*"+", "-",*/ "both" ];

    bam.fetch(chromosome_name).unwrap();

    println!("chromosome\tposition\tbase\tcontext\tis_apobec\tstrand\tdepth\terror_rate\tobs_qual\tmean_base_qual\tcount_A\tcount_C\tcount_G\tcount_T");
    let mut ps = PileupStats::new();
    let mut pileups = bam.pileup();
    pileups.set_max_depth(max_depth);

    for p in pileups {
        let pileup = p.unwrap();
        let contig = String::from_utf8_lossy(header_view.tid2name(pileup.tid() as u32));
        let pos = pileup.pos() as usize;
        let ref_base = chromosome_bytes[pos];
        
        //println!("{}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());
        ps.fill_pileup(pileup.alignments());

        let context = if pos > 0 && pos < chromosome_bytes.len() - 3 {
            String::from_utf8_lossy(& chromosome_bytes[(pos - 1)..(pos + 2)] ).into_owned()
        } else {
                "-".to_string()
        };

        let is_apobec = if pos > 0 && pos < chromosome_bytes.len() - 3 {
            let fwd_context = String::from_utf8_lossy(& chromosome_bytes[(pos - 1)..(pos + 1)] ).into_owned();
            let rev_context = String::from_utf8( dna::revcomp( &chromosome_bytes[pos..(pos + 2)] ) ).unwrap();
            fwd_context == "TC" || rev_context == "TC"
        } else {
            false
        };

        for os in output_strand {
            let mut counts = [ 0, 0, 0, 0 ];
            if os == "both" || os == "+" {
                for i in 0..4 { counts[i] += ps.get(i as u32, 0) }
            }

            if os == "both" || os == "-" {
                for i in 0..4 { counts[i] += ps.get(i as u32, 1) }
            }

            let mq = match os {
                "+" => ps.mean_qual[0],
                "-" => ps.mean_qual[1],
                "both" => ps.mean_qual[2],
                _ => PHREDProb(0.0)
            };

            let total = counts[0] + counts[1] + counts[2] + counts[3];
            let total_ref = counts[base2index(ref_base as char) as usize];
            let diff_rate = (total - total_ref) as f64 / total as f64;

            let qual = *PHREDProb::from(Prob(diff_rate));
            println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{}", 
                    contig, pileup.pos(), ref_base as char, context, is_apobec as u8, os,
                    total, diff_rate, qual, *mq,
                    counts[0], counts[1], counts[2], counts[3]);
        }
    }
}

