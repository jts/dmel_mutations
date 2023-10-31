//---------------------------------------------------------
// Copyright 2023 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
extern crate clap;
use std::fs::File;
use std::collections::HashMap;
use clap::{Arg, App, SubCommand, value_t};
use rust_htslib::{bam, bam::Read, bam::record::Aux, bam::ext::BamRecordExtensions};
use rust_htslib::{bcf::Reader as BcfReader, bcf::Read as BcfRead};
use bio::stats::{Prob, PHREDProb};
use bio::alphabets::dna;
use csv::ReaderBuilder;
use rand::Rng;
use itertools::Itertools;

mod pileup_stats;
mod annotation;
use crate::annotation::MitoAnnotation;
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
                .arg(Arg::with_name("annotation")
                    .short("a")
                    .long("annotation")
                    .required(true)
                    .takes_value(true)
                    .help("annotation VCF"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .subcommand(SubCommand::with_name("read-statistics")
                .about("calculate read-level summary")
                .arg(Arg::with_name("genome")
                    .short("g")
                    .long("genome")
                    .required(true)
                    .takes_value(true)
                    .help("the reference genome"))
                .arg(Arg::with_name("annotation")
                    .short("a")
                    .long("annotation")
                    .required(true)
                    .takes_value(true)
                    .help("annotation VCF"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .subcommand(SubCommand::with_name("run-statistics")
                .about("calculate run-level statistics")
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .subcommand(SubCommand::with_name("read-matrix")
                .about("calculate mutation matrix for each read")
                .arg(Arg::with_name("genome")
                    .short("g")
                    .long("genome")
                    .required(true)
                    .takes_value(true)
                    .help("the reference genome"))
                .arg(Arg::with_name("sites")
                    .short("s")
                    .long("sites")
                    .required(true)
                    .takes_value(true)
                    .help("the sites to use to calculate APOBEC genotype"))
                .arg(Arg::with_name("output-style")
                    .short("o")
                    .long("output-style")
                    .required(false)
                    .takes_value(true)
                    .help("the output format"))
                .arg(Arg::with_name("max-reads")
                    .short("m")
                    .long("max-reads")
                    .required(false)
                    .takes_value(true)
                    .help("the maximum number of reads to process"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .subcommand(SubCommand::with_name("apobec-edit-distance")
                .about("calculate edit distance between apobec sites")
                .arg(Arg::with_name("alleles")
                    .short("a")
                    .long("alleles")
                    .required(true)
                    .takes_value(true)
                    .help("file containing alleles for each read"))
                .arg(Arg::with_name("from")
                    .short("a")
                    .long("from")
                    .required(true)
                    .takes_value(true)
                    .help("either a filename with the duplex read names or 'random'")))
        .subcommand(SubCommand::with_name("apobec-vcf")
                .about("generate a VCF file containing putative APOBEC mutations")
                .arg(Arg::with_name("input-fasta")
                    .required(true)
                    .index(1)
                    .help("the input fasta file to process")))
        .subcommand(SubCommand::with_name("process-vcf")
                .about("convert VEP's VCF format to a more usable format")
                .arg(Arg::with_name("input-bam")
                    .short("b")
                    .long("bam")
                    .required(true)
                    .takes_value(true)
                    .help("bam file containing chrM-to-chrM mapping"))
                .arg(Arg::with_name("input-vcf")
                    .required(true)
                    .index(1)
                    .help("the input vcf file to process")))
        .get_matches();

    if let Some(matches) = matches.subcommand_matches("pileup-frequency") {
        calculate_pileup_frequency(matches.value_of("input-bam").unwrap(), matches.value_of("genome").unwrap(), matches.value_of("annotation").unwrap());
    }

    if let Some(matches) = matches.subcommand_matches("read-statistics") {
        calculate_read_statistics(matches.value_of("input-bam").unwrap(), matches.value_of("genome").unwrap(), matches.value_of("annotation").unwrap());
    } 
    
    if let Some(matches) = matches.subcommand_matches("run-statistics") {
        calculate_run_statistics(matches.value_of("input-bam").unwrap());
    } 
    
    if let Some(matches) = matches.subcommand_matches("read-matrix") {
        let max_reads = value_t!(matches, "max-reads", usize).unwrap_or(1000000);
        read_to_matrix(matches.value_of("input-bam").unwrap(), 
                       matches.value_of("genome").unwrap(), 
                       matches.value_of("sites").unwrap(), 
                       matches.value_of("output-style").unwrap_or("string"),
                       max_reads);
    }
    
    if let Some(matches) = matches.subcommand_matches("apobec-edit-distance") {
        apobec_edit_distance(matches.value_of("alleles").unwrap(), matches.value_of("from").unwrap());
    }
    
    if let Some(matches) = matches.subcommand_matches("apobec-vcf") {
        apobec_vcf(matches.value_of("input-fasta").unwrap());
    }
    
    if let Some(matches) = matches.subcommand_matches("process-vcf") {
        process_vcf(matches.value_of("input-vcf").unwrap(), matches.value_of("input-bam").unwrap());
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

fn is_apobec(chromosome: &[u8], pos: usize) -> bool {
    let is_apobec = if pos > 0 && pos < chromosome.len() - 3 {
        let fwd_context = String::from_utf8_lossy(& chromosome[(pos - 1)..(pos + 1)] ).into_owned();
        let rev_context = String::from_utf8( dna::revcomp( &chromosome[pos..(pos + 2)] ) ).unwrap();
        fwd_context == "TC" || rev_context == "TC"
    } else {
        false
    };
    return is_apobec;
}

fn calculate_pileup_frequency(input_bam: &str, reference_genome: &str, annotation_vcf: &str) {
    let annotation = MitoAnnotation::from_vcf(annotation_vcf);
    let mut bam = bam::IndexedReader::from_path(input_bam).expect("Could not read input bam file:");
    let header_view = bam.header().clone();
    let max_depth = 1000000;

    let chromosome_name = "chrM";
    let tid = header_view.tid(chromosome_name.as_bytes()).unwrap();
    let chromosome_bytes = get_chromosome_sequence(reference_genome, &header_view, tid).into_bytes();
    let output_strand = [ /*"+", "-",*/ "both" ];

    bam.fetch(chromosome_name).unwrap();

    println!("chromosome\tposition\tbase\tcontext\tis_apobec\tstrand\tcsq\tdepth\terror_rate\tobs_qual\tmean_base_qual\tcount_A\tcount_C\tcount_G\tcount_T");
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
            let pos1 = pileup.pos() + 1;
            let class = annotation.get_class(pos1 as usize);
            println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{}", 
                    contig, pos1, ref_base as char, context, is_apobec as u8, os,
                    class, total, diff_rate, qual, *mq,
                    counts[0], counts[1], counts[2], counts[3]);
        }
    }
}

fn calculate_read_statistics(input_bam: &str, reference_genome: &str, annotation_vcf: &str) {
    let annotation = MitoAnnotation::from_vcf(annotation_vcf);

    let mut bam = bam::IndexedReader::from_path(input_bam).expect("Could not read input bam file:");
    let header_view = bam.header().clone();

    let chromosome_name = "chrM";
    let tid = header_view.tid(chromosome_name.as_bytes()).unwrap();
    let chromosome_bytes = get_chromosome_sequence(reference_genome, &header_view, tid).into_bytes();

    let mutation_classes = MitoAnnotation::output_order();

    let mut header_fields = Vec::from( [ "read_name",
                                     "chromosome",
                                     "position",
                                     "total_aligned",
                                     "total_mismatch",
                                     "mismatch_rate",
                                     "apobec_aligned",
                                     "apobec_mismatch",
                                     "apobec_rate" ] );

    header_fields.extend_from_slice(mutation_classes.as_slice());
    let header_str = header_fields.iter().format("\t");
    println!("{header_str}");

    bam.fetch(chromosome_name).unwrap();
    for r in bam.records() {
        let record = r.unwrap();
        let contig = String::from_utf8_lossy(header_view.tid2name(record.tid() as u32));
        
        if record.seq().len() == 0 {
            continue;
        }
        let mut total_match = 0;
        let mut total_aligned = 0;

        let mut apobec_match = 0;
        let mut apobec_aligned = 0;

        let mut apobec_variant_classes = HashMap::new();
        for k in &mutation_classes { apobec_variant_classes.insert(String::from(*k), 0); }

        for ap in record.aligned_pairs() {
            
            let qpos = ap[0];
            let rpos = ap[1];
            let read_base = record.seq()[qpos as usize];
            let ref_base = chromosome_bytes[rpos as usize];

            total_match += (read_base == ref_base) as usize;
            total_aligned += 1;

            let is_apobec = is_apobec(&chromosome_bytes, rpos as usize);
            if is_apobec {

                let mutated_base = if ref_base as char == 'C' {
                    'T' as u8
                } else {
                    'A' as u8
                };

                if read_base == ref_base || read_base == mutated_base {
                    let is_match = (read_base == ref_base) as usize;

                    apobec_match += is_match;
                    apobec_aligned += 1;

                    if read_base == mutated_base {

                        let rpos1 = rpos + 1; // 1-based coord for lookup
                        let class = annotation.get_class(rpos1 as usize);
                        if let Some(x) = apobec_variant_classes.get_mut(&class) {
                            *x += 1;
                        }
                    }
                }
            }
        }

        let qname = std::str::from_utf8(record.qname()).unwrap();
        let mismatch_rate = (total_aligned - total_match) as f64 / total_aligned as f64;
        let apobec_rate = (apobec_aligned - apobec_match) as f64 / apobec_aligned as f64;

        let class_counts: Vec<_> = mutation_classes.iter().map(|x| apobec_variant_classes.get(&String::from(*x)).unwrap() ).collect();
        let class_str = class_counts.iter().format("\t");

        let f = format!("{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{:.3}\t{}", 
            qname, contig, record.pos(), total_aligned, total_aligned - total_match, mismatch_rate, apobec_aligned, apobec_aligned - apobec_match, apobec_rate, class_str);
        println!("{f}");
    }
}

fn calculate_run_statistics(input_bam: &str) {
    let mut bam = bam::Reader::from_path(input_bam).expect("Could not read input bam file:");
    let header_view = bam.header().clone();

    let mut total_reads = 0;
    let mut total_duplex = 0;
    let mut total_bases = 0;
    
    let mut full_length_mt_reads = 0;
    let mut full_length_mt_reads_duplex = 0;

    for r in bam.records() {
        let record = r.unwrap();
     
        if record.is_secondary() || record.is_supplementary() || record.seq().len() == 0 {
            continue;
        }

        total_reads += 1;
        total_bases += record.seq().len();
        let mut duplex_tag = 0;
        if let Ok(value) = record.aux(b"dx") {
            if let Aux::U8(v) = value {
                duplex_tag = v as usize;
            }
        }
        total_duplex += duplex_tag;

        if record.is_unmapped() { 
            continue;
        }

            
        let contig = String::from_utf8_lossy(header_view.tid2name(record.tid() as u32));
        if contig == "chrM" && record.pos() <= 50 && record.seq().len() > 18500 {
            full_length_mt_reads += 1;
            full_length_mt_reads_duplex += duplex_tag;
        }
    }

    println!("filename\ttotal_bases\ttotal_reads\ttotal_full_length_mtdna\tmtdna_percent\tmtdna_duplex_percent\tnon_mtdna_duplex_percent");
    
    let percent_mtdna = (full_length_mt_reads as f64 / total_reads as f64) * 100.0;

    let mt_duplex_percent = (full_length_mt_reads_duplex as f64 / full_length_mt_reads as f64) * 100.0;
    let non_mt_duplex_percent = ((total_duplex - full_length_mt_reads_duplex) as f64 / (total_reads - full_length_mt_reads) as f64) * 100.0;
    let total_mb = total_bases as f64 / 1000000.0;
    println!("{input_bam}\t{total_mb:.1}\t{total_reads}\t{full_length_mt_reads}\t{percent_mtdna:.2}\t{mt_duplex_percent:.2}\t{non_mt_duplex_percent:.2}");

}

fn read_to_matrix(input_bam: &str, reference_genome: &str, typing_path: &str, output_format: &str, max_reads: usize) {
    let mut bam = bam::IndexedReader::from_path(input_bam).expect("Could not read input bam file:");
    let header_view = bam.header().clone();

    let typing_file = File::open(typing_path).unwrap();
    let mut rdr = ReaderBuilder::new()
                       .delimiter(b'\t')
                       .from_reader(typing_file);
    let mut typing_positions = Vec::new();
    let mut typing_map = HashMap::new();
    let mut idx = 0;
    for result in rdr.records() {
        let record = result.unwrap();
        let position = record.get(1).unwrap().parse::<usize>().unwrap();
        typing_map.insert(position, idx);
        typing_positions.push(position);
        idx += 1;
    }
    
    let min_base_quality = 20.0;

    let chromosome_name = "chrM";
    let tid = header_view.tid(chromosome_name.as_bytes()).unwrap();
    let chromosome_bytes = get_chromosome_sequence(reference_genome, &header_view, tid).into_bytes();

    bam.fetch(chromosome_name).unwrap();

    let mut read_idx = 0;
    if output_format == "string" {
        println!("read_idx\tread_name\tchromosome\talleles");
    } else {
        println!("read_idx\tread_name\tchromosome\tposition\tallele");
    }

    for r in bam.records() {
        let record = r.unwrap();
        let contig = String::from_utf8_lossy(header_view.tid2name(record.tid() as u32));
        
        if record.seq().len() == 0 {
            continue;
        }

        let mut apobec_type = vec![b'N'; typing_map.len()];
        for ap in record.aligned_pairs() {
            
            let qpos = ap[0];
            let rpos = ap[1];
            let read_base = record.seq()[qpos as usize];
            let read_base_qual = record.seq()[qpos as usize];
            let ref_base = chromosome_bytes[rpos as usize];
            let qual = PHREDProb(read_base_qual as f64);

            if *qual < min_base_quality { continue; }

            let is_apobec = is_apobec(&chromosome_bytes, rpos as usize);
            if is_apobec {

                let mutated_base = if ref_base as char == 'C' {
                    'T' as u8
                } else {
                    'A' as u8
                };

                if read_base == ref_base || read_base == mutated_base {
                    let is_match = (read_base == ref_base) as usize;

                    if let Some(type_idx) = typing_map.get( &(rpos as usize) ) {
                        apobec_type[*type_idx as usize] = (1 - is_match) as u8 + 48;
                    }
                }
            }
        }

        let qname = std::str::from_utf8(record.qname()).unwrap();
        
        if output_format == "string" {
            let apobec_str = apobec_type.iter().map(|x| (*x as char).to_string()).collect::<Vec<String>>().join("");
            println!("{}\t{}\t{}\t{}", read_idx, qname, contig, apobec_str);
        } else {
            for idx in 0..apobec_type.len() {
              println!("{}\t{}\t{}\t{}\t{}", read_idx, qname, contig, idx, apobec_type[idx] as char);
            }
        }
        
        read_idx += 1;
        if read_idx > max_reads { break; }
    }
}

fn apobec_edit_distance(read_matrix_tsv: &str, from: &str) {
    let matrix_file = File::open(read_matrix_tsv).unwrap();
    let mut rdr = ReaderBuilder::new()
                       .delimiter(b'\t')
                       .from_reader(matrix_file);
    let mut read_names = Vec::new();
    let mut read_map = HashMap::new();
    for result in rdr.records() {
        let record = result.unwrap();
        let name = record.get(1).unwrap().parse::<String>().unwrap();
        let alleles = record.get(3).unwrap().parse::<String>().unwrap();
        read_map.insert(name.clone(), alleles.into_bytes());
        read_names.push(name);
    }
    let mut n = read_names.len();

    let mut pairs = Vec::new();
    if from == "random" {
        let mut rng = rand::thread_rng();
        for _i in 0..1000 {
            let ai = rng.gen_range(0..n);
            let bi = rng.gen_range(0..n);
            let an = &read_names[ai];
            let bn = &read_names[bi];
            pairs.push( (an.clone(), bn.clone()) );
        }
    } else if from == "allpairs" {
        n = 500;
        for i in 0..n {
            for j in (i+1)..n {
                pairs.push( (read_names[i].clone(), read_names[j].clone()) );
            }
        }
    } else {
        let pairs_file = File::open(from).unwrap();
        let mut pair_rdr = ReaderBuilder::new()
                              .has_headers(false)
                              .delimiter(b'\t')
                              .from_reader(pairs_file);
        for result in pair_rdr.records() {
            let record = result.unwrap();
            let an = record.get(0).unwrap().parse::<String>().unwrap();
            let bn = record.get(1).unwrap().parse::<String>().unwrap();
            pairs.push( (an, bn) );
        }
    };
    
    println!("read_name_a\tread_name_b\ttotal_sites\tapobec_edit_distance");
    for (an, bn) in pairs {

        if ! (read_map.contains_key(&an) && read_map.contains_key(&bn)) {
            continue;
        }
        let a = read_map.get(&an).unwrap();
        let b = read_map.get(&bn).unwrap();

        let mut t = 0;
        let mut d = 0;
        for j in 0..a.len() {
            if a[j] == b'N' || b[j] == b'N' {
                continue;
            }

            d += (a[j] != b[j]) as usize;
            t += 1;
        }

        println!("{an}\t{bn}\t{t}\t{d}")
    }
}

fn apobec_vcf(fasta: &str) {
    let faidx = rust_htslib::faidx::Reader::from_path(fasta).expect("Could not read reference genome:");
    let chromosome_name = "chrM";
    let mut chromosome_sequence = faidx.fetch_seq_string(&chromosome_name, 0, 1000000).unwrap();
    chromosome_sequence.make_ascii_uppercase();
    let chromosome_bytes = chromosome_sequence.clone().into_bytes();

    println!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    for rpos in 0..chromosome_sequence.len() {
        let ref_base = chromosome_bytes[rpos as usize];
        let is_apobec = is_apobec(&chromosome_bytes, rpos as usize);
        if is_apobec {

            let mutated_base = if ref_base as char == 'C' {
                'T' as u8
            } else {
                'A' as u8
            };

            println!("{chromosome_name}\t{}\t.\t{}\t{}\t.\t.\t.", rpos + 1, ref_base as char, mutated_base as char);
        }
    }
}

fn process_vcf(vcf: &str, bam: &str) {

    // build a map from dm6 mtDNA position to our assembled reference
    let mut bam = bam::Reader::from_path(bam).expect("Could not read input bam file:");
    let mut position_map = HashMap::new();
    for r in bam.records() {
        let record = r.unwrap();
        for ap in record.aligned_pairs() {
            let qpos = ap[0];
            let rpos = ap[1];
            position_map.insert(rpos, qpos);   
        }
    }

    // Maps that collapse different VEP consequences into a reduced set of classes
    let csq_map = HashMap::from([
            ("upstream_gene_variant", "noncoding"),
            ("downstream_gene_variant", "noncoding"),
            ("incomplete_terminal_codon_variant&coding_sequence_variant", "nonsyn"),
            ("missense_variant", "nonsyn"),
            ("non_coding_transcript_exon_variant", "tRNA"),
            ("start_lost", "lof"),
            ("stop_gained", "lof"),
            ("stop_retained_variant", "syn"),
            ("synonymous_variant", "syn")
    ]);

    // Define the precedence (for selection) for each consequence class
    let precedence_map = MitoAnnotation::precedence_map();

    let mut bcf = BcfReader::from_path(vcf).expect("Error opening file.");
    println!("##fileformat=VCFv4.1");
    println!("##INFO=<ID=assembly_position,Number=1,Type=Integer,Description=\"Position of variant on chrM assembly\">");
    println!("##INFO=<ID=CSQClass,Number=1,Type=String,Description=\"Variant consequence type\">");
    println!("##INFO=<ID=CSQGene,Number=1,Type=String,Description=\"Annotated gene for variant\">");
    println!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    for record_result in bcf.records() {
        let record = record_result.expect("Fail to read record");
        let csq_array = record.info(b"CSQ").string().unwrap().unwrap();

        let mut selected_class = None;
        let mut max_precedence = precedence_map.len();
        let mut selected_gene = None;
        for idx in 0..csq_array.len() {
            
            let parts = std::str::from_utf8(csq_array[idx]).unwrap().split("|").collect::<Vec<&str>>();
            let csq = parts[1];
            let gene = parts[3];
            let class = csq_map.get(csq).unwrap();
            let precedence = precedence_map.get(class).unwrap();

            if precedence < &max_precedence {
                max_precedence = *precedence;
                selected_class = Some(class);
                selected_gene = Some(gene);
            }
        }

        let dm6_pos = record.pos(); // 0 based
        let assembly_position1 = if position_map.contains_key(&dm6_pos) { (*position_map.get(&dm6_pos).unwrap()) + 1 } else { -1 };

        // output, convert to 1-based
        println!("{}\t{}\t.\t{}\t{}\t.\t.\tassembly_position={};CSQClass={};CSQGene={}", 
               "chrM", record.pos() + 1, record.alleles()[0][0] as char, record.alleles()[1][0] as char, assembly_position1, selected_class.unwrap(), selected_gene.unwrap());
    }
}

