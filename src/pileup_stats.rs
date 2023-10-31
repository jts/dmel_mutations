//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
use bio::stats::{Prob, PHREDProb};

pub struct PileupStats
{
    // this is indexed by base (dim four), strand (two)
    base_counts: [u32; 8],
    pub mean_qual: [PHREDProb; 3],
}

impl PileupStats {
    #[inline(always)]
    pub fn get(&self, b: u32, s: u32) -> u32 {
        let i = (s * 4) + b;
        return self.base_counts[ i as usize ];
    }

    pub fn increment(&mut self, b: u32, s: u32) -> () {
        let i = (s * 4) + b;
        self.base_counts[i as usize] += 1;
    }

    pub fn new() -> PileupStats {
        let ps = PileupStats { base_counts: [0; 8], mean_qual: [PHREDProb(0.0); 3] };
        return ps;
    }

    pub fn clear(& mut self) -> () {
        self.base_counts = [0; 8];
        self.mean_qual = [PHREDProb(0.0); 3];
    }

    pub fn fill_pileup(& mut self, alignments: rust_htslib::bam::pileup::Alignments<'_>) -> () {
        self.clear();

        let mut sum_p = [ 0.0, 0.0, 0.0 ];
        let mut count = [ 0, 0, 0];
        for a in alignments {
            if a.record().seq().len() == 0 {
                continue;
            }

            if let Some(qpos) = a.qpos() {
                let read_base = a.record().seq()[qpos] as char;
                let bi = base2index(read_base) as u32;
                let si = a.record().is_reverse() as u32;
                self.increment(bi, si);
                
                let qual = PHREDProb(a.record().qual()[qpos] as f64);
                sum_p[si as usize] += *Prob::from(qual);
                count[si as usize] += 1;

                // count both strands
                sum_p[2] += *Prob::from(qual);
                count[2] += 1;
            }
        }

        for si in 0..=2 {
            let mean_p = sum_p[si]/ count[si] as f64;
            self.mean_qual[si] = PHREDProb::from(Prob(mean_p));
        }
    }
}

pub fn base2index(base: char) -> i32 {
    return match base {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        _ => -1
    };
}

