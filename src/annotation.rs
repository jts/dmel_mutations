use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufRead};

pub struct VariantAnnotation
{
    pub dm6_position: usize,
    pub class: String
}

pub struct MitoAnnotation
{
    pub map: HashMap<usize, VariantAnnotation>
}

impl MitoAnnotation
{
    pub fn from_vcf(vcf: &str) -> MitoAnnotation {
        let file = File::open(vcf).unwrap();
        let reader = BufReader::new(file);
        let mut mt = MitoAnnotation { map: HashMap::new() };

        for line in reader.lines() {
            let record = line.unwrap();
            if record.chars().nth(0).unwrap() == '#' { continue; }
            let fields: Vec<_> = record.split("\t").collect();
            //if fields[0][0] == "#" { continue; } // skip headers
            let info_map : HashMap<String, String> = HashMap::from_iter(fields[7].split(";").map(|a| {
                    let (key, value) = a.split_once("=").unwrap();
                        (key.to_string(), value.to_string())
                }));
            
            if info_map["assembly_position"] == "-1" { continue; }
            let assembly_position = info_map["assembly_position"].parse::<usize>().unwrap();
            let csq = info_map["CSQClass"].clone();
            //println!("{} {}", info_map["assembly_position"], csq);
            
            let va = VariantAnnotation {
                dm6_position: fields[1].parse::<usize>().unwrap(),
                class: csq
            };
            
            mt.map.insert(assembly_position as usize, va);
        }
        mt
    }

    pub fn get_class(self: &MitoAnnotation, position1: usize) -> String {
        if let Some(va) = self.map.get(&(position1)) {
            return va.class.clone();
        } else {
            return "unknown".to_string();
        }
    }

    pub fn precedence_map() -> HashMap<&'static str, usize> {
        // Define the precedence (for selection) for each consequence class
        return HashMap::from([
            ("lof", 0),
            ("nonsyn", 1),
            ("tRNA", 2),
            ("syn", 3),
            ("noncoding", 4),
            ("misc", 5)
        ]);
    }

    pub fn output_order() -> Vec<&'static str> {
        return Vec::from([ "lof", "nonsyn", "tRNA", "syn", "noncoding", "misc", "unknown" ]);
    }
}
