use std::env;
extern crate needletail;
use needletail::{parse_fastx_file, Sequence};

fn main() {
    let args: Vec<String> = env::args().collect();
    let filepath = &args[1];
    let pattern_string = &args[2];
    let pattern = pattern_string.as_bytes();
    let patter_length = pattern.len() as u8;
    
    let mut n_bases = 0;
    let mut n_valid_kmers = 0;
    let mut reader = parse_fastx_file(&filepath).expect("valid path/file");
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        // keep track of the total number of bases
        n_bases += seqrec.num_bases();
        // normalize to make sure all the bases are consistently capitalized and
        // that we remove the newlines since this is FASTA
        let norm_seq = seqrec.normalize(false);
        // we make a reverse complemented copy of the sequence first for
        // `canonical_kmers` to draw the complemented sequences from.
        let rc = norm_seq.reverse_complement();
        // now we keep track of the number of AAAAs (or TTTTs via
        // canonicalization) in the file; note we also get the position (i.0;
        // in the event there were `N`-containing kmers that were skipped)
        // and whether the sequence was complemented (i.2) in addition to
        // the canonical kmer (i.1)
        for (_, kmer, _) in norm_seq.canonical_kmers(patter_length, &rc) {
            if kmer == pattern {
                n_valid_kmers += 1;
            }
        }
    }
    println!("There are {} bases in your file.", n_bases);
    println!("There are {n_valid_kmers} {pattern_string} in your file.");
}
