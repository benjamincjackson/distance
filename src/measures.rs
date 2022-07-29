use crate::fastaio::EncodedFastaRecord;

// We can return this for all the distance-generating functions instead of 
// switching on whether they return a float or an integer measure
#[derive(Clone, Debug, PartialEq)]
pub enum FloatInt {
    Float(f64),
    Int(i64),
}

// Conventional snp-distance. Compares every site in the alignment. Might be faster than snp2
// for high diversity datasets.
// This is -m n_high in the CLI
pub fn snp(query: &EncodedFastaRecord, target: &EncodedFastaRecord) -> FloatInt {
    let mut d: i64 = 0;
    for i in 0..target.seq.len() {
        if query.seq[i] & target.seq[i] < 16 {
            d += 1;
        }
    }
    FloatInt::Int(d)
}

// Reduced snp-distance. Only compares sites that differ from the alignment(s)'s consensus in
// either record. Is fast in low-diversity datasets.
// This is -m n in the CLI
pub fn snp2(query: &EncodedFastaRecord, target: &EncodedFastaRecord) -> FloatInt {
    let mut d: i64 = 0;

    for idx in query.differences.iter() {
        if (query.seq[*idx] & target.seq[*idx]) < 16 {
            d += 1;
        }
    }

    let mut start = 0;
    for idx in target.differences.iter() {
        
        // if this site is different from the consensus in seq1 too, we've already tested it, so we must skip it here
        let result = query
                .differences[start..]
                .binary_search(idx);
        if result.is_ok() {
            start = result.unwrap(); // we can incrementally search a smaller slice of query.distances in future iterations if we find a match. 
            continue
        }

        // otherwise we check for a nucleotide difference
        if (query.seq[*idx] & target.seq[*idx]) < 16 {
            d += 1;
        }
    }
    
    FloatInt::Int(d)
}

// Number of nucleotide differences *per site*
pub fn raw(query: &EncodedFastaRecord, target: &EncodedFastaRecord) -> FloatInt {
    let mut d = 0.0;
    let mut n = 0.0;
    for i in 0..target.seq.len() {
        if query.seq[i] & 8 == 8 && query.seq[i] == target.seq[i] {
            d += 1.0;
        } else if query.seq[i] & target.seq[i] < 16 {
            d += 1.0;
            n += 1.0;
        }
    }

    FloatInt::Float(n / d)
}

// Jukes and Cantor's (1969) evolutionary distance
pub fn jc69(query: &EncodedFastaRecord, target: &EncodedFastaRecord) -> FloatInt {
    let temp = raw(query, target);
    let mut p: f64 = 0.0;
    match temp {
        FloatInt::Float(f) => p = f,
        _ => (),
    }

    FloatInt::Float(-0.75 * (1.0 - (4_f64 / 3_f64) * p).ln())
}

// Kimura's (1980) evolutionary distance
pub fn k80(query: &EncodedFastaRecord, target: &EncodedFastaRecord) -> FloatInt {
    let mut count_L: usize = 0;
    let mut ts: usize = 0;
    let mut tv: usize = 0;

    for i in 0..target.seq.len() {
        if (query.seq[i] & 8) == 8 && query.seq[i] == target.seq[i] { // are the bases certainly the same
            count_L += 1;
        } else if (query.seq[i] & target.seq[i]) < 16 { // they are certainly different
            if (query.seq[i] & 55) == 0 && (target.seq[i] & 55) == 0 { // both are purines, this is a transition
                ts += 1;
                count_L += 1;
            } else if (query.seq[i] & 199) == 0 && (target.seq[i] & 199) == 0 { // both are pyramidines, this is a transition
                ts += 1;
                count_L += 1;
            } else if ((query.seq[i] & 55) == 0 && (target.seq[i] & 199) == 0) 
                      || ((query.seq[i] & 199) == 0 && (target.seq[i] & 55) == 0) { // one of each, this is a transversion
                tv += 1; 
                count_L += 1;
            }
        }
    }

    let P = ts as f64 / count_L as f64;
    let Q = tv as f64 / count_L as f64;
    
    FloatInt::Float(-0.5 * ((1.0 - 2.0*P - Q) * (1.0 - 2.0*Q).sqrt()).ln())
}

// Tamura and Nei's (1993) evolutionary distance
pub fn tn93(query: &EncodedFastaRecord, target: &EncodedFastaRecord) -> FloatInt {
    // Total ATGC length of the two sequences
    let L: usize = query.count_A
        + query.count_T
        + query.count_G
        + query.count_C
        + target.count_A
        + target.count_T
        + target.count_G
        + target.count_C;

    // estimates of the equilibrium base contents from the pair's sequence data
    let g_A: f64 = (target.count_A as f64 + query.count_A as f64) / L as f64;
    let g_C: f64 = (target.count_C as f64 + query.count_C as f64) / L as f64;
    let g_G: f64 = (target.count_G as f64 + query.count_G as f64) / L as f64;
    let g_T: f64 = (target.count_T as f64 + query.count_T as f64) / L as f64;

    let g_R: f64 = (target.count_A as f64
        + query.count_A as f64
        + target.count_G as f64
        + query.count_G as f64)
        / L as f64;

    let g_Y: f64 = (target.count_C as f64
        + query.count_C as f64
        + target.count_T as f64
        + query.count_T as f64)
        / L as f64;

    // tidies up the equations a bit, after ape
    let k1: f64 = 2.0 * g_A * g_G / g_R;
    let k2: f64 = 2.0 * g_T * g_C / g_Y;
    let k3: f64 = 2.0 * (g_R * g_Y - g_A * g_G * g_Y / g_R - g_T * g_C * g_R / g_Y);

    let mut count_P1: usize = 0; // count of transitional differences between purines (A ⇄ G)
    let mut count_P2: usize = 0; // count of transitional differences between pyramidines (C ⇄ T)

    let mut count_d: usize = 0; // total number of differences
    let mut count_L: usize = 0; // total length of resolved comparison

    for i in 0..target.seq.len() {
        if query.seq[i] & 8 == 8 && query.seq[i] == target.seq[i] {
            // are the bases certainly the same
            count_L += 1;
        } else if (query.seq[i] & target.seq[i]) < 16
            && query.seq[i] & 8 == 8
            && target.seq[i] & 8 == 8
        {
            // are the bases different (and known for sure)
            count_d += 1;
            count_L += 1;
            if (query.seq[i] | target.seq[i]) == 200 {
                // 1 if one of the bases is adenine and the other one is guanine, 0 otherwise
                count_P1 += 1;
            } else if (query.seq[i] | target.seq[i]) == 56 {
                // 1 if one of the bases is cytosine and the other one is thymine, 0 otherwise
                count_P2 += 1;
            }
        }
    }

    // estimated rates from this pairwise comparison
    let P1: f64 = count_P1 as f64 / count_L as f64; // rate of changes which are transitional differences between purines (A ⇄ G)
    let P2: f64 = count_P2 as f64 / count_L as f64; // rate of changes which are transitional differences between pyramidines (C ⇄ T)
    let Q: f64 = (count_d - (count_P1 + count_P2)) as f64 / count_L as f64; // rate of changes which are transversional differences  (A ⇄ C || A ⇄ T || G ⇄ T || C ⇄ G) (i.e. everything else)

    // tidies up the equations a bit, after ape
    let w1: f64 = 1.0 - P1 / k1 - Q / (2.0 * g_R);
    let w2: f64 = 1.0 - P2 / k2 - Q / (2.0 * g_Y);
    let w3: f64 = 1.0 - Q / (2.0 * g_R * g_Y);

    let mut d = -k1 * w1.ln() - k2 * w2.ln() - k3 * w3.ln();
    if d == 0.0 {
        d = 0.0
    }

    FloatInt::Float(d)
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fasta::Record;
    use num_traits::Float;
    use crate::fastaio::*;

    #[test]
    fn test_snp() {
        let target_unencoded = Record::with_attrs("target", None, b"ATGATG");
        let target = encode(&target_unencoded).unwrap();

        let query_unencoded = Record::with_attrs("query", None, b"ATTATT");
        let query = encode(&query_unencoded).unwrap();

        let result = snp(&target, &query);
        assert_eq!(result, FloatInt::Int(2));
    }

    #[test]
    fn test_snp2() {
        let target_unencoded = Record::with_attrs("target", None, b"ATGATG");
        let target = encode(&target_unencoded).unwrap();

        let query_unencoded = Record::with_attrs("query", None, b"ATTATT");
        let query = encode(&query_unencoded).unwrap();
        
        let mut v = vec![vec![target, query]];
        let c = consensus(&v);

        v[0][0].get_differences(&c);
        v[0][1].get_differences(&c);

        let result = snp2(&v[0][0], &v[0][1]);
        assert_eq!(result, FloatInt::Int(2));
    }

    #[test]
    fn test_raw() {
        let target_unencoded = Record::with_attrs("target", None, b"ATGATG");
        let target = encode(&target_unencoded).unwrap();

        let query_unencoded = Record::with_attrs("query", None, b"ATTATT");
        let query = encode(&query_unencoded).unwrap();
        
        let v = vec![vec![target, query]];

        let result = raw(&v[0][0], &v[0][1]);
        assert_eq!(result, FloatInt::Float(2 as f64 / 6 as f64));
    }

    #[test]
    fn test_jc69() {
        let target_unencoded = Record::with_attrs("target", None, b"ATGATG");
        let target = encode(&target_unencoded).unwrap();

        let query_unencoded = Record::with_attrs("query", None, b"ATTATT");
        let query = encode(&query_unencoded).unwrap();
        
        let v = vec![vec![target, query]];

        let result = jc69(&v[0][0], &v[0][1]);
        assert_eq!(result, FloatInt::Float(-0.75 * (1.0 - (4_f64 / 3_f64) * (1.0 / 3.0)).ln()));
    }

    #[test]
    fn test_k80() {
        let target_unencoded = Record::with_attrs("target", None, b"ATGATG");
        let target = encode(&target_unencoded).unwrap();

        let query_unencoded = Record::with_attrs("query", None, b"ATTATT");
        let query = encode(&query_unencoded).unwrap();
        
        let v = vec![vec![target, query]];

        let result = k80(&v[0][0], &v[0][1]);

        let P = 0.0 / 6.0; // transitions
        let Q = 2.0 / 6.0; // transversions
        let desired_result = FloatInt::Float(-0.5 * ((1.0 - 2.0*P - Q) * (1.0 - 2.0*Q).sqrt()).ln());

        assert_eq!(result, desired_result);
    }

    #[test]
    fn test_tn93() {
        let target_unencoded = Record::with_attrs("target", None, b"ATGATGATGATGCCC");
        let target = encode(&target_unencoded).unwrap();

        let query_unencoded = Record::with_attrs("query", None, b"ATTATTATGATGCCC");
        let query = encode(&query_unencoded).unwrap();
        
        let mut v = vec![vec![target, query]];

        v[0][0].count_bases();
        v[0][1].count_bases();

        let result = tn93(&v[0][0], &v[0][1]);

        let g_A = 8.0 / 30.0;
        let g_T = 10.0 / 30.0;
        let g_C = 6.0 / 30.0;
        let g_G = 6.0 / 30.0;

        let g_R: f64 = (8.0 + 6.0) / 30.0;
        let g_Y: f64 = (7.0 + 9.0) / 30.0;

        let k1: f64 = 2.0 * g_A * g_G / g_R;
        let k2: f64 = 2.0 * g_T * g_C / g_Y;
        let k3: f64 = 2.0 * (g_R * g_Y - g_A * g_G * g_Y / g_R - g_T * g_C * g_R / g_Y);

        // estimated rates from this pairwise comparison
        let P1: f64 = 0.0 / 15.0; // rate of changes which are transitional differences between purines (A ⇄ G)
        let P2: f64 = 0.0 as f64 / 15.0; // rate of changes which are transitional differences between pyramidines (C ⇄ T)
        let Q: f64 = (2.0 - (0.0 + 0.0)) as f64 / 15.0; // rate of changes which are transversional differences  (A ⇄ C || A ⇄ T || G ⇄ T || C ⇄ G) (i.e. everything else)

        // tidies up the equations a bit, after ape
        let w1: f64 = 1.0 - P1 / k1 - Q / (2.0 * g_R);
        let w2: f64 = 1.0 - P2 / k2 - Q / (2.0 * g_Y);
        let w3: f64 = 1.0 - Q / (2.0 * g_R * g_Y);

        let d = -k1 * w1.ln() - k2 * w2.ln() - k3 * w3.ln();
        let desired_result = FloatInt::Float(d);

        assert_eq!(result, desired_result);

    }
}