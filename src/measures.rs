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
    let mut d = 0u32;
    // u8 counters let the compiler use narrow SIMD lanes; 255 is the most one can hold without overflow
    for (qc, tc) in query.seq.chunks(255).zip(target.seq.chunks(255)) {
        let mut dc = 0u8;
        for (&a, &b) in qc.iter().zip(tc) {
            dc += ((a & b) < 16) as u8;
        }
        d += dc as u32;
    }

    FloatInt::Int(d as i64)
}

// Reduced snp-distance. Only compares sites that differ from the alignment(s)'s consensus in
// either record. Is fast in low-diversity datasets.
// This is -m n in the CLI
pub fn snp_consensus(query: &EncodedFastaRecord, target: &EncodedFastaRecord) -> FloatInt {
    let mut d: i64 = 0;

    for idx in query.differences.iter() {
        if (query.seq[*idx] & target.seq[*idx]) < 16 {
            d += 1;
        }
    }

    let mut start = 0;
    for idx in target.differences.iter() {
        // if this site is different from the consensus in seq1 too, we've already tested it, so we must skip it here
        let result = query.differences[start..].binary_search(idx);
        if let Ok(pos) = result {
            start = pos; // we can incrementally search a smaller slice of query.distances in future iterations if we find a match.
            continue;
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
    let mut d = 0u32;
    let mut n = 0u32;
    for (qc, tc) in query.seq.chunks(255).zip(target.seq.chunks(255)) {
        let mut dc = 0u8;
        let mut nc = 0u8;
        for (&a, &b) in qc.iter().zip(tc) {
            let same = ((a & 8 != 0) & (a == b)) as u8;
            let disjoint = ((a & b) < 16) as u8;
            dc += same + disjoint;
            nc += disjoint;
        }
        d += dc as u32;
        n += nc as u32;
    }

    FloatInt::Float(n as f64 / d as f64)
}

// Jukes and Cantor's (1969) evolutionary distance
pub fn jc69(query: &EncodedFastaRecord, target: &EncodedFastaRecord) -> FloatInt {
    let temp = raw(query, target);
    let mut p: f64 = 0.0;
    if let FloatInt::Float(f) = temp {
        p = f
    }
    FloatInt::Float(-0.75 * (1.0 - (4.0 / 3.0) * p).ln())
}

// Returns (count_L, ts, tv) for k80
fn k80_counts(query: &EncodedFastaRecord, target: &EncodedFastaRecord) -> (usize, usize, usize) {
    let mut count_L: usize = 0;
    let mut ts: usize = 0;
    let mut tv: usize = 0;

    // Branch-free per-site classification; see `raw` for the chunking rationale.
    // Within a chunk each of these counters grows by at most 1 per site, so u8 is safe.
    for (qc, tc) in query.seq.chunks(255).zip(target.seq.chunks(255)) {
        let mut lc = 0u8;
        let mut tsc = 0u8;
        let mut tvc = 0u8;
        for (&a, &b) in qc.iter().zip(tc) {
            // the bases are certainly the same
            let same = ((a & 8 != 0) & (a == b)) as u8;
            // the bases are certainly different
            let disjoint = (a & b) < 16;

            let a_pur = a & 55 == 0;
            let b_pur = b & 55 == 0;
            let a_pyr = a & 199 == 0;
            let b_pyr = b & 199 == 0;

            // both purines or both pyrimidines: transition
            let is_ts = disjoint & ((a_pur & b_pur) | (a_pyr & b_pyr));
            // one of each: transversion
            let is_tv = disjoint & !is_ts & ((a_pur & b_pyr) | (a_pyr & b_pur));

            lc += same + is_ts as u8 + is_tv as u8;
            tsc += is_ts as u8;
            tvc += is_tv as u8;
        }
        count_L += lc as usize;
        ts += tsc as usize;
        tv += tvc as usize;
    }

    (count_L, ts, tv)
}

// Kimura's (1980) evolutionary distance
pub fn k80(query: &EncodedFastaRecord, target: &EncodedFastaRecord) -> FloatInt {
    let (count_L, ts, tv) = k80_counts(query, target);

    let P = ts as f64 / count_L as f64;
    let Q = tv as f64 / count_L as f64;

    FloatInt::Float(-0.5 * ((1.0 - 2.0 * P - Q) * (1.0 - 2.0 * Q).sqrt()).ln())
}

// Returns (count_L, count_d, count_P1, count_P2) for tn93:
// count_L:  total length of resolved comparison
// count_d:  total number of differences
// count_P1: transitional differences between purines (A <-> G)
// count_P2: transitional differences between pyrimidines (C <-> T)
fn tn93_counts(
    query: &EncodedFastaRecord,
    target: &EncodedFastaRecord,
) -> (usize, usize, usize, usize) {
    let mut count_P1: usize = 0;
    let mut count_P2: usize = 0;
    let mut count_d: usize = 0;
    let mut count_L: usize = 0;

    // Branch-free per-site classification; see `raw` for the chunking rationale.
    // Within a chunk each of these counters grows by at most 1 per site, so u8 is safe.
    for (qc, tc) in query.seq.chunks(255).zip(target.seq.chunks(255)) {
        let mut lc = 0u8;
        let mut dc = 0u8;
        let mut p1c = 0u8;
        let mut p2c = 0u8;
        for (&a, &b) in qc.iter().zip(tc) {
            // the bases are certainly the same
            let same = ((a & 8 != 0) & (a == b)) as u8;
            // the bases are different (and known for sure)
            let diff = (((a & b) < 16) & (a & 8 != 0) & (b & 8 != 0)) as u8;
            // one base is adenine and the other guanine
            let p1 = diff & ((a | b) == 200) as u8;
            // one base is cytosine and the other thymine
            let p2 = diff & ((a | b) == 56) as u8;

            lc += same + diff;
            dc += diff;
            p1c += p1;
            p2c += p2;
        }
        count_L += lc as usize;
        count_d += dc as usize;
        count_P1 += p1c as usize;
        count_P2 += p2c as usize;
    }

    (count_L, count_d, count_P1, count_P2)
}

// Tamura and Nei's (1993) evolutionary distance
pub fn tn93(query: &EncodedFastaRecord, target: &EncodedFastaRecord) -> FloatInt {
    let query_base_counts = query.base_counts.clone().unwrap();
    let target_base_counts = target.base_counts.clone().unwrap();

    // Total ATGC length of the two sequences
    let L: usize = query_base_counts.A
        + query_base_counts.T
        + query_base_counts.G
        + query_base_counts.C
        + target_base_counts.A
        + target_base_counts.T
        + target_base_counts.G
        + target_base_counts.C;

    // estimates of the equilibrium base contents from the pair's sequence data
    let g_A: f64 = (target_base_counts.A as f64 + query_base_counts.A as f64) / L as f64;
    let g_C: f64 = (target_base_counts.C as f64 + query_base_counts.C as f64) / L as f64;
    let g_G: f64 = (target_base_counts.G as f64 + query_base_counts.G as f64) / L as f64;
    let g_T: f64 = (target_base_counts.T as f64 + query_base_counts.T as f64) / L as f64;

    let g_R: f64 = (target_base_counts.A as f64
        + query_base_counts.A as f64
        + target_base_counts.G as f64
        + query_base_counts.G as f64)
        / L as f64;

    let g_Y: f64 = (target_base_counts.C as f64
        + query_base_counts.C as f64
        + target_base_counts.T as f64
        + query_base_counts.T as f64)
        / L as f64;

    // tidies up the equations a bit, after ape
    let k1: f64 = 2.0 * g_A * g_G / g_R;
    let k2: f64 = 2.0 * g_T * g_C / g_Y;
    let k3: f64 = 2.0 * (g_R * g_Y - g_A * g_G * g_Y / g_R - g_T * g_C * g_R / g_Y);

    let (count_L, count_d, count_P1, count_P2) = tn93_counts(query, target);

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
    use crate::fastaio::*;
    use bio::io::fasta::Reader;
    use num_traits::Float;

    const TARGET_FASTA: &[u8] = b">target
ATGATGATGATGCCC
";

    const QUERY_FASTA: &[u8] = b">query
ATTATTATGATGCCC
";

    fn setup() -> (EncodedFastaRecord, EncodedFastaRecord) {
        let target_reader = Reader::new(TARGET_FASTA);
        let target = encode(&target_reader.records().next().unwrap().unwrap()).unwrap();
        let query_reader = Reader::new(QUERY_FASTA);
        let query = encode(&query_reader.records().next().unwrap().unwrap()).unwrap();

        (target, query)
    }

    // Straightforward branching reference implementations, used to check that the
    // branch-free versions give identical counts.
    fn ref_snp(q: &[u8], t: &[u8]) -> i64 {
        let mut d = 0;
        for i in 0..t.len() {
            if q[i] & t[i] < 16 {
                d += 1;
            }
        }
        d
    }

    fn ref_raw(q: &[u8], t: &[u8]) -> (usize, usize) {
        let (mut d, mut n) = (0, 0);
        for i in 0..t.len() {
            if q[i] & 8 == 8 && q[i] == t[i] {
                d += 1;
            } else if q[i] & t[i] < 16 {
                d += 1;
                n += 1;
            }
        }
        (d, n)
    }

    fn ref_k80(q: &[u8], t: &[u8]) -> (usize, usize, usize) {
        let (mut count_l, mut ts, mut tv) = (0, 0, 0);
        for i in 0..t.len() {
            if (q[i] & 8) == 8 && q[i] == t[i] {
                count_l += 1;
            } else if (q[i] & t[i]) < 16 {
                if (q[i] & 55) == 0 && (t[i] & 55) == 0 {
                    ts += 1;
                    count_l += 1;
                } else if (q[i] & 199) == 0 && (t[i] & 199) == 0 {
                    ts += 1;
                    count_l += 1;
                } else if ((q[i] & 55) == 0 && (t[i] & 199) == 0)
                    || ((q[i] & 199) == 0 && (t[i] & 55) == 0)
                {
                    tv += 1;
                    count_l += 1;
                }
            }
        }
        (count_l, ts, tv)
    }

    fn ref_tn93(q: &[u8], t: &[u8]) -> (usize, usize, usize, usize) {
        let (mut count_l, mut count_d, mut p1, mut p2) = (0, 0, 0, 0);
        for i in 0..t.len() {
            if q[i] & 8 == 8 && q[i] == t[i] {
                count_l += 1;
            } else if (q[i] & t[i]) < 16 && q[i] & 8 == 8 && t[i] & 8 == 8 {
                count_d += 1;
                count_l += 1;
                if (q[i] | t[i]) == 200 {
                    p1 += 1;
                } else if (q[i] | t[i]) == 56 {
                    p2 += 1;
                }
            }
        }
        (count_l, count_d, p1, p2)
    }

    fn record(seq: Vec<u8>) -> EncodedFastaRecord {
        let mut r = EncodedFastaRecord::new();
        r.seq = seq;
        r
    }

    // Every valid encoded value (see encoding.rs)
    const CODES: [u8; 17] = [
        136, 72, 40, 24, 192, 160, 144, 96, 80, 48, 224, 176, 208, 112, 240, 244, 242,
    ];

    // Sequence pairs covering every combination of codes, a long run of identical
    // sites (to exercise the 255-site chunk boundaries), and pseudo-random data.
    fn test_pairs() -> Vec<(Vec<u8>, Vec<u8>)> {
        let mut pairs = Vec::new();

        let mut q = Vec::new();
        let mut t = Vec::new();
        for &a in CODES.iter() {
            for &b in CODES.iter() {
                q.push(a);
                t.push(b);
            }
        }
        pairs.push((q, t));

        for len in [0, 1, 254, 255, 256, 510, 511, 1000] {
            pairs.push((vec![136; len], vec![136; len]));
            pairs.push((vec![136; len], vec![24; len]));
            pairs.push((vec![136; len], vec![72; len]));
            pairs.push((vec![136; len], vec![40; len]));
        }

        let mut state: u64 = 0x2545F4914F6CDD1D;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 32) as usize
        };
        for len in [100, 255, 1000, 2551] {
            let q: Vec<u8> = (0..len).map(|_| CODES[next() % CODES.len()]).collect();
            let t: Vec<u8> = (0..len).map(|_| CODES[next() % CODES.len()]).collect();
            pairs.push((q, t));
        }

        pairs
    }

    #[test]
    fn test_snp_matches_reference() {
        for (q, t) in test_pairs() {
            let expected = ref_snp(&q, &t);
            assert_eq!(snp(&record(q), &record(t)), FloatInt::Int(expected));
        }
    }

    #[test]
    fn test_raw_matches_reference() {
        for (q, t) in test_pairs() {
            let (d, n) = ref_raw(&q, &t);
            let expected = n as f64 / d as f64;
            match raw(&record(q), &record(t)) {
                FloatInt::Float(f) => assert!(f == expected || (f.is_nan() && expected.is_nan())),
                _ => panic!("raw should return a float"),
            }
        }
    }

    #[test]
    fn test_k80_counts_match_reference() {
        for (q, t) in test_pairs() {
            assert_eq!(
                k80_counts(&record(q.clone()), &record(t.clone())),
                ref_k80(&q, &t)
            );
        }
    }

    #[test]
    fn test_tn93_counts_match_reference() {
        for (q, t) in test_pairs() {
            assert_eq!(
                tn93_counts(&record(q.clone()), &record(t.clone())),
                ref_tn93(&q, &t)
            );
        }
    }

    #[test]
    fn test_snp() {
        let (target, query) = setup();
        let result = snp(&target, &query);
        assert_eq!(result, FloatInt::Int(2));
    }

    #[test]
    fn test_snp_consensus() {
        let (target, query) = setup();

        let mut v = vec![vec![target, query]];
        let c = consensus(&v);

        v[0][0].get_differences(&c);
        v[0][1].get_differences(&c);

        let result = snp_consensus(&v[0][0], &v[0][1]);
        assert_eq!(result, FloatInt::Int(2));
    }

    #[test]
    fn test_raw() {
        let (target, query) = setup();
        let result = raw(&query, &target);
        assert_eq!(result, FloatInt::Float(2.0 / 15.0));
    }

    #[test]
    fn test_jc69() {
        let (target, query) = setup();
        let result = jc69(&query, &target);
        assert_eq!(
            result,
            FloatInt::Float(-0.75 * (1.0 - (4.0 / 3.0) * (2.0 / 15.0)).ln())
        );
    }

    #[test]
    fn test_k80() {
        let (target, query) = setup();

        let result = k80(&query, &target);

        let P = 0.0 / 15.0; // transitions
        let Q = 2.0 / 15.0; // transversions
        let desired_result =
            FloatInt::Float(-0.5 * ((1.0 - 2.0 * P - Q) * (1.0 - 2.0 * Q).sqrt()).ln());

        assert_eq!(result, desired_result);
    }

    #[test]
    fn test_tn93() {
        let (target, query) = setup();

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
