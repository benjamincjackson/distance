use crate::fastaio::EncodedFastaRecord;

pub fn tn93(query: &EncodedFastaRecord, target: &EncodedFastaRecord) -> f64 {
	
	// Total ATGC length of the two sequences
	let L: usize = query.count_A + query.count_T + query.count_G + query.count_C + target.count_A + target.count_T + target.count_G + target.count_C;

	// estimates of the equilibrium base contents from the pair's sequence data
	let g_A: f64 = (target.count_A as f64 + query.count_A as f64) / L as f64;
	let g_C: f64 = (target.count_C as f64 + query.count_C as f64) / L as f64;
	let g_G: f64 = (target.count_G as f64 + query.count_G as f64) / L as f64;
	let g_T: f64 = (target.count_T as f64 + query.count_T as f64) / L as f64;

	let g_R: f64 = (target.count_A as f64 + query.count_A as f64 + target.count_G as f64 + query.count_G as f64) / L as f64;
	let g_Y: f64 =  (target.count_C as f64 + query.count_C as f64 + target.count_T as f64 + query.count_T as f64) / L as f64;

	// tidies up the equations a bit, after ape
	let k1: f64 = 2.0 * g_A * g_G / g_R;
	let k2: f64 = 2.0 * g_T * g_C / g_Y;
	let k3: f64 = 2.0 * (g_R*g_Y - g_A*g_G*g_Y/g_R - g_T*g_C*g_R/g_Y);

	let mut count_P1: usize = 0; // count of transitional differences between purines (A ⇄ G)
	let mut count_P2: usize = 0; // count of transitional differences between pyramidines (C ⇄ T)

	let mut count_d: usize = 0; // total number of differences
	let mut count_L: usize = 0; // total length of resolved comparison

	for i in 0..target.seq.len() {
		if query.seq[i]&8 == 8 && query.seq[i] == target.seq[i] { // at the bases certainly the same
			count_L += 1;
		} else if (query.seq[i]&target.seq[i]) < 16 && query.seq[i]&8 == 8 && target.seq[i]&8 == 8 { // are the bases different (and known for sure)
			count_d += 1;
			count_L += 1;
			if (query.seq[i] | target.seq[i]) == 200 { // 1 if one of the bases is adenine and the other one is guanine, 0 otherwise
				count_P1 += 1;
			} else if (query.seq[i] |  target.seq[i]) == 56 { // 1 if one of the bases is cytosine and the other one is thymine, 0 otherwise
				count_P2 += 1;
			}
		} 
	}

	// estimated rates from this pairwise comparison
	let P1: f64 = count_P1 as f64 / count_L as f64;				   // rate of changes which are transitional differences between purines (A ⇄ G)
	let P2: f64 = count_P2 as f64 / count_L as f64;				   // rate of changes which are transitional differences between pyramidines (C ⇄ T)
	let Q: f64 = (count_d-(count_P1+count_P2)) as f64 / count_L as f64;  // rate of changes which are transversional differences  (A ⇄ C || A ⇄ T || C ⇄ A || C ⇄ G) (i.e. everything else)

	// tidies up the equations a bit, after ape
	let w1: f64 = 1.0 - P1/k1 - Q/(2.0*g_R);
	let w2:f64 = 1.0 - P2/k2 - Q/(2.0*g_Y);
	let w3:f64 = 1.0 - Q/(2.0*g_R*g_Y);

	-k1 * w1.ln() - k2 * w2.ln() - k3 * w3.ln() // d!
}