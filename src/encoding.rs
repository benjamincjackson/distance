// An array whose indices are integer representations of characters and 
// who contents are Emmanual Paradis' bitwise coding scheme. Used to encode
// the sequences for fast comparisons
pub fn encoding_array() -> [u8; 256] {
    let mut a: [u8; 256] = [0; 256];

    a['A' as usize] = 136;
    a['a' as usize] = 136;
    a['G' as usize] = 72;
    a['g' as usize] = 72;
    a['C' as usize] = 40;
    a['c' as usize] = 40;
    a['T' as usize] = 24;
    a['t' as usize] = 24;
    a['R' as usize] = 192;
    a['r' as usize] = 192;
    a['M' as usize] = 160;
    a['m' as usize] = 160;
    a['W' as usize] = 144;
    a['w' as usize] = 144;
    a['S' as usize] = 96;
    a['s' as usize] = 96;
    a['K' as usize] = 80;
    a['k' as usize] = 80;
    a['Y' as usize] = 48;
    a['y' as usize] = 48;
    a['V' as usize] = 224;
    a['v' as usize] = 224;
    a['H' as usize] = 176;
    a['h' as usize] = 176;
    a['D' as usize] = 208;
    a['d' as usize] = 208;
    a['B' as usize] = 112;
    a['b' as usize] = 112;
    a['N' as usize] = 240;
    a['n' as usize] = 240;
    a['-' as usize] = 244;
    a['?' as usize] = 242;

    a
}
