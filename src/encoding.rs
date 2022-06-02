// to do - just make the (immutable) table straight away
pub fn encoding_array() -> [u8; 256] {
    let mut a: [u8; 256] = [0; 256];

    a['A' as usize] = 136; // A
    a['a' as usize] = 136; // a
    a['G' as usize] = 72; // G
    a['g' as usize] = 72; // g
    a['C' as usize] = 40; // C
    a['c' as usize] = 40; // c
    a['T' as usize] = 24; // T
    a['t' as usize] = 24; // t
    a['R' as usize] = 192; // R
    a['r' as usize] = 192; // r
    a['M' as usize] = 160; // M
    a['m' as usize] = 160; // m
    a['W' as usize] = 144; // W
    a['w' as usize] = 144; // w
    a['S' as usize] = 96; // S
    a['s' as usize] = 96; // s
    a['K' as usize] = 80; // K
    a['k' as usize] = 80; // k
    a['Y' as usize] = 48; // Y
    a['y' as usize] = 48; // y
    a['V' as usize] = 224; // V
    a['v' as usize] = 224; // v
    a['H' as usize] = 176; // H
    a['h' as usize] = 176; // h
    a['D' as usize] = 208; // D
    a['d' as usize] = 208; // d
    a['B' as usize] = 112; // B
    a['b' as usize] = 112; // b
    a['N' as usize] = 240; // N
    a['n' as usize] = 240; // n
    a['-' as usize] = 244; // -
    a['?' as usize] = 242; // ?

    a
}
