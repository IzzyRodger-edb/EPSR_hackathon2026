//! Random number generator — translation of `ran1`, `initran1` and `saveran1`
//! from `auxiliary_routines.f` (common block `/ran1com/` in `ran1com.inc`).
//!
//! This is the Park–Miller "minimal standard" generator with a Bays–Durham
//! shuffle (Numerical Recipes `ran1`).  Its full state (`mcall`, `idum`, `iy`,
//! `iv(32)`) is stored at the end of every `.ato` file, so a simulation can be
//! continued — and reproduced — exactly.  The single-precision arithmetic of the
//! original (`AM = 1./IM` in `real*4`) is preserved.

pub const IA: i32 = 16807;
pub const IM: i32 = 2147483647;
pub const IQ: i32 = 127773;
pub const IR: i32 = 2836;
pub const NTAB: usize = 32;
pub const NDIV: i32 = 1 + (IM - 1) / NTAB as i32;

/// `AM = 1./IM` evaluated in single precision.
fn am() -> f32 {
    1.0f32 / (IM as f32)
}
/// `RNMX = 1.-EPS` with `EPS = 1.2e-7`.
fn rnmx() -> f32 {
    1.0f32 - 1.2e-7f32
}

#[derive(Clone, Debug, Default)]
pub struct Ran1 {
    /// Number of calls made (written to the `.ato` file; informative only).
    pub mcall: i32,
    pub idum: i32,
    pub iy: i32,
    pub iv: [i32; NTAB],
}

impl Ran1 {
    /// Uniform deviate in (0,1) — `FUNCTION ran1()`.
    pub fn ran1(&mut self) -> f32 {
        if self.idum <= 0 || self.iy == 0 {
            self.idum = (-self.idum).max(1);
            let mut j = NTAB + 8;
            while j >= 1 {
                let k = self.idum / IQ;
                self.idum = IA.wrapping_mul(self.idum - k * IQ).wrapping_sub(IR.wrapping_mul(k));
                if self.idum < 0 {
                    self.idum = self.idum.wrapping_add(IM);
                }
                if j <= NTAB {
                    self.iv[j - 1] = self.idum;
                }
                j -= 1;
            }
            self.iy = self.iv[0];
        }
        let k = self.idum / IQ;
        self.idum = IA.wrapping_mul(self.idum - k * IQ).wrapping_sub(IR.wrapping_mul(k));
        if self.idum < 0 {
            self.idum = self.idum.wrapping_add(IM);
        }
        let j = 1 + self.iy / NDIV;
        self.iy = self.iv[(j - 1) as usize];
        self.iv[(j - 1) as usize] = self.idum;
        let r = (am() * self.iy as f32).min(rnmx());
        self.mcall = self.mcall.wrapping_add(1);
        r
    }

    /// `initran1`: restore state from a list-directed record
    /// `mcall, idum, iy, (iv(i), i=1,ntab)`; on error the generator is reset.
    pub fn from_record(vals: Option<(i32, i32, i32, [i32; NTAB])>) -> Ran1 {
        match vals {
            Some((mcall, idum, iy, iv)) => Ran1 { mcall, idum, iy, iv },
            None => Ran1 { mcall: 0, idum: 0, iy: 0, iv: [0; NTAB] },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn deterministic_sequence() {
        let mut r = Ran1 { idum: -1, ..Default::default() };
        let a = r.ran1();
        let b = r.ran1();
        assert!(a > 0.0 && a < 1.0 && b > 0.0 && b < 1.0 && a != b);
        let mut r2 = Ran1 { idum: -1, ..Default::default() };
        assert_eq!(r2.ran1(), a);
    }
}
