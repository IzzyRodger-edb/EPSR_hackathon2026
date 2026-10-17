//! Fortran intrinsic semantics needed for bit-identical arithmetic.
//!
//! gfortran (at -O2, without -ffast-math) evaluates single-precision
//! expressions strictly left to right in IEEE `f32`.  The only places where the
//! generated code is not obvious are integer powers:
//!
//! * `x**2`, `x**3`, `x**4` with a *constant* exponent are expanded inline to
//!   `x*x`, `(x*x)*x`, `(x*x)*(x*x)` — use [`pow2`], [`pow3`], [`pow4`];
//! * `x**n` with a *variable* integer exponent calls libgcc's `__powisf2`
//!   (binary exponentiation) — use [`powi_var`];
//! * `x**y` with a real exponent calls `powf` — use `f32::powf`.
//!
//! Transcendental functions (`exp`, `sin`, `cos`, `atan`, `log`, `sinh`, `sqrt`)
//! map to the same glibc `libm` routines from both gfortran and Rust's `std`.

/// `x**2` (constant exponent).
#[inline(always)]
pub fn pow2(x: f32) -> f32 {
    x * x
}
/// `x**3` (constant exponent): `(x*x)*x`.
#[inline(always)]
pub fn pow3(x: f32) -> f32 {
    x * x * x
}
/// `x**4` (constant exponent): `(x*x)*(x*x)`.
#[inline(always)]
pub fn pow4(x: f32) -> f32 {
    let x2 = x * x;
    x2 * x2
}
/// `x**2` for `real*8`.
#[inline(always)]
pub fn dpow2(x: f64) -> f64 {
    x * x
}

/// `x**m` with a run-time integer exponent: exact port of libgcc `__powisf2`.
pub fn powi_var(x: f32, m: i32) -> f32 {
    let mut n: u32 = if m < 0 { (m as i64).unsigned_abs() as u32 } else { m as u32 };
    let mut x = x;
    let mut y: f32 = if n % 2 == 1 { x } else { 1.0 };
    loop {
        n >>= 1;
        if n == 0 {
            break;
        }
        x = x * x;
        if n % 2 == 1 {
            y = y * x;
        }
    }
    if m < 0 {
        1.0 / y
    } else {
        y
    }
}

/// Fortran `NINT` for default real: round half away from zero.
///
/// Implemented inline (rather than `f32::round`, which is a libm call on
/// baseline x86-64) because it sits in the innermost loop of every energy
/// evaluation (minimum image).  `x - trunc(x)` is exact in floating point, so
/// the result equals `lroundf(x)` for every finite `x` within `i32` range.
#[inline(always)]
pub fn nint(x: f32) -> i32 {
    let i = x as i32;
    let frac = x - i as f32;
    if frac >= 0.5 {
        i + 1
    } else if frac <= -0.5 {
        i - 1
    } else {
        i
    }
}

/// Fortran `INT` for default real: truncation toward zero.
#[inline(always)]
pub fn int(x: f32) -> i32 {
    x as i32
}

/// Fortran `SIGN(a, b)`.
#[inline(always)]
pub fn sign(a: f32, b: f32) -> f32 {
    a.abs().copysign(b)
}

// ---- gfortran intrinsic SIND / COSD -----------------------------------------
// A Fortran unit that calls `sind`/`cosd` without declaring them EXTERNAL gets
// libgfortran's intrinsics, which differ from EPSR's own `sind`/`cosd`
// (auxiliary_routines.f, see `ato::sind`).  Ported from libgfortran.


const PIO180H: f32 = 0.0174560546875; // 0x3c8f0000
const PIO180L: f32 = -2.7621674689726206e-06; // 0xb6395dad

#[inline]
fn d2r(t: f32) -> f32 {
    t.mul_add(PIO180H, t * PIO180L)
}

/// Integer degree value of `z` if `z` is integral (`(int)z == z`).
fn integral(z: f32) -> Option<i32> {
    let n = z as i32;
    if n as f32 - z == 0.0 {
        Some(n)
    } else {
        None
    }
}

/// gfortran intrinsic `SIND(x)` for `real(4)`.
pub fn sind_intrinsic(x: f32) -> f32 {
    if !x.is_finite() {
        return x - x;
    }
    let ax = x.abs();
    if ax < 0.03125 {
        return x.mul_add(PIO180H, x * PIO180L);
    }
    let mut s = 1.0f32.copysign(x);
    let z = ax % 360.0;
    if let Some(n) = integral(z) {
        if n % 30 == 0 {
            if n % 180 == 0 {
                if n == 180 {
                    s = -s;
                }
                return s * 0.0;
            }
            if n % 90 == 0 {
                return if n == 90 { s } else { -s };
            }
            if n % 60 == 0 {
                return if n > 179 { s * -0.8660253882408142 } else { s * 0.8660253882408142 };
            }
            return if n > 179 { s * -0.5 } else { s * 0.5 };
        }
    }
    if z <= 180.0 {
        if z <= 90.0 {
            if z <= 45.0 {
                s * d2r(z).sin()
            } else {
                s * d2r(90.0 - z).cos()
            }
        } else if z <= 135.0 {
            s * d2r(z - 90.0).cos()
        } else {
            s * d2r(180.0 - z).sin()
        }
    } else {
        s = -s;
        if z <= 270.0 {
            if z <= 225.0 {
                s * d2r(z - 180.0).sin()
            } else {
                s * d2r(270.0 - z).cos()
            }
        } else if z <= 315.0 {
            s * d2r(z - 270.0).cos()
        } else {
            s * d2r(360.0 - z).sin()
        }
    }
}

/// gfortran intrinsic `COSD(x)` for `real(4)`.
pub fn cosd_intrinsic(x: f32) -> f32 {
    if !x.is_finite() {
        return x - x;
    }
    let ax = x.abs();
    if ax <= 0.0078125 {
        return if x == 0.0 { 1.0 } else { 1.0 - 7.888609052210118e-31 };
    }
    let z = ax % 360.0;
    if let Some(n) = integral(z) {
        if n % 30 == 0 {
            if n % 180 == 0 {
                return if n == 180 { -1.0 } else { 1.0 };
            }
            if n % 90 == 0 {
                return 0.0;
            }
            if n % 60 == 0 {
                return if n == 60 || n == 300 { 0.5 } else { -0.5 };
            }
            return if n == 30 || n == 330 { 0.8660253882408142 } else { -0.8660253882408142 };
        }
    }
    if z <= 180.0 {
        if z <= 90.0 {
            if z <= 45.0 {
                d2r(z).cos()
            } else {
                d2r(90.0 - z).sin()
            }
        } else if z <= 135.0 {
            -d2r(z - 90.0).sin()
        } else {
            -d2r(180.0 - z).cos()
        }
    } else if z <= 270.0 {
        if z <= 225.0 {
            -d2r(z - 180.0).cos()
        } else {
            -d2r(270.0 - z).sin()
        }
    } else if z <= 315.0 {
        d2r(z - 270.0).sin()
    } else {
        d2r(360.0 - z).cos()
    }
}


/// The value gfortran folds `4.0*atan(1.0)` to.
pub const PI4ATAN: f32 = std::f32::consts::PI;

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn powi_matches_repeated_squaring() {
        assert_eq!(powi_var(2.0, 10), 1024.0);
        assert_eq!(powi_var(1.5, 0), 1.0);
        assert_eq!(powi_var(2.0, -2), 0.25);
        assert_eq!(nint(2.5), 3);
        assert_eq!(nint(-2.5), -3);
        assert_eq!(int(-2.7), -2);
        assert_eq!(PI4ATAN, 4.0 * 1.0f32.atan());
    }
}

#[cfg(test)]
mod nint_tests {
    use super::nint;
    #[test]
    fn nint_matches_round_half_away() {
        let mut x: f32 = -50.0;
        while x < 50.0 {
            assert_eq!(nint(x), x.round() as i32, "x = {x}");
            x += 0.0137;
        }
        for v in [0.5f32, -0.5, 1.5, -1.5, 2.5, 0.49999997, -0.49999997, 8388607.5, -8388607.5, 1.0e9, -1.0e9] {
            assert_eq!(nint(v), v.round() as i32, "x = {v}");
        }
    }
}
