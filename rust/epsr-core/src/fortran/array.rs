//! Fortran-style arrays: 1-based indexing, column-major storage.
//!
//! The EPSR Fortran source indexes every array from 1 and stores multi-dimensional
//! arrays in column-major order.  To keep the translation line-by-line comparable
//! with the original (and so make numerical parity easy to audit) the Rust port
//! uses these thin wrappers instead of re-basing every index expression.
//!
//! * [`F1`] — `real x(n)`            → `x[i]`
//! * [`F2`] — `real a(n1,n2)`        → `a[(i, j)]`
//! * [`F3`] — `real a(n1,n2,n3)`     → `a[(i, j, k)]`
//!
//! Indices may be given as `i32` (the natural type of a translated Fortran
//! `integer`) or `usize`.  Out-of-range access panics (equivalent to compiling
//! the Fortran with `-fbounds-check`).
//!
//! The underlying storage is exposed through `as_slice()` so that a column of an
//! `F2` can be handed to routines that, in Fortran, received an array section or
//! the address of `a(1,j)`.

use std::ops::{Index, IndexMut};

/// Convert a 1-based Fortran index to a 0-based offset, panicking on index < 1.
#[inline(always)]
fn off(i: i64, n: usize, what: &str) -> usize {
    debug_assert!(i >= 1 && (i as usize) <= n, "{what}: index {i} out of 1..={n}");
    if i < 1 || i as usize > n {
        panic!("{what}: index {i} out of range 1..={n}");
    }
    (i - 1) as usize
}

/// Types usable as a Fortran index.
pub trait FIdx: Copy {
    fn fi(self) -> i64;
}
impl FIdx for i32 {
    #[inline(always)]
    fn fi(self) -> i64 {
        self as i64
    }
}
impl FIdx for usize {
    #[inline(always)]
    fn fi(self) -> i64 {
        self as i64
    }
}
impl FIdx for i64 {
    #[inline(always)]
    fn fi(self) -> i64 {
        self
    }
}

/// One-dimensional, 1-based array (`dimension(n)`).
#[derive(Clone, Debug, Default)]
pub struct F1<T> {
    pub data: Vec<T>,
}

impl<T: Clone + Default> F1<T> {
    /// Allocate `n` default-initialised elements (Fortran COMMON storage is zeroed).
    pub fn new(n: usize) -> Self {
        F1 { data: vec![T::default(); n] }
    }
    /// Equivalent of `reallocate1d_*`: grow (never shrink) keeping old contents.
    pub fn grow(&mut self, n: usize) {
        if self.data.len() < n {
            self.data.resize(n, T::default());
        }
    }
    /// Fill every element with `v` (Fortran whole-array assignment `a = v`).
    pub fn fill(&mut self, v: T) {
        for x in self.data.iter_mut() {
            *x = v.clone();
        }
    }
}

impl<T> F1<T> {
    pub fn len(&self) -> usize {
        self.data.len()
    }
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
    pub fn as_slice(&self) -> &[T] {
        &self.data
    }
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data
    }
}

impl<T, I: FIdx> Index<I> for F1<T> {
    type Output = T;
    #[inline(always)]
    fn index(&self, i: I) -> &T {
        let o = off(i.fi(), self.data.len(), "F1");
        unsafe { self.data.get_unchecked(o) }
    }
}
impl<T, I: FIdx> IndexMut<I> for F1<T> {
    #[inline(always)]
    fn index_mut(&mut self, i: I) -> &mut T {
        let o = off(i.fi(), self.data.len(), "F1");
        unsafe { self.data.get_unchecked_mut(o) }
    }
}

/// Two-dimensional, 1-based, column-major array (`dimension(n1,n2)`).
#[derive(Clone, Debug, Default)]
pub struct F2<T> {
    pub n1: usize,
    pub n2: usize,
    pub data: Vec<T>,
}

impl<T: Clone + Default> F2<T> {
    pub fn new(n1: usize, n2: usize) -> Self {
        F2 { n1, n2, data: vec![T::default(); n1 * n2] }
    }
    /// Reshape keeping the (i,j) contents that still fit (like `reallocate2d_*`).
    pub fn grow(&mut self, n1: usize, n2: usize) {
        if n1 <= self.n1 && n2 <= self.n2 {
            return;
        }
        let m1 = n1.max(self.n1);
        let m2 = n2.max(self.n2);
        let mut nd = vec![T::default(); m1 * m2];
        for j in 0..self.n2 {
            for i in 0..self.n1 {
                nd[i + j * m1] = self.data[i + j * self.n1].clone();
            }
        }
        self.n1 = m1;
        self.n2 = m2;
        self.data = nd;
    }
    pub fn fill(&mut self, v: T) {
        for x in self.data.iter_mut() {
            *x = v.clone();
        }
    }
}

impl<T> F2<T> {
    /// Contiguous column `j` (Fortran `a(:,j)`), 0-based slice.
    pub fn col(&self, j: usize) -> &[T] {
        let o = off(j as i64, self.n2, "F2 col");
        &self.data[o * self.n1..(o + 1) * self.n1]
    }
    pub fn col_mut(&mut self, j: usize) -> &mut [T] {
        let o = off(j as i64, self.n2, "F2 col");
        let n1 = self.n1;
        &mut self.data[o * n1..(o + 1) * n1]
    }
}

impl<T, I: FIdx, J: FIdx> Index<(I, J)> for F2<T> {
    type Output = T;
    #[inline(always)]
    fn index(&self, (i, j): (I, J)) -> &T {
        let oi = off(i.fi(), self.n1, "F2 dim1");
        let oj = off(j.fi(), self.n2, "F2 dim2");
        unsafe { self.data.get_unchecked(oi + oj * self.n1) }
    }
}
impl<T, I: FIdx, J: FIdx> IndexMut<(I, J)> for F2<T> {
    #[inline(always)]
    fn index_mut(&mut self, (i, j): (I, J)) -> &mut T {
        let oi = off(i.fi(), self.n1, "F2 dim1");
        let oj = off(j.fi(), self.n2, "F2 dim2");
        let n1 = self.n1;
        unsafe { self.data.get_unchecked_mut(oi + oj * n1) }
    }
}

/// Three-dimensional, 1-based, column-major array (`dimension(n1,n2,n3)`).
#[derive(Clone, Debug, Default)]
pub struct F3<T> {
    pub n1: usize,
    pub n2: usize,
    pub n3: usize,
    pub data: Vec<T>,
}

impl<T: Clone + Default> F3<T> {
    pub fn new(n1: usize, n2: usize, n3: usize) -> Self {
        F3 { n1, n2, n3, data: vec![T::default(); n1 * n2 * n3] }
    }
    pub fn fill(&mut self, v: T) {
        for x in self.data.iter_mut() {
            *x = v.clone();
        }
    }
}

impl<T, I: FIdx, J: FIdx, K: FIdx> Index<(I, J, K)> for F3<T> {
    type Output = T;
    #[inline(always)]
    fn index(&self, (i, j, k): (I, J, K)) -> &T {
        let oi = off(i.fi(), self.n1, "F3 dim1");
        let oj = off(j.fi(), self.n2, "F3 dim2");
        let ok = off(k.fi(), self.n3, "F3 dim3");
        unsafe { self.data.get_unchecked(oi + self.n1 * (oj + self.n2 * ok)) }
    }
}
impl<T, I: FIdx, J: FIdx, K: FIdx> IndexMut<(I, J, K)> for F3<T> {
    #[inline(always)]
    fn index_mut(&mut self, (i, j, k): (I, J, K)) -> &mut T {
        let oi = off(i.fi(), self.n1, "F3 dim1");
        let oj = off(j.fi(), self.n2, "F3 dim2");
        let ok = off(k.fi(), self.n3, "F3 dim3");
        let (n1, n2) = (self.n1, self.n2);
        unsafe { self.data.get_unchecked_mut(oi + n1 * (oj + n2 * ok)) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn column_major_layout() {
        let mut a: F2<i32> = F2::new(3, 2);
        a[(1, 1)] = 11;
        a[(3, 1)] = 31;
        a[(1, 2)] = 12;
        assert_eq!(a.data, vec![11, 0, 31, 12, 0, 0]);
        assert_eq!(a.col(2)[0], 12);
    }
    #[test]
    fn grow_keeps_contents() {
        let mut a: F2<i32> = F2::new(2, 2);
        a[(2, 2)] = 5;
        a.grow(4, 3);
        assert_eq!(a[(2, 2)], 5);
        assert_eq!(a[(4, 3)], 0);
    }
    #[test]
    #[should_panic]
    fn zero_index_panics() {
        let a: F1<f32> = F1::new(3);
        let _ = a[0i32];
    }
}
