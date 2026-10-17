# Changes made to the Rust engine for the Python wrapper

None. The wrapper uses only the public API of `rust/epsr-core` as it stands; no file
under `rust/` was modified.

If a change ever becomes necessary, it must be additive (a new accessor appended to
`rust/epsr-core/src/api.rs`), must not rename or move anything, and must be listed here.
