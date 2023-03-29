use glam::DVec3;

/// Small utility struct for orthonormal basis
pub struct Onb {
    /// `x` (or `y`?) direction
    pub u: DVec3,
    /// `y` (or `x`?) direction
    pub v: DVec3,
    /// `z` direction
    pub w: DVec3,
}

impl Onb {
    /// Creates a new ONB.
    ///
    /// # Arguments
    /// * `dir` - Direction of `z` axis. Not necessarily normalized.
    pub fn new(dir: DVec3) -> Self {
        let w = dir.normalize();
        let (u, v) = w.any_orthonormal_pair();
        Self { u, v, w }
    }

    /// Translate from the ONB to world basis
    ///
    /// # Arguments
    /// * `v` - The vector in ONB basis.
    pub fn to_world(&self, v: DVec3) -> DVec3 {
        v.x * self.u + v.y * self.v + v.z * self.w
    }

    /// Translate from world basis to the ONB
    ///
    /// # Arguments
    /// * `v` - The vector in canonical basis
    pub fn to_local(&self, v: DVec3) -> DVec3 {
        DVec3::new(
            v.dot(self.u),
            v.dot(self.v),
            v.dot(self.w),
        )
    }
}
