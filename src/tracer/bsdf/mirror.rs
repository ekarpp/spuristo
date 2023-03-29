use super::*;

/// Perfect mirror
pub struct Mirror_BSDF {
    /// Point of impact
    xo: DVec3,
    /// Perfect reflection direction
    wi: DVec3,
}

impl Mirror_BSDF {
    pub fn new(ro: &Ray, ho: &Hit) -> Self {
        let v = -ro.dir;
        let no = ho.norm;
        let wi = 2.0 * v.project_onto(no) - v;
        let xo = ho.p;

        Self {
            xo,
            wi,
        }
    }
}

impl BSDF for Mirror_BSDF {
    fn eval(&self, _ri: &Ray) -> DVec3 { DVec3::ONE }

    fn sample_ray(&self, _rand_sq: DVec2) -> Ray {
        Ray::new(self.xo, self.wi)
    }

    fn prob_for(&self, ri: &Ray) -> f64 {
        let wi = ri.dir;
        if self.wi.dot(wi) >= 1.0 - EPSILON { 1.0 } else { 0.0 }
    }
}
