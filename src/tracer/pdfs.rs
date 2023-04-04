use crate::rand_utils;
use crate::tracer::bxdfs;
use crate::tracer::microfacet::MfDistribution;
use crate::tracer::object::Object;
use crate::tracer::onb::Onb;
use crate::tracer::ray::Ray;
use crate::EPSILON;
use glam::{DVec2, DVec3};
use std::f64::consts::PI;

/// Assumes that each generation and evaluation has same starting point. DO AS ENUM?
pub trait Pdf {
    /// Generates a random direction according to the sampling strategy
    ///
    /// # Arguments
    /// * `rand_sq` - Random point on the unit square.
    fn sample_ray(&self, rand_sq: DVec2) -> Option<Ray>;
    /// Computes the probability of the given direction.
    ///
    /// # Arguments
    /// * `ri` - Ray to compute probability for
    fn value_for(&self, ri: &Ray) -> f64;
}

/// TODO
pub struct IsotropicPdf {
    xo: DVec3,
}

impl IsotropicPdf {
    pub fn new(xo: DVec3) -> Self {
        Self { xo }
    }
}

impl Pdf for IsotropicPdf {
    fn sample_ray(&self, rand_sq: DVec2) -> Option<Ray> {
        let wi = rand_utils::square_to_sphere(rand_sq);
        Some( Ray::new(self.xo, wi) )
    }

    fn value_for(&self, _ri: &Ray) -> f64 {
        let d: f64 = 0.1; //hi.object.density();
                          // hi.t = 1.5
        d * (-d * 1.5).exp()
    }
}

/// Randomly samples a direction towards a point on the object that is visible
pub struct ObjectPdf<'a> {
    /// Object to do sampling from
    object: &'a dyn Object,
    /// Point from where the object should be visible
    xo: DVec3,
}

impl<'a> ObjectPdf<'a> {
    pub fn new(object: &'a dyn Object, xo: DVec3) -> Self {
        Self { object, xo }
    }
}

impl Pdf for ObjectPdf<'_> {
    fn sample_ray(&self, rand_sq: DVec2) -> Option<Ray> {
        Some( self.object.sample_towards(self.xo, rand_sq) )
    }

    fn value_for(&self, ri: &Ray) -> f64 {
        self.object.sample_towards_pdf(ri)
    }
}

/// Delta distribution PDF. Always samples the same ray. For glass/mirror.
pub struct DeltaPdf {
    xo: DVec3,
    wi: DVec3,
}

impl DeltaPdf {
    pub fn new(xo: DVec3, wi: DVec3) -> Self {
        Self { xo, wi }
    }
}

impl Pdf for DeltaPdf {
    fn sample_ray(&self, _rand_sq: DVec2) -> Option<Ray> {
        Some( Ray::new(self.xo, self.wi) )
    }

    fn value_for(&self, ri: &Ray) -> f64 {
        let wi = ri.dir;
        if wi.dot(self.wi) >= 1.0 - EPSILON {
            1.0
        } else {
            0.0
        }
    }
}

/// PDF for microfacet distribution.
pub struct MfdPdf {
    /// Point of impact
    xo: DVec3,
    /// Direction from point of impact to viewer
    v: DVec3,
    /// Macrosurface geometric normal. Same hemisphere as `v`.
    ng: DVec3,
    /// Probability to sample ray from NDF
    ndf_sample_prob: f64,
    /// ONB for macrosurface normal
    uvw: Onb,
    /// The microfacet distribution of the surface
    mfd: MfDistribution,
}

// add cos pdf. add reflection pdf. reflection pdf has delta pdf?
impl MfdPdf {
    pub fn new(xo: DVec3, v: DVec3, ng: DVec3, albedo: DVec3, mfd: MfDistribution) -> Self {
        // refraction needs v and wh to be in same hemisphere so we do this
        let w = if v.dot(ng) < 0.0 { -ng } else { ng };
        let uvw = Onb::new(w);

        Self {
            xo,
            v,
            uvw,
            ndf_sample_prob: mfd.probability_ndf_sample(albedo),
            ng,
            mfd,
        }
    }
}

// NEEDS A REFACTOR. split to 3? different pdfs? (cos, refract, ndf?)
impl Pdf for MfdPdf {
    /// Sample microsurface normal from the distribution. Mirror direction from
    /// camera around the normal. GGX uses VNDF sampling, Beckmann NDF sampling
    fn sample_ray(&self, rand_sq: DVec2) -> Option<Ray> {
        let wi = if rand_utils::rand_f64() < self.ndf_sample_prob {
            let local_v = self.uvw.to_local(self.v);
            let local_wh = self.mfd.sample_normal(local_v, rand_sq).normalize();
            let local_wi = bxdfs::reflect(local_v, local_wh);

            if local_wi.z < 0.0 {
                // bad sample, do something else?
                return None;
            }

            self.uvw.to_world(local_wi)
        } else if !self.mfd.is_transparent() {
            self.uvw
                .to_world(rand_utils::square_to_cos_hemisphere(rand_sq))
        } else {
            let inside = self.ng.dot(self.v) < 0.0;
            let eta_ratio = if inside {
                self.mfd.get_rfrct_idx()
            } else {
                1.0 / self.mfd.get_rfrct_idx()
            };
            let local_v = self.uvw.to_local(self.v);
            let local_wh = self.mfd.sample_normal(local_v, rand_sq).normalize();

            let wh = self.uvw.to_world(local_wh).normalize();

            match bxdfs::refract(eta_ratio, self.v, wh) {
                None => bxdfs::reflect(self.v, wh),
                Some(wi) => wi,
            }
        };

        Some( Ray::new(self.xo, wi) )
    }

    /// Read it directly from the NFD and do change of variables
    /// from `wh` to `wi`.
    fn value_for(&self, ri: &Ray) -> f64 {
        let wi = ri.dir;
        let wh = (self.v + wi).normalize();

        let wh_dot_v = self.v.dot(wh);

        // probability to sample wh w.r.t. to v
        let ndf = self.mfd.sample_normal_pdf(wh, self.v, self.ng)
            // jacobian
            / (4.0 * wh_dot_v);

        // transmission / scatter probability
        let st = if !self.mfd.is_transparent() {
            let cos_theta = self.uvw.w.dot(wi);
            if cos_theta > 0.0 {
                cos_theta / PI
            } else {
                0.0
            }
        } else {
            // refraction

            // check if same hemisphere w.r.t. ng?

            let inside = self.ng.dot(self.v) < 0.0;
            let eta_ratio = if inside {
                1.0 / self.mfd.get_rfrct_idx()
            } else {
                self.mfd.get_rfrct_idx()
            };
            let wh = -(self.v + wi * eta_ratio).normalize();
            let wh_dot_wi = wi.dot(wh);
            let wh_dot_v = wh.dot(self.v);

            if wh_dot_wi * wh_dot_v > 0.0 {
                // same hemisphere w.r.t wh, zero probability for refraction
                0.0
            } else {
                // wh and ng need to be in same hemisphere, hemisphere of v makes
                // no difference.
                self.mfd.sample_normal_pdf(wh, self.v, self.ng)
                    // jacobian
                    * (eta_ratio * eta_ratio * wh_dot_wi).abs()
                    / (wh_dot_v + eta_ratio * wh_dot_wi).powi(2)
            }
        };

        self.ndf_sample_prob * ndf + (1.0 - self.ndf_sample_prob) * st
    }
}
