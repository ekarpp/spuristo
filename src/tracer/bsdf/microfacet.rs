use super::*;
use crate::tracer::microfacet::MfDistribution;

/// Threshold for roughness of microfacet at which we switch to delta pdf
const DELTA_THRESHOLD: f64 = 0.001;

/// Microfacet BSDF
pub struct Microfacet_BSDF {
    /// Distribution of the microfacet
    mfd: MfDistribution,
    /// Albedo at the point
    albedo: DVec3,
    /// Basis for tangent space
    onb: Onb,
    /// Point of impact
    xo: DVec3,
    /// Direction towards viewer from `xo`, in tangent space
    v: DVec3,
}

impl Microfacet_BSDF {
    pub fn new(ro: &Ray, ho: &Hit, mfd: MfDistribution, albedo: DVec3) -> Self {
        let no = ho.norm;
        let onb = Onb::new(no);
        let v = -ro.dir;
        let v = onb.to_local(v);
        let xo = ro.origin;

        Self {
            mfd,
            albedo,
            onb,
            xo,
            v,
        }
    }

    /// Refract direction with Snell-Descartes law.
    ///
    /// # Arguments
    /// * `eta_ratio` - Ratio of refraction indices. `from / to`
    /// * `wh` - Normal for reflraction
    fn refract(&self, eta_ratio: f64, wh: DVec3) -> DVec3 {
        /* Snell-Descartes law */
        let cos_to = wh.dot(self.v);
        let sin2_to = 1.0 - cos_to * cos_to;
        let sin2_ti = eta_ratio * eta_ratio * sin2_to;

        if sin2_ti >= 1.0 {
            /* total internal reflection */
            2.0 * self.v.project_onto(wh) - self.v
        } else {
            let cos_ti = (1.0 - sin2_ti).sqrt();
            -self.v * eta_ratio + (eta_ratio * cos_to - cos_ti) * wh
        }
    }
}

impl BSDF for Microfacet_BSDF {
    fn eval(&self, ri: &Ray) -> DVec3 {
        let wi = self.onb.to_local(ri.dir);

        let no_dot_wi = wi.z;
        let no_dot_v = self.v.z;

        let ro_inside = no_dot_v < 0.0;
        let ri_inside = no_dot_wi < 0.0;
        if ro_inside == ri_inside {
            let wh = (wi + self.v).normalize();
            let no_dot_wh = wh.z;

            let d = self.mfd.d(wh);
            let f = self.mfd.f(self.v, wh, self.albedo);
            let g = self.mfd.g(self.v, wi);

            let specular = d * f * g / (4.0 * no_dot_v * no_dot_wi);

            // BRDF: specular + diffuse, where
            // specular = D(wh) * F(v, wh) * G(v, wi) / (4.0 * (no • v) * (no • wi))
            // diffuse = normalized_disney_term * albedo / π
            // normalized_disney_term = (1.0 + α^2 * (1.0 / 1.51 - 1.0))
            // * (1.0 + (F_90 - 1.0) * (1.0 - (no • v))^5)
            // * (1.0 + (F_90 - 1.0) * (1.0 - (no • wi))^5)
            // F_90 = 0.5 * α^2 + 2.0 * (no • wh)^2 * α^2

            // transparent materials don't have a diffuse term
            if self.mfd.is_transparent() {
                specular
            } else {
                let diffuse = (DVec3::ONE - f) * self.albedo
                    * self.mfd.disney_diffuse(no_dot_v, no_dot_wh, no_dot_wi)
                    / PI;
                diffuse + specular
            }
        } else {
            let eta_ratio = if ro_inside {
                1.0 / self.mfd.get_rfrct_idx()
            } else {
                self.mfd.get_rfrct_idx()
            };

            let wh = (wi * eta_ratio + self.v).normalize();
            let wh = if wh.dot(self.v) < 0.0 { -wh } else { wh };

            let wh_dot_wi = wh.dot(wi);
            let wh_dot_v = wh.dot(self.v);

            let d = self.mfd.d(wh);
            let f = self.mfd.f(self.v, wh, self.albedo);
            let g = self.mfd.g(self.v, wi);

            // BTDF:
            // albedo * abs[(wh • wi) * (wh • v)/((no • wi) * (no • v))]
            // * D(wh) * (1 - F(v, wh)) * G(v, wi) /  (η_r * (wh • wi) + (wh • v))^2

            (wh_dot_wi * wh_dot_v / (no_dot_wi * no_dot_v)).abs()
                * self.albedo * d * (DVec3::ONE - f) * g
                / (eta_ratio * wh_dot_wi + wh_dot_v).powi(2)
        }

    }

    fn sample_ray(&self, rand_sq: DVec2) -> Ray {
        let prob_ndf = self.mfd.probability_ndf_sample(self.albedo);

        let wi = if rand_utils::rand_f64() < prob_ndf {
            let wh = self.mfd.sample_normal(self.v, rand_sq).normalize();
            // reflect v around wh
            2.0 * self.v.project_onto(wh) - self.v
        } else if !self.mfd.is_transparent() {
            rand_utils::square_to_cos_hemisphere(rand_sq)
        } else {
            let inside = self.v.z < 0.0;
            let eta_ratio = if inside {
                self.mfd.get_rfrct_idx()
            } else {
                1.0 / self.mfd.get_rfrct_idx()
            };
            let wh = if self.mfd.get_roughness() <= DELTA_THRESHOLD {
                DVec3::Z
            } else {
                self.mfd.sample_normal(self.v, rand_sq).normalize()
            };

            self.refract(eta_ratio, wh)
        };
        let wi = self.onb.to_world(wi);

        Ray::new(self.xo, wi)
    }

    fn prob_for(&self, ri: &Ray) -> f64 {
        let wi = self.onb.to_local(ri.dir);
        let wh = (self.v + wi).normalize();
        let wh_dot_no = wh.z;
        let wh_dot_v = self.v.dot(wh);
        let v_dot_no = self.v.z;
        // probability to sample wh w.r.t. to wo. mirror??
        //let ndf = self.mfd.d(wh) * wh_dot_no.abs() / (4.0 * wh_dot_v);
        let ndf = self.mfd.g1(self.v) * wh_dot_v.max(0.0)
            * self.mfd.d(wh) / (4.0 * wh_dot_v * v_dot_no);

        // transmission / scatter probability
        let st = if !self.mfd.is_transparent() {
            let cos_theta = wi.z;
            if cos_theta > 0.0 {
                cos_theta * PI.recip()
            } else {
                0.0
            }
        } else if self.v.dot(wi) > 0.0 {
            // in the same hemisphere, zero probability for transmission
            0.0
        } else if self.mfd.get_roughness() <= DELTA_THRESHOLD {
            if 1.0 - self.v.dot(wi) < EPSILON {
                1.0
            } else {
                0.0
            }
        } else {
            let inside = self.v.z < 0.0;
            let eta_ratio = if inside {
                1.0 / self.mfd.get_rfrct_idx()
            } else {
                self.mfd.get_rfrct_idx()
            };
            let wh = (self.v + wi * eta_ratio).normalize();
            let wh_dot_wi = wi.dot(wh);
            let wh_dot_v = wh.dot(self.v);

            self.mfd.d(wh) * wh_dot_no.abs() * (eta_ratio * eta_ratio * wh_dot_wi).abs()
                / (wh_dot_v + eta_ratio * wh_dot_wi).powi(2)
        };

        let prob_ndf = ndf * ndf / (ndf * ndf + st * st);

        prob_ndf * ndf + (1.0 - prob_ndf) * st

    }
}
