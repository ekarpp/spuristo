use crate::rand_utils;
use crate::samplers::JitteredSampler;
use crate::tracer::hit::Hit;
use crate::tracer::pdfs::{ObjectPdf, Pdf};
use crate::tracer::bsdf::{BSDF, Mirror_BSDF, Microfacet_BSDF};
use crate::tracer::ray::Ray;
use crate::tracer::scene::Scene;
use glam::{DVec2, DVec3};
use std::fmt;

mod bd_path_trace;
mod direct_light;
mod path_trace;

/// How many shadow rays per vertex in path tracer? Preferably square for
/// jittered sampler.
const SHADOW_SPLITS: u32 = 1;

/// Russian roulette probability for the path tracer.
/// Terminates a path at each step with this probability.
/// Computed values are multiplied by the reciprocal of the inverse probability.
const PATH_TRACE_RR: f64 = 0.2;

/// Enum to choose which integrator to use
pub enum Integrator {
    /// Implements the path tracing algorithm with
    /// Russian Roulette (With probability `p` terminate each path.
    /// Multiply contributions by reciprocal of `1-p`) and
    /// next event estimation (Importance sample light at each impact).
    PathTrace,
    /// Naive integrator that importance samples light once.
    DirectLight,
    /// Bidirectional path tracing. Not implemented.
    BDPathTrace,
}

impl fmt::Display for Integrator {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PathTrace => write!(f, "path tracing"),
            Self::DirectLight => write!(f, "direct light integration"),
            Self::BDPathTrace => write!(f, "bidirectional path tracing"),
        }
    }
}

impl Integrator {
    /// Calls the corresponding integration function
    pub fn integrate(&self, s: &Scene, r: &Ray) -> DVec3 {
        match self {
            Self::PathTrace => path_trace::integrate(s, r, 1.0),
            Self::DirectLight => direct_light::integrate(s, r),
            Self::BDPathTrace => bd_path_trace::integrate(s, r),
        }
    }
}

/// Shoots a shadow ray towards random light from `ho`.
fn shadow_ray(scene: &Scene, ro: &Ray, ho: &Hit, bsdf: &dyn BSDF, rand_sq: DVec2) -> DVec3 {
    let material = ho.object.material();

    if !material.is_diffuse() {
        DVec3::ZERO
    } else {
        let xo = ho.p;
        let no = ho.norm;
        let light = scene.uniform_random_light();

        let pdf_light = ObjectPdf::new(light, xo);
        let ri = pdf_light.sample_ray(rand_sq);
        let wi = ri.dir;

        match scene.hit_light(&ri, light) {
            None => DVec3::ZERO,
            Some(_) => {
                let p_light = pdf_light.value_for(&ri);
                let p_scatter = bsdf.prob_for(&ri);

                let weight = p_light * p_light
                    / (p_light * p_light + p_scatter * p_scatter);

                bsdf.eval(&ri)
                    * no.dot(wi).abs()
                    //* weight
                    / p_light
            }
        }
    }
}
