use glam::{DVec2, DVec3};
use std::f64::consts::PI;
use crate::EPSILON;
use crate::rand_utils;
use crate::tracer::hit::Hit;
use crate::tracer::onb::Onb;
use crate::tracer::ray::Ray;

pub use microfacet::Microfacet_BSDF;
pub use mirror::Mirror_BSDF;

/// Microfacet BSDF
mod microfacet;
/// Perfect mirror BSDF
mod mirror;

/// Trait for different BSDFs
pub trait BSDF {
    /// Evaluate the BSDF at for scattered ray `ri`
    fn eval(&self, ri: &Ray) -> DVec3;
    /// Samples a ray from the BSDF
    fn sample_ray(&self, rand_sq: DVec2) -> Ray;
    /// Probability that `ri` got sampled
    fn prob_for(&self, ri: &Ray) -> f64;
}
