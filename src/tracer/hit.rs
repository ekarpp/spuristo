use crate::efloat;
use crate::tracer::material::Material;
use crate::tracer::object::Sampleable;
use crate::tracer::ray::Ray;
use glam::{DVec2, DVec3};

/// Stores information about a hit between a ray and an object
pub struct Hit<'a> {
    /// The `t` value of ray at which the hit occurred
    pub t: f64,
    /// Material of the object which got hit
    pub material: &'a Material,
    /// 3D point where object was hit
    pub p: DVec3,
    /// Optional reference to light if we hit one
    pub light: Option<&'a dyn Sampleable>,
    /// Floating point error bounds of the impact point
    pub fp_error: DVec3,
    /// Normal of the surface used for shading calculations
    pub ns: DVec3,
    /// Geometric normal of the surface
    pub ng: DVec3,
    /// Texture coordinates in `\[0,1\]^2`
    pub uv: DVec2,
    /// Are we on the backface?
    pub backface: bool,
}

impl<'a> Hit<'a> {
    /// # Arguments
    ///
    /// * `t` - Value of ray at which hit occurred
    /// * `material` - Material of the object which got hit
    /// * `wo` - Direction towards to point of intersection
    /// * `xi` - Point in world space at which object got hit
    /// * `fp_error` - Floating point error bounds for `xi`
    /// * `ns` - Shading normal of the object at the point of impact
    /// * `ng` - Geometric normal of the object at the point of impact
    /// * `uv` - Texture coordinates in `\[0,1\]^2`
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        t: f64,
        material: &'a Material,
        wo: DVec3,
        xi: DVec3,
        fp_error: DVec3,
        ns: DVec3,
        ng: DVec3,
        uv: DVec2,
    ) -> Option<Self> {
        let backface = wo.dot(ng) > 0.0;
        Some(Self {
            t,
            material,
            backface,
            p: xi,
            light: None,
            fp_error,
            ns,
            ng,
            uv,
        })
    }

    /// Generates a ray at point of impact. Would be better to use accurate
    /// error bounds instead of `EPSILON`.
    pub fn generate_ray(&self, wi: DVec3) -> Ray {
        let scaled_err = self.fp_error.dot(self.ns.abs());

        let offset = if wi.dot(self.ng) >= 0.0 {
            self.ns * scaled_err
        } else {
            -self.ns * scaled_err
        };

        let xi = self.p + offset;

        let move_double = |v: f64, n: f64| {
            if n > 0.0 {
                efloat::next_double(v)
            } else if n < 0.0 {
                efloat::previous_double(v)
            } else {
                v
            }
        };

        let xi = DVec3::new(
            move_double(xi.x, offset.x),
            move_double(xi.y, offset.y),
            move_double(xi.z, offset.z),
        );

        Ray::new(
            xi,
            wi
        )
    }

    /// Did we hit a medium?
    pub fn is_medium(&self) -> bool {
        matches!(self.material, Material::Volumetric(..))
    }

    /// Did we hit a light?
    pub fn is_light(&self) -> bool {
        matches!(self.material, Material::Light(..))
    }

    /// Are we on delta material?
    pub fn is_delta(&self) -> bool {
        self.material.is_delta()
    }
}
