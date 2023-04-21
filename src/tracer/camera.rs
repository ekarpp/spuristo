use crate::tracer::ray::Ray;
use crate::tracer::onb::Onb;
use glam::{DVec2, DVec3, IVec2};

pub use perspective_camera::PerspectiveCamera;
pub use orthographic_camera::OrthographicCamera;

/// Perspective
mod perspective_camera;
/// Orthographic camera that preserves angles
mod orthographic_camera;

fn _camera_basis(origin: DVec3, towards: DVec3, up: DVec3) -> Onb {
    assert!(origin != towards);

    let forward = (towards - origin).normalize();
    let right = forward.cross(up).normalize();
    let down = forward.cross(right);

    // x = right, y = down, z = towards
    Onb::new_from_basis(right, down, forward)
}

/// Camera trait. Might want to make this an enum
pub trait Camera: Sync {
    /// Returns a ray pointing towards a point on the image plane.
    ///
    /// # Arguments
    /// * `xy` - Position in `\[0,width\] x \[0,height\]` raster space
    fn generate_ray(&self, raster_xy: DVec2) -> Ray;

    /// Returns the image resolution
    fn get_resolution(&self) -> IVec2;
}
