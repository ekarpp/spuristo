pub use camera::Camera;
pub use integrator::Integrator;
pub use material::Material;
pub use object::{Bounded, Disk, Instance, Instanceable, KdTree, Mesh, Object};
pub use object::{Cone, Cube, Cylinder, Plane, Rectangle, Sphere, Triangle};
pub use scene::Scene;
pub use texture::Texture;

/// Different BSDFs.
mod bxdfs;
mod bsdf;
/// Abstraction for a camera
mod camera;
/// Abstraction for a hit between a ray and an object.
mod hit;
/// Integrator to estimate the irradiance at each point
mod integrator;
/// Material of an object that defines how it behaves with rays
mod material;
/// MFDistribution
mod microfacet;
/// Abstractions for objects in the 3D world
mod object;
/// Utility struct for orthonormal basis.
mod onb;
/// Implementation of different probability density functions for sampling.
mod pdfs;
/// Abstractions for rays.
mod ray;
/// Scene that describes the 3D world to render.
mod scene;
/// Textures that can be given to some materials
mod texture;
