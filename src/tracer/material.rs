use glam::DVec3;
use crate::tracer::bsdf::{BSDF, Mirror_BSDF, Microfacet_BSDF};
use crate::tracer::hit::Hit;
use crate::tracer::ray::Ray;
use crate::tracer::bxdfs;
use crate::tracer::texture::Texture;
use crate::tracer::microfacet::MfDistribution;

/// Describes which material an object is made out of
pub enum Material {
    /// Glossy
    Microfacet(Texture, MfDistribution),
    /// Emits light
    Light(Texture),
    /// Perfect mirror
    Mirror,
    /// Isotropic medium
    Isotropic(Texture),
    /// Not specified. Used with objects that are built on top of other objects.
    Blank,
}

impl Material {

    /// Metallic microfacet material
    pub fn metal(texture: Texture, roughness: f64) -> Self {
        Self::Microfacet(texture, MfDistribution::metallic(roughness))
    }

    /// Specular microfacet material
    pub fn specular(texture: Texture, roughness: f64) -> Self {
        Self::Microfacet(texture, MfDistribution::specular(roughness))
    }

    /// Diffuse material
    pub fn diffuse(texture: Texture) -> Self {
        Self::Microfacet(texture, MfDistribution::diffuse())
    }

    /// Transparent material
    pub fn transparent(texture: Texture, rfrct_idx: f64, roughness: f64) -> Self {
        Self::Microfacet(
            texture,
            MfDistribution::transparent(rfrct_idx, roughness)
        )
    }

    /// How specular is the material? 1.0 fully, 0.0 none
    pub fn specularity(&self) -> f64 {
        match self {
            Self::Mirror => 1.0,
            Self::Microfacet(_, mfd) => mfd.specularity(),
            _ => 0.0,
        }
    }

    /// duno if this is correct. need to check, hack for now
    pub fn is_transparent(&self) -> bool {
        match self {
            Self::Mirror => true,
            Self::Microfacet(_, mfd) => mfd.is_transparent(),
            _ => false,
        }
    }

    /// Is the material diffuse? I.e. do shadow rays have effect on it.
    /// (CORRECT TERM?)
    pub fn is_diffuse(&self) -> bool {
        matches!(self, Self::Microfacet(..))
    }

    /// How much light emitted at `h`?
    pub fn emit(&self, h: &Hit) -> DVec3 {
        match self {
            Self::Light(t) => t.albedo_at(h.p),
            _ => DVec3::ZERO,
        }
    }
/*
    /// What is the color at `ri.origin`?
    /// Returns BSDF for `ro` at `ho`
    pub fn bsdf(&self, ro: &Ray, ho: &Hit) -> DVec3 {
        let xo = ho.p;
        match self {
            Self::Isotropic(t) => t.albedo_at(xo),
            Self::Mirror => DVec3::ONE,
            Self::Microfacet(t, mfd) => {
                bxdfs::bsdf_microfacet(ro, ri, no, t.albedo_at(xo), mfd)
            }
            _ => DVec3::ZERO,
        }
    }
*/
    /// Returns the BSDF for `ro` at `ho`
    pub fn bsdf(&self, ro: &Ray, ho: &Hit) -> Option<Box<dyn BSDF>> {
        match self {
            Self::Mirror => Some( Box::new(Mirror_BSDF::new(ro, ho)) ),
            Self::Microfacet(tex, mfd) => {
                let no = ho.norm;
                let xo = ho.p;
                let albedo = tex.albedo_at(xo);
                Some( Box::new(Microfacet_BSDF::new(ro, ho, *mfd, albedo)) )
            }
            _ => None,
        }
    }
}
