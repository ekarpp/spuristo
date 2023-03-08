use crate::{DVec3, DMat3};
use crate::perlin::Perlin;
use crate::rand_utils;

use crate::tracer::object::{Object, Plane, Rectangle, Cuboid};
use crate::tracer::object::sphere::Sphere;
use crate::tracer::hit::Hit;
use crate::tracer::ray::Ray;
use crate::tracer::material::Material;
use crate::tracer::texture::Texture;

#[cfg(test)]
mod scene_tests;

pub struct Scene {
    pub ambient: DVec3,
    light: Sphere,
    objects: Vec<Box<dyn Object>>,
}

const LIGHT_R: f64 = 0.1;

impl Scene {
    pub fn new(l: DVec3, amb: DVec3, objs: Vec<Box<dyn Object>>) -> Self {
        Self {
            ambient: amb,
            /* use Material::Light later on */
            light: *Sphere::new(l, LIGHT_R, Material::Blank),
            objects: objs,
        }
    }

    pub fn size(&self) -> usize { self.objects.len() }

    pub fn hit(&self, r: &Ray) -> Option<Hit> {
        self.objects.iter().map(|obj| obj.hit(r))
            .fold(None, |closest, hit| {
                if closest.is_none() || (hit.is_some() && hit < closest) {
                    hit
                } else {
                    closest
                }
            })
    }

    pub fn to_light(&self, p: DVec3) -> DVec3 {
        self.light.origin - p + LIGHT_R * rand_utils::rand_unit_sphere()
    }

    pub fn in_light(&self, p: DVec3) -> bool {
        let ray_to_light = Ray::new(
            p,
            /* for now, just choose a random point in the light sphere */
            self.to_light(p),
            0,
        );

        let no_block_light = |obj: &&Box<dyn Object>| -> bool {
            obj.hit(&ray_to_light).filter(|hit| {
                !hit.object.is_translucent()
                    /* check if object is behind light */
                    && (hit.p - ray_to_light.origin).length_squared() <
                    (self.light.origin - ray_to_light.origin).length_squared()
            }).is_none()
        };

        self.objects.iter().take_while(no_block_light).count()
            == self.objects.len()
    }

    pub fn hit_light(&self, r: &Ray) -> bool {
        let no_block_light = |obj: &&Box<dyn Object>| -> bool {
            obj.hit(r).filter(|hit| {
                !hit.object.is_translucent()
                    /* check if object is behind light */
                    && (hit.p - r.origin).length_squared() <
                    (self.light.origin - r.origin).length_squared()
            }).is_none()
        };

        self.objects.iter().take_while(no_block_light).count()
            == self.objects.len()
    }

    pub fn default() -> Self {
        let l = DVec3::new(-0.3, 0.2, -0.1);
        Self::new(
            l,
            DVec3::splat(0.15),
            vec![
                // floor
                Plane::new(
                    DVec3::new(0.0, -0.5, 0.0),
                    DVec3::new(0.0, 1.0, 0.0),
                    Material::Phong(Texture::Checkerboard(
                        Box::new(Texture::Checkerboard(
                            Box::new(Texture::Solid(DVec3::ZERO)),
                            Box::new(Texture::Solid(DVec3::ONE)),
                            4.0,
                        )),
                        /* share same perlin between all textures?
                         * could make cool checkers that way */
                        Box::new(Texture::Marble(
                            Perlin::new(DVec3::splat(192.0) / 255.9)
                        )),
                        1.0,
                    )),
                ),
                // right
                Plane::new(
                    DVec3::new(3.0, 0.0, -3.0),
                    DVec3::new(-1.0, 0.0, 1.0),
                    Material::Phong(Texture::Solid(DVec3::new(0.0, 0.0, 1.0))),
                ),
                Rectangle::new(
                    DMat3::from_cols(
                        DVec3::new(1.2, 0.2, -0.8),
                        DVec3::new(0.8, 0.6, -0.4),
                        DVec3::new(0.4, 0.6, -0.8),
                    ),
                    Material::Mirror,
                ),
                Cuboid::new(
                    DMat3::from_cols(
                        DVec3::new(-0.6, -0.5, -0.7),
                        DVec3::new(-0.5, -0.5, -0.7),
                        DVec3::new(-0.5, -0.5, -0.8),
                    ),
                    DMat3::from_cols(
                        DVec3::new(-0.6, -0.4, -0.7),
                        DVec3::new(-0.5, -0.4, -0.7),
                        DVec3::new(-0.5, -0.4, -0.8),
                    ),
                    Material::Phong(Texture::Checkerboard(
                        Box::new(Texture::Solid(DVec3::new(1.0, 0.0, 1.0))),
                        Box::new(Texture::Solid(
                            DVec3::new(50.0, 205.0, 50.0) / 255.9
                        )),
                        9.0,
                    )),
                ),
                // left
                Plane::new(
                    DVec3::new(-3.0, 0.0, -3.0),
                    DVec3::new(1.0, 0.0, 1.0),
                    Material::Phong(Texture::Solid(DVec3::new(1.0, 0.0, 0.0))),
                ),
                // behind
                Plane::new(
                    DVec3::new(0.0, 0.0, 1.0),
                    DVec3::new(0.0, 0.0, -1.0),
                    Material::Phong(Texture::Solid(DVec3::new(1.0, 0.0, 1.0))),
                ),
                Sphere::new(
                    DVec3::new(0.0, 0.0, -1.0),
                    0.5,
                    Material::Phong(Texture::Solid(
                        DVec3::new(136.0, 8.0, 8.0) / 255.9
                    )),
                ),
                Sphere::new(
                    DVec3::new(-0.9, 0.0, -1.0),
                    0.1,
                    Material::Mirror,
                ),
                Sphere::new(
                    DVec3::new(-0.4, -0.12, -0.5),
                    0.1,
                    Material::Glass,
                ),
                Sphere::new(
                    DVec3::new(0.4, -0.2, -0.5),
                    0.1,
                    Material::Phong(Texture::Marble(Perlin::new(
                        DVec3::new(255.0, 182.0, 193.0) / 255.9
                    ))),
                ),
            ]
        )
    }
}
