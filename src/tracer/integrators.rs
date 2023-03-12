use crate::DVec3;
use crate::rand_utils;
use crate::consts::PATH_TRACE_MAX_DEPTH;
use crate::pdfs::{Pdf, ObjectPdf, CosPdf};
use crate::tracer::hit::Hit;
use crate::tracer::ray::Ray;
use crate::tracer::scene::Scene;

pub trait Integrator {
    fn integrate(&self, r: &Ray, depth: usize) -> DVec3;
}

pub struct PathTracingIntegrator {
    scene: Scene,
}

impl Integrator for PathTracingIntegrator {
    fn integrate(&self, r: &Ray, depth: usize) -> DVec3 {
        if depth > PATH_TRACE_MAX_DEPTH {
            return DVec3::ZERO;
        }

        match self.scene.hit(r) {
            None => DVec3::new(0.0, 1.0, 0.0),
            Some(h) => {
                let material = h.object.material();

                match material.bsdf(&h, r) {
                    Some(sr) => self.integrate(&sr, depth + 1),
                    None => {
                        let pdf_scatter = CosPdf::new(h.norm);
                        let r = Ray::new(
                            h.p,
                            pdf_scatter.generate_dir(
                                rand_utils::rand_unit_square()
                            ),
                        );

                        material.emit(&h)
                            + material.albedo_at(h.p)
                            * self.integrate(&r, depth + 1)
                            * h.norm.dot(r.dir.normalize())
                    }
                }
            }
        }
    }
}

impl PathTracingIntegrator {
    pub fn new(s: Scene) -> Self {
        Self {
            scene: s,
        }
    }
}

pub struct DirectLightingIntegrator {
    scene: Scene,
}

impl Integrator for DirectLightingIntegrator {
    fn integrate(&self, r: &Ray, depth: usize) -> DVec3 {
        /* put this here incase mirrors reflect from each other until infinity */
        if depth > PATH_TRACE_MAX_DEPTH {
            return DVec3::ZERO;
        }

        match self.scene.hit(r) {
            None => DVec3::new(0.0, 1.0, 0.0),
            Some(h) => {
                let material = h.object.material();
                match material.bsdf(&h, r) {
                    Some(sr) => self.integrate(&sr, depth + 1),
                    None => {
                        material.emit(&h)
                            + material.albedo_at(h.p)
                            * self.light_at(&h)

                    }
                }
            }
        }
    }
}

impl DirectLightingIntegrator {
    pub fn new(s: Scene) -> Self {
        Self {
            scene: s,
        }
    }

    fn light_at(&self, h: &Hit) -> f64 {
        let light = self.scene.uniform_random_light();

        let pdf_light = ObjectPdf::new(
            light,
            h.p,
        );
        let pdf_scatter = CosPdf::new(
            h.norm,
        );
        let r = Ray::new(
            h.p,
            pdf_light.generate_dir(rand_utils::rand_unit_square()),
        );
        match self.scene.hit_light(&r, &light) {
            None => 0.0,
            Some(_lh) => {
                5.0 * pdf_scatter.pdf_val(r.dir)
                    / pdf_light.pdf_val(r.dir)
            }
        }
    }
}
