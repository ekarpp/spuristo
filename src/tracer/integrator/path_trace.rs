use super::*;

pub fn integrate(
             scene: &Scene,
             r: &Ray,
             depth: usize,
             last_specular: bool
) -> DVec3 {
    if depth > 2 && rand_utils::rand_f64() < PATH_TRACE_RR {
        return DVec3::ZERO;
    }

    match scene.hit(r) {
        None => DVec3::new(0.0, 0.0, 1.0),
        Some(h) => {
            let material = h.object.material();

            match material.bsdf(&h, r) {
                None => if last_specular {
                    material.emit(&h)
                } else {
                    DVec3::ZERO
                },
                Some((sr, pdf)) => {
                    let is_specular = match material {
                        Material::Mirror | Material::Glass => true,
                        _ => false,
                    };

                    shadow_ray(scene, &h, rand_utils::rand_unit_square())
                        + material.brdf(h.p)
                        * integrate(scene, &sr, depth + 1, is_specular)
                        * h.norm.dot(sr.dir.normalize()).abs()
                        /* something fishy here. pdf cancels one cosine... */
                        * h.norm.dot(sr.dir.normalize()).abs()
                        / ((1.0 - PATH_TRACE_RR)
                           * pdf)
                }
            }
        }
    }
}
