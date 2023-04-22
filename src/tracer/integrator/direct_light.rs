use super::*;

pub fn integrate(scene: &Scene, ro: Ray) -> DVec3 {
    _integrate(scene, ro, 0)
}

const MAX_RECURSION: usize = 50;

fn _integrate(scene: &Scene, ro: Ray, depth: usize) -> DVec3 {
    match scene.hit(&ro) {
        None => DVec3::new(0.0, 0.0, 0.0),
        Some(ho) => {
            let material = ho.object.material();
            match material.bsdf_pdf(&ho, &ro) {
                None => material.emit(&ho),
                Some(scatter_pdf) => {
                    if material.is_specular() {
                        if depth > MAX_RECURSION {
                            return DVec3::ZERO;
                        }

                        match scatter_pdf.sample_direction(rand_utils::unit_square()) {
                            None => DVec3::ZERO,
                            Some(wi) => {
                                let ri = ho.generate_ray(wi);
                                let wi = ri.dir;
                                let wo = ro.dir;

                                let p_scatter = scatter_pdf.value_for(&ri);

                                if p_scatter <= 0.0 {
                                    // return something better?
                                    return DVec3::ZERO;
                                }

                                let ns = ho.ns;

                                material.bsdf_f(wo, wi, &ho)
                                    * scene.transmittance(ho.t)
                                    * ns.dot(wi).abs()
                                    * _integrate(scene, ri, depth + 1)
                                    / p_scatter
                            }
                        }
                    } else {
                        JitteredSampler::new(SHADOW_SPLITS)
                            .fold(DVec3::ZERO, |sum, rand_sq| {
                                shadow_ray(
                                    scene,
                                    &ro,
                                    &ho,
                                    scatter_pdf.as_ref(),
                                    rand_sq
                                ) + sum
                            })
                            * scene.transmittance(ho.t)
                            / SHADOW_SPLITS as f64
                    }
                }
            }
        }
    }
}
