use super::*;

pub fn integrate(scene: &Scene, ro: Ray, last_specular: bool) -> DVec3 {
    match scene.hit(&ro) {
        None => DVec3::ZERO,
        Some(ho) => {
            let material = ho.object.material();

            match material.bsdf_pdf(&ho, &ro) {
                None => {
                    if last_specular {
                        material.emit(&ho)
                    } else {
                        DVec3::ZERO
                    }
                }
                Some(scatter_pdf) => {
                    // jittered sampler
                    let shadow = if material.is_delta() {
                        DVec3::ZERO
                    } else {
                        JitteredSampler::new(SHADOW_SPLITS)
                            .map(|rand_sq| {
                                shadow_ray(
                                    scene,
                                    &ro,
                                    &ho,
                                    scatter_pdf.as_ref(),
                                    rand_sq
                                )
                            })
                        .sum::<DVec3>()
                        / SHADOW_SPLITS as f64
                    };

                    if rand_utils::rand_f64() < PATH_TRACE_RR {
                        return shadow;
                    }


                    match scatter_pdf.sample_ray(rand_utils::unit_square()) {
                        None => shadow,
                        Some(ri) => {
                            let wi = ri.dir;
                            let p_scatter = scatter_pdf.value_for(&ri);

                            // resample bad samples?
                            if p_scatter <= 0.0 {
                                return shadow;
                            }
                            let ng = ho.ng;
                            let ns = ho.ns;

                            shadow
                                + material.bsdf_f(&ro, &ri, ns, ng)
                                * ns.dot(wi).abs()
                                * integrate(scene, ri, material.is_specular())
                                / (p_scatter * (1.0 - PATH_TRACE_RR))
                        }
                    }
                }
            }
        }
    }
}
