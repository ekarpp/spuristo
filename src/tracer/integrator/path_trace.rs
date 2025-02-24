use super::*;

pub fn integrate(scene: &Scene, mut ro: Ray, raster_xy: Vec2) -> FilmSample {
    let mut last_specular = true;
    let mut radiance = Color::BLACK;
    let mut gathered = Color::WHITE;
    let mut depth = 0;

    while let Some(ho) = scene.hit(&ro) {
        let material = ho.material;
        gathered *= scene.transmittance(ho.t);
        let wo = ro.dir;

        match material.bsdf_sample(wo, &ho, rand_utils::unit_square()) {
            None => {
                if last_specular {
                    radiance += gathered * material.emit(&ho)
                }
                break;
            }
            Some(wi) => {
                if !material.is_delta() {
                    radiance += gathered
                        * shadow_ray(
                            scene,
                            &ro,
                            &ho,
                            rand_utils::unit_square(),
                            rand_utils::unit_square()
                        );
                }

                let ri = ho.generate_ray(wi);
                let wo = ro.dir;
                let wi = ri.dir;

                let p_scatter = material.bsdf_pdf(wo, wi, &ho, false);
                // resample bad sample?
                if p_scatter <= 0.0 {
                    break;
                }

                let bsdf = material.bsdf_f(wo, wi, Transport::Radiance, &ho);
                let bsdf = if ho.is_medium() {
                    // assume that mediums get sampled perfectly
                    // according to the BSDF and thus cancel out PDF
                    bsdf * p_scatter
                } else {
                    bsdf
                };

                let ns = ho.ns;
                gathered *= bsdf * material.shading_cosine(wi, ns)
                    / p_scatter;

                // russian roulette
                if depth > 3 {
                    let luminance = gathered.luminance();
                    let rr_prob = (1.0 - luminance).max(0.05);
                    if rand_utils::rand_float() < rr_prob {
                        break;
                    }
                    gathered /= 1.0 - rr_prob;
                }

                last_specular = material.is_specular();
                depth += 1;
                ro = ri;
            }
        }
    }

    FilmSample::new(radiance, raster_xy, false)
}
