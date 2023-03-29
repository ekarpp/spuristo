use super::*;

pub fn integrate(scene: &Scene, ro: &Ray, last_specularity: f64) -> DVec3 {
    match scene.hit(ro) {
        None => DVec3::ZERO,
        Some(ho) => {
            let material = ho.object.material();

            match material.bsdf(ro, &ho) {
                None => last_specularity * material.emit(&ho),
                Some(bsdf) => {
                    // jittered sampler
                    let shadow = JitteredSampler::new(SHADOW_SPLITS)
                        .map(|rand_sq| shadow_ray(scene, ro, &ho, bsdf.as_ref(), rand_sq))
                        .sum::<DVec3>()
                        / SHADOW_SPLITS as f64;

                    if rand_utils::rand_f64() < PATH_TRACE_RR {
                        return shadow;
                    }

                    let no = ho.norm;
                    let ri = bsdf.sample_ray(rand_utils::unit_square());
                    let wi = ri.dir;
                    let p_scatter = bsdf.prob_for(&ri);

                    // scatter with 0 probability
                    if p_scatter.is_nan() || p_scatter < 0.0001 {
                        return shadow;
                    }

                    // correct?
                    let cos_theta = if material.is_transparent() {
                        1.0
                    } else {
                        no.dot(wi).abs()
                    };

                    shadow + bsdf.eval(&ri)
                        * cos_theta
                        * integrate(scene, &ri, material.specularity())
                        / (p_scatter * (1.0 - PATH_TRACE_RR))
                }
            }
        }
    }
}
