use super::*;

pub fn integrate(scene: &Scene, ro: &Ray) -> DVec3 {
    /* two mirrors next to each other might cause issues... */

    match scene.hit(ro) {
        None => DVec3::new(0.0, 0.0, 0.0),
        Some(ho) => {
            let material = ho.object.material();
            match material.bsdf(ro, &ho) {
                None => material.emit(&ho),
                Some(bsdf) => {
                    //  ¯\_(ツ)_/¯
                    if material.specularity() > 0.92 {
                        let no = ho.norm;
                        let ri = bsdf.sample_ray(rand_utils::unit_square());
                        let wi = ri.dir;
                        println!("{}", wi);

                        // correct?
                        let cos_theta = if material.is_transparent() {
                            1.0
                        } else {
                            no.dot(wi).abs()
                        };

                        bsdf.eval(&ri) * cos_theta * integrate(scene, &ri)
                            / bsdf.prob_for(&ri)
                    } else {
                        JitteredSampler::new(SHADOW_SPLITS).fold(DVec3::ZERO, |acc, rand_sq| {
                            acc + shadow_ray(scene, ro, &ho, bsdf.as_ref(), rand_sq)
                        }) / SHADOW_SPLITS as f64
                    }
                }
            }
        }
    }
}
