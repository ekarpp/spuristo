use super::*;

/* util functions */
mod util {
    use super::*;

    pub fn reflect(v: Direction, no: Normal) -> Direction {
        2.0 * v.project_onto(no) - v
    }

    pub fn refract(eta_ratio: Float, v: Direction, no: Normal) -> Option<Direction> {
        /* Snell-Descartes law */
        let cos_to = no.dot(v);
        let cos_to_abs = cos_to.abs();
        let sin2_to = 1.0 - cos_to * cos_to;
        let sin2_ti = eta_ratio * eta_ratio * sin2_to;

        if sin2_ti >= 1.0 {
            /* total internal reflection */
            // we don't do it?
            None
        } else {
            let cos_ti = (1.0 - sin2_ti).sqrt();
            let n = if v.z < 0.0 { -no } else { no };
            Some( -v * eta_ratio + (eta_ratio * cos_to_abs - cos_ti) * n )
        }
    }

    pub fn reflect_coeff(
        wo: Direction,
        wi: Direction,
        mfd: &MfDistribution,
    ) -> Float {
        let v = -wo;
        let cos_theta_v = v.z.abs();
        let cos_theta_wi = wi.z.abs();
        let wh = (wi + v).normalize();

        let d = mfd.d(wh);
        let f = mfd.f(v, wh);
        let g = mfd.g(v, wi, wh);

        // need reflection color, its in the .mtl files somewhere
        d * f * g / (4.0 * cos_theta_v * cos_theta_wi)
    }

}

/*
 * MICROFACET REFLECTION
 */
pub fn reflection_f(
    wo: Direction,
    wi: Direction,
    mfd: &MfDistribution,
    albedo: Color
) -> Color {
    let v = -wo;
    let wh = (wi + v).normalize();
    let f = mfd.f(v, wh);

    albedo * ((1.0 - f) / crate::PI + util::reflect_coeff(wo, wi, mfd))
}

pub fn reflection_sample(
    wo: Direction,
    mfd: &MfDistribution,
    rand_sq: Vec2
) -> Option<Direction> {
    let v = -wo;
    let wh = mfd.sample_normal(v, rand_sq).normalize();
    let wi = util::reflect(v, wh);

    if wi.z <= 0.0 {
        // bad sample, do something else?
        None
    } else {
        Some( wi )
    }
}

pub fn reflection_pdf(
    wo: Direction,
    wi: Direction,
    mfd: &MfDistribution,
) -> Float {
    let v = -wo;
    let wh = (v + wi).normalize();
    let wh_dot_v = v.dot(wh);

    mfd.sample_normal_pdf(wh, v) / (4.0 * wh_dot_v)
}

/*
 * MICROFACET DIFFUSE
 * Just lambertian
 */
pub fn diffuse_f(
    wo: Direction,
    wi: Direction,
    mfd: &MfDistribution,
    albedo: Color,
) -> Color {
    let v = -wo;
    let wh = (v + wi).normalize();

    /*
    let cos_theta_v = v.z;
    let cos_theta_wi = wi.z;
    let cos_theta_wh = wh.z;
    let disney = mfd.disney_diffuse(cos_theta_v, cos_theta_wi, cos_theta_wh);
    */
    let f = mfd.f(v, wh);

    albedo * ((1.0 - f) / crate::PI + util::reflect_coeff(wo, wi, mfd))
}

/*
 * MICROFACET TRANSMISSION
 */

pub fn transmission_f(
    wo: Direction,
    wi: Direction,
    mfd: &MfDistribution,
    albedo: Color,
    mode: Transport,
) -> Color {
    let v = -wo;
    let cos_theta_v = v.z;
    let cos_theta_wi = wi.z;

    /* on the same hemisphere */
    if cos_theta_v * cos_theta_wi > 0.0 {
        return Color::BLACK;
    }

    let eta_ratio = if cos_theta_v < 0.0 {
        1.0 / mfd.eta()
    } else {
        mfd.eta()
    };
    let scale = match mode {
        Transport::Radiance => eta_ratio * eta_ratio,
        Transport::Importance => 1.0,
    };
    let wh = (wi * eta_ratio + v).normalize();
    let wh = if wh.z < 0.0 { -wh } else { wh };

    let wh_dot_v = wh.dot(v);
    let wh_dot_wi = wh.dot(wi);

    /* same hemisphere w.r.t. wh */
    if wh_dot_v * wh_dot_wi > 0.0 {
        return Color::BLACK;
    }

    let d = mfd.d(wh);
    let f = mfd.f(v, wh);
    let g = mfd.g(v, wi, wh);

    albedo * (scale * d * (1.0 - f) * g
        * (wh_dot_wi * wh_dot_v / (cos_theta_wi * cos_theta_v)).abs()
        / (eta_ratio * wh_dot_wi + wh_dot_v).powi(2)
        + util::reflect_coeff(wo, wi, mfd))
}

pub fn transmission_sample(
    wo: Direction,
    mfd: &MfDistribution,
    rand_sq: Vec2
) -> Option<Direction> {
    let v = -wo;
    let inside = v.z < 0.0;
    // v needs to be "outside" for normal sampling
    let wh = if inside {
        mfd.sample_normal(-v, rand_sq).normalize()
    } else {
        mfd.sample_normal(v, rand_sq).normalize()
    };


    let eta_ratio = if inside {
        mfd.eta()
    } else {
        1.0 / mfd.eta()
    };

    util::refract(eta_ratio, v, wh)
}

pub fn transmission_pdf(
    wo: Direction,
    wi: Direction,
    mfd: &MfDistribution,
    swap_dir: bool,
) -> Float {
    let v = -wo;
    let (v, wi) = if swap_dir { (wi, v) } else { (v, wi) };
    let v_inside = v.z < 0.0;
    let wi_inside = wi.z < 0.0;

    if v_inside == wi_inside {
        /* same hemisphere */
        0.0
    } else {
        let eta_ratio = if v_inside {
            1.0 / mfd.eta()
        } else {
            mfd.eta()
        };

        let wh = -(v + wi * eta_ratio).normalize();
        let wh_dot_wi = wi.dot(wh);
        let wh_dot_v = wh.dot(v);

        if wh_dot_wi * wh_dot_v > 0.0 {
            /* same hemisphere w.r.t. wh */
            0.0
        } else {
            let wh = if wh.z < 0.0 { -wh } else { wh };
            mfd.sample_normal_pdf(wh, v)
                * (eta_ratio * eta_ratio * wh_dot_wi).abs()
                / (wh_dot_v + eta_ratio * wh_dot_wi).powi(2)
        }
    }
}
