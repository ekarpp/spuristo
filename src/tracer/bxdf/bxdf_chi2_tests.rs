use super::*;
use statrs::distribution::{ ChiSquared, ContinuousCDF };

const NUM_SAMPLES: usize = 1_000_000;
const THETA_BINS: usize = 80;
const PHI_BINS: usize = 2 * THETA_BINS;
const CHI2_RUNS: usize = 5;
const CHI2_SLEVEL: Float = 0.05;
const CHI2_MIN_FREQ: Float = 5.0;

fn mfd(roughness: Float) -> MfDistribution {
    MfDistribution::new(roughness, 1.5, 0.0, true)
}

#[test]
fn lambertian_chi2() {
    let bxdf = BxDF::Lambertian;
    for _ in 0..CHI2_RUNS {
        let wo = -rand_utils::square_to_cos_hemisphere(rand_utils::unit_square());
        assert!(chi2_pass(wo, &bxdf));
    }
}

#[test]
fn conductor50_chi2() {
    let mfd = mfd(0.5);
    let bxdf = BxDF::MfConductor(mfd);
    for _ in 0..CHI2_RUNS {
        let wo = -rand_utils::square_to_cos_hemisphere(rand_utils::unit_square());
        assert!(chi2_pass(wo, &bxdf));
    }
}

#[test]
fn conductor25_chi2() {
    let mfd = mfd(0.25);
    let bxdf = BxDF::MfConductor(mfd);
    for _ in 0..CHI2_RUNS {
        let wo = -rand_utils::square_to_cos_hemisphere(rand_utils::unit_square());
        assert!(chi2_pass(wo, &bxdf));
    }
}

#[test]
fn dielectric50_chi2() {
    let mfd = mfd(0.5);
    let bxdf = BxDF::MfDielectric(mfd);
    for _ in 0..CHI2_RUNS {
        let wo = -rand_utils::square_to_cos_hemisphere(rand_utils::unit_square());
        assert!(chi2_pass(wo, &bxdf));
    }
}

#[test]
fn dielectric25_chi2() {
    let mfd = mfd(0.25);
    let bxdf = BxDF::MfDielectric(mfd);
    for _ in 0..CHI2_RUNS {
        let wo = -rand_utils::square_to_cos_hemisphere(rand_utils::unit_square());
        assert!(chi2_pass(wo, &bxdf));
    }
}

fn chi2_pass(wo: Direction, bxdf: &BxDF) -> bool {
    let actual_freq = sample_frequencies(wo, &bxdf);
    let expected_freq = compute_frequencies(wo, &bxdf);

    // degrees of freedom
    let mut dof = 0;
    // test statistic = sum_i (O_i - E_i)^2 / E_i
    let mut stat = 0.0;

    let mut pooled_samples = 0;
    let mut pooled_computed = 0.0;

    for phi_bin in 0..PHI_BINS {
        for theta_bin in 0..THETA_BINS {
            let idx = phi_bin + theta_bin * PHI_BINS;

            if expected_freq[idx] == 0.0 {
                if actual_freq[idx] as Float > NUM_SAMPLES as Float * 1e-5 {
                    println!(
                        "Found sampled value of {} where expectation was zero.",
                        actual_freq[idx]
                    );
                    return false;
                }
            } else if expected_freq[idx] < CHI2_MIN_FREQ {
                pooled_samples += actual_freq[idx];
                pooled_computed += expected_freq[idx];
            } else if pooled_computed < CHI2_MIN_FREQ {
                // make sure the pooled is high enough
                pooled_samples += actual_freq[idx];
                pooled_computed += expected_freq[idx];
            } else {
                let delta = actual_freq[idx] as Float - expected_freq[idx];
                stat += (delta * delta) / expected_freq[idx];
                dof += 1;
            }
        }
    }

    if pooled_samples + pooled_computed as usize > 0 {
        let delta = pooled_samples as Float - pooled_computed;
        stat += (delta * delta) / pooled_computed as Float;
        dof += 1;
    }

    // we assume that all parameters are known in the pearson chi2 test of goodness of fit
    dof -= 1;

    if dof <= 0 {
        println!("Got 0 DoF!");
        false
    } else {
        // null hypothesis = our sampling follows the PDF

        // todo: code the CDF, should be simple but boring
        let chi2 = ChiSquared::new(dof as Float).unwrap();
        /* p-value, probability to get test statistic at least as extreme,
         * assuming null hypothesis holds
         */
        let pval = 1.0 - chi2.cdf(stat);
        println!("test statistic: {} p-value: {}", stat, pval);

        // we are possibly running multiple chi2 tests. apply Šidák correction
        let alpha = 1.0 - (1.0 - CHI2_SLEVEL).powf(1.0 / CHI2_RUNS as Float);

        pval >= alpha
    }
}

fn sample_frequencies(wo: Direction, bxdf: &BxDF) -> [usize; THETA_BINS*PHI_BINS] {
    let mut samples = [0; THETA_BINS*PHI_BINS];

    let theta_factor = THETA_BINS as Float / crate::PI;
    let phi_factor = PHI_BINS as Float / (2.0 * crate::PI);

    for _ in 0..NUM_SAMPLES {
        match bxdf.sample(wo, rand_utils::unit_square()) {
            None => (),
            Some(wi) => {
                let theta = spherical_utils::theta(wi);
                let phi = spherical_utils::phi(wi);

                let theta_bin = ((theta * theta_factor) as usize).max(0).min(THETA_BINS - 1);
                let phi_bin = ((phi * phi_factor) as usize).max(0).min(PHI_BINS - 1);

                samples[phi_bin + theta_bin * PHI_BINS] += 1;
            }
        }
    }

    return samples;
}

fn compute_frequencies(wo: Direction, bxdf: &BxDF) -> [Float; THETA_BINS*PHI_BINS] {
    let mut samples = [0.0; THETA_BINS*PHI_BINS];
    let mut ig = 0.0;
    let theta_factor = crate::PI / THETA_BINS as Float;
    let phi_factor = (2.0 * crate::PI) / PHI_BINS as Float;

    for theta_bin in 0..THETA_BINS {
        let theta0 = theta_bin as Float * theta_factor;
        let theta1 = theta0 + theta_factor;
        for phi_bin in 0..PHI_BINS {
            let phi0 = phi_bin as Float * phi_factor;
            let phi1 = phi0 + phi_factor;
            let f = |theta: Float, phi: Float| {
                let wi = Direction::new(
                    theta.sin() * phi.cos(),
                    theta.sin() * phi.sin(),
                    theta.cos(),
                );
                // pdf in solid angle, change to spherical coordinates
                bxdf.pdf(wo, wi, false) * theta.sin()
            };
            let integral = simpson2d(f, theta0, phi0, theta1, phi1);
            ig += integral;
            samples[phi_bin + theta_bin * PHI_BINS] = integral * NUM_SAMPLES as Float;
        }
    }
    println!("integral: {}", ig);
    return samples;
}

struct SimpsonIteration {
    a: Float,
    b: Float,
    c: Float,
    fa: Float,
    fb: Float,
    fc: Float,
    integrand: Float,
}

impl SimpsonIteration {
    pub fn new(a: Float, b: Float, c: Float, fa: Float, fb: Float, fc: Float, integrand: Float) -> Self {
        Self { a, b, c, fa, fb, fc, integrand }
    }
}

fn simpson2d<F>(f: F, x0: Float, y0: Float, x1: Float, y1: Float) -> Float where
    F: Fn(Float, Float) -> Float {
    let g = |y: Float| {
        let h = |x: Float| f(x, y);
        simpson(h, x0, x1)
    };
    simpson(g, y0, y1)
}

fn simpson<F>(f: F, x0: Float, x1: Float) -> Float where
    F: Fn(Float) -> Float {

    fn adaptive<G>(f: &G, it: SimpsonIteration, eps: Float, depth: usize) -> Float where
        G: Fn(Float) -> Float {
        let SimpsonIteration { a, b, c, fa, fb, fc, integrand } = it;
        // split to sub intervals
        let (d, e) = (0.5 * (a + b), 0.5 * (b + c));
        let (fd, fe) = (f(d), f(e));

        // integrate subintervals
        let h = c - a;
        let i0 = h * (fa + 4.0 * fd + fb) / 12.0;
        let i1 = h * (fb + 4.0 * fe + fc) / 12.0;
        let i01 = i0 + i1;

        if depth <= 0 || (i01 - integrand).abs() < 15.0 * eps {
            i01 + (i01 - integrand) / 15.0
        } else {
            let it1 = SimpsonIteration::new(a, d, b, fa, fd, fb, i0);
            let it2 = SimpsonIteration::new(b, e, c, fb, fe, fc, i1);

            adaptive(f, it1, 0.5 * eps, depth - 1) + adaptive(f, it2, 0.5 * eps, depth - 1)
        }
    }

    let (a, b, c) = (x0, 0.5 * (x0 + x1), x1);
    let (fa, fb, fc) = (f(a), f(b), f(c));
    let integrand = (fa + 4.0 * fb + fc) * (c - a) / 6.0;
    let it = SimpsonIteration { a, b, c, fa, fb, fc, integrand };

    adaptive(&f, it, 1e-6, 6)
}
