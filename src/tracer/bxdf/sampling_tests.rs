use super::*;
use crate::simpson_integration;

const NUM_BINS: usize = 10;
const NUM_SAMPLES: usize = 10_000_000;

// used for numerically integrating PDF over whole space
const THETA_BINS: usize = 80;
const PHI_BINS: usize = 2 * THETA_BINS;

#[test]
fn lambertian_sampling() {
    let bxdf = BxDF::Lambertian;

    assert_bins(do_sampling(bxdf));
}

#[test]
fn conductor75_sampling() {
    let mfd = MfDistribution::new(0.75, 1.0, 0.0, true);
    let bxdf = BxDF::MfConductor(mfd);

    assert_bins(do_sampling(bxdf));
}

#[test]
fn conductor50_sampling() {
    let mfd = MfDistribution::new(0.5, 1.0, 0.0, true);
    let bxdf = BxDF::MfConductor(mfd);

    assert_bins(do_sampling(bxdf));
}

#[test]
fn conductor25_sampling() {
    let mfd = MfDistribution::new(0.25, 1.0, 0.0, true);
    let bxdf = BxDF::MfConductor(mfd);

    assert_bins(do_sampling(bxdf));
}

#[test]
fn conductor10_sampling() {
    let mfd = MfDistribution::new(0.10, 1.0, 0.0, true);
    let bxdf = BxDF::MfConductor(mfd);

    assert_bins(do_sampling(bxdf));
}

fn print_bins(bins: &Vec<Vec<Float>>) {
    println!("bin values:");
    for i in 0..NUM_BINS {
        for j in 0..NUM_BINS {
            print!("{:.3} ", bins[i][j]);
        }
        println!("");
    }
}

fn assert_bins(bins: Vec<Vec<Float>>) {
    // print bins first for debugging
    print_bins(&bins);
    let mut sum = 0.0;
    println!("bin deltas:");
    for i in 0..NUM_BINS {
        for j in 0..NUM_BINS {
            sum += bins[i][j];
            let delta = (2.0 * crate::PI - bins[i][j]).abs();
            print!("{:.3} ", delta);
            // rather lax here...
            assert!(delta < 5e-2 * NUM_BINS as Float);
        }
        println!("");
    }

    sum /= (NUM_BINS * NUM_BINS) as Float;
    let delta = (2.0 * crate::PI - sum).abs();
    println!("delta: {}", delta);
    assert!(delta < 1e-2);
}

/* ripped from pbrt. take v = normal and sample NUM_SAMPLES times.
 * if sample is ok add 1 / pdf to a *bin*. each bin should converge to 2*PI.
 * split bins using spherical coordinates of the sampled vector.
 */
fn do_sampling(bxdf: BxDF) -> Vec<Vec<Float>> {
    let mut bins: Vec<Vec<Float>> = vec!();
    let mut num_failed: usize = 0;

    for i in 0..NUM_BINS {
        bins.push(vec!());
        for _ in 0..NUM_BINS {
            bins[i].push(0.0);
        }
    }
    // let v face directly the normal
    let wo = Direction::NEG_Z;
    for _ in 0..NUM_SAMPLES {
        let wi = bxdf.sample(wo, rand_utils::unit_square());
        match wi {
            None => num_failed += 1,
            Some(wi) => {
                let pdf = bxdf.pdf(wo, wi, false);
                if pdf == 0.0 {
                    num_failed += 1;
                    continue;
                }
                let wi_phi = spherical_utils::phi(wi) / (2.0 * crate::PI);
                let phi_idx = (wi_phi * NUM_BINS as Float) as usize;
                let wi_cos_theta = spherical_utils::cos_theta(wi);
                let theta_idx = (wi_cos_theta * NUM_BINS as Float) as usize;
                // can be out of bounds in rare cases. fix if it happens.
                bins[theta_idx][phi_idx] += 1.0 / pdf;
            }
        }
    }
    let good_samples = NUM_SAMPLES - num_failed;

    let integral = integrate_sphere(wo, &bxdf);
    println!("integral over whole space: {}", integral);
    /* scale bins properly based on samples and integral over whole space of PDF.
     * (PDF for sampling microfacets does not always integrate to 1 over whole space)
     */
    for i in 0..NUM_BINS {
        for j in 0..NUM_BINS {
            bins[i][j] *= integral * (NUM_BINS * NUM_BINS) as Float / good_samples as Float;
        }
    }
    println!("failed samples: {:.2} %", 100.0 * num_failed as Float / NUM_SAMPLES as Float);
    bins
}

fn integrate_sphere(wo: Direction, bxdf: &BxDF) -> Float {
    let mut integral = 0.0;

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
            integral += simpson_integration::simpson2d(f, theta0, phi0, theta1, phi1);
        }
    }

    integral
}
