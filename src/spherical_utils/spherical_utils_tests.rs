use super::*;
use crate::rand_utils;

const NUM_SAMPLES: usize = 100_000;

#[test]
fn cos_phi_test() {
    for _ in 0..NUM_SAMPLES {
        let w = rand_utils::square_to_sphere(rand_utils::unit_square()).normalize();
        assert!((cos_phi(w) - phi(w).cos()) < crate::EPSILON);
    }
}

#[test]
fn sin_phi_test() {
    for _ in 0..NUM_SAMPLES {
        let w = rand_utils::square_to_sphere(rand_utils::unit_square()).normalize();
        assert!((sin_phi(w) - phi(w).sin()) < crate::EPSILON);
    }
}

#[test]
fn cos_theta_test() {
    for _ in 0..NUM_SAMPLES {
        let w = rand_utils::square_to_sphere(rand_utils::unit_square()).normalize();
        assert!((cos_theta(w) - theta(w).cos()) < crate::EPSILON);
    }
}

#[test]
fn cos2_theta_test() {
    for _ in 0..NUM_SAMPLES {
        let w = rand_utils::square_to_sphere(rand_utils::unit_square()).normalize();
        assert!((cos2_theta(w) - theta(w).cos().powi(2)) < crate::EPSILON);
    }
}

#[test]
fn sin_theta_test() {
    for _ in 0..NUM_SAMPLES {
        let w = rand_utils::square_to_sphere(rand_utils::unit_square()).normalize();
        assert!((sin_theta(w) - theta(w).sin()) < crate::EPSILON);
    }
}

#[test]
fn sin2_theta_test() {
    for _ in 0..NUM_SAMPLES {
        let w = rand_utils::square_to_sphere(rand_utils::unit_square()).normalize();
        assert!((sin2_theta(w) - theta(w).sin().powi(2)) < crate::EPSILON);
    }
}

#[test]
fn tan_theta_test() {
    for _ in 0..NUM_SAMPLES {
        let w = rand_utils::square_to_cos_hemisphere(rand_utils::unit_square()).normalize();
        if cos_theta(w) >= crate::EPSILON {
            assert!((tan_theta(w) - theta(w).tan()) < 10.0 * crate::EPSILON);
        }
    }
}

#[test]
fn tan_nan_theta_test() {
    assert!(tan_theta(Direction::X).is_nan());
    assert!(tan_theta(Direction::Y).is_nan());
    assert!(tan_theta(Direction::splat(0.999 * crate::EPSILON)).is_nan());
}

#[test]
fn tan2_theta_test() {
    for _ in 0..NUM_SAMPLES {
        let w = rand_utils::square_to_cos_hemisphere(rand_utils::unit_square()).normalize();
        if cos2_theta(w) >= crate::EPSILON {
            assert!((tan2_theta(w) - theta(w).tan().powi(2)) < 100.0 * crate::EPSILON);
        }
    }
}

#[test]
fn tan2_nan_theta_test() {
    assert!(tan2_theta(Direction::X).is_nan());
    assert!(tan2_theta(Direction::Y).is_nan());
    assert!(tan2_theta(Direction::splat(0.999 * crate::EPSILON.sqrt())).is_nan());
}
