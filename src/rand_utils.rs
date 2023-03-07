use crate::DVec3;
use rand::Rng;
use rand::prelude::SliceRandom;

/* should figure better way to rng creation */

/* n normalized random dvec3 in a vec */
pub fn rand_vec_dvec3(n: usize) -> Vec<DVec3> {
    let mut rng = rand::thread_rng();
    (0..n).map(|_| {
        (DVec3::new(rng.gen(), rng.gen(), rng.gen()) - DVec3::splat(0.5))
            .normalize()
    }).collect()
}

/* random permutation of 0..n */
pub fn perm_n(n: usize) -> Vec<usize> {
    let mut v: Vec<usize> = (0..n).collect();
    v.shuffle(&mut rand::thread_rng());
    v
}
