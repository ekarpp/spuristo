use super::*;
use crate::tracer::texture::Texture;

const NUM_RAYS: usize = 10000;

fn get_mat() -> Material {
    Material::diffuse(Texture::Solid(DVec3::ONE))
}

#[test]
fn no_self_intersect() {
    let mesh = TriangleMesh {
        vertices: vec![DVec3::ZERO, DVec3::X, DVec3::X + DVec3::NEG_Z],
        normals: vec![],
        uvs: vec![],
    };

    let t = Triangle::new(Arc::new(mesh), (0, 1, 2), None, None);

    let r = Ray::new(DVec3::new(0.5, 0.0, -0.25), DVec3::Y);
    assert!(t.hit(&r, 0.0, INFINITY).is_none());
}

#[test]
fn no_hit_behind() {
    let mesh = TriangleMesh {
        vertices: vec![DVec3::ZERO, DVec3::X, DVec3::X + DVec3::NEG_Z],
        normals: vec![],
        uvs: vec![],
    };

    let t = Triangle::new(Arc::new(mesh), (0, 1, 2), None, None);

    let r = Ray::new(DVec3::new(0.5, 0.1, -0.25), DVec3::Y);
    assert!(t.hit(&r, 0.0, INFINITY).is_none());
}

#[test]
fn sampled_rays_hit() {
    let mesh = TriangleMesh {
        vertices: vec![DVec3::ZERO, DVec3::X, DVec3::X + DVec3::Y],
        normals: vec![],
        uvs: vec![],
    };

    let tri = Triangle::new(Arc::new(mesh), (0, 1, 2), None, None);

    let xo = 5.0 * DVec3::Z;

    for _ in 0..NUM_RAYS {
        let wi = tri.sample_towards(xo, rand_utils::unit_square());
        let ri = Ray::new(xo, wi);
        let (p, _) = tri.sample_towards_pdf(&ri);

        assert!(p > 0.0);
    }
}
