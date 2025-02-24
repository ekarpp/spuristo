use lumo::tracer::*;
use lumo::*;

// By CosmoWenmann (https://www.thingiverse.com/thing:3974391) licensed under CC BY-NC-SA 4.0
const TEX_FILE: &str = "aem_aem21300_3dsl01_mo08-03_p_img.png";
const NEFE_URL: &str = "https://cdn.thingiverse.com/assets/c7/e1/b6/f6/12/SPK_Nefertiti_Scan_FOIA_Response.zip";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut scene = Scene::default();

    /* floor */
    scene.add(Plane::new(
        Vec3::NEG_Y,
        Vec3::Y,
        Material::diffuse(Texture::Solid(Color::BLACK))
    ));

    /* roof */
    scene.add(Plane::new(
        Vec3::Y,
        Vec3::NEG_Y,
        Material::diffuse(Texture::Solid(Color::BLACK))
    ));

    /* left wall */
    scene.add(Plane::new(
        Vec3::NEG_X,
        Vec3::X,
        Material::diffuse(Texture::Solid(Color::BLACK)),
    ));

    /* right wall */
    scene.add(Plane::new(
        Vec3::X,
        Vec3::NEG_X,
        Material::diffuse(Texture::Solid(Color::BLACK))
    ));

    /* front */
    scene.add(Plane::new(
        2.0 * Vec3::NEG_Z,
        Vec3::Z,
        Material::diffuse(Texture::Solid(Color::BLACK))
    ));

    /* back */
    scene.add(Plane::new(
        Vec3::Z,
        Vec3::NEG_Z,
        Material::diffuse(Texture::Solid(Color::BLACK))
    ));

    /* bust */
    scene.add(
        Cube::new(Material::metal(
            Texture::Solid(Color::new(61, 45, 36)),
            0.0,
            2.5,
            0.0,
        ))
        .translate(-0.5, -0.5, -0.5)
        .scale(0.45, 0.5, 0.45)
        .translate(0.0, -0.75, -1.45),
    );

    /* statue */
    if cfg!(debug_assertions) {
	scene.add(
            Cylinder::new(
		0.6,
		0.1,
		Material::diffuse(Texture::Solid(Color::new(255, 0, 0))))
		.translate(0.0, -0.5, -1.45)
	);
    } else {
	scene.add(
	    parser::mesh_from_url(
                NEFE_URL,
		Material::microfacet(
                    Texture::Image(parser::texture_from_url(NEFE_URL, TEX_FILE)?),
                    0.8,
                    1.5,
                    0.0,
                    false,
                    true,
                )
            )?
		.to_unit_size()
		.to_origin()
		.scale(0.5, 0.5, 0.5)
		.rotate_x(-PI / 2.0)
		.translate(0.0, -0.25, -1.45),
	);
    }

    let xy_rect = Mat3::from_cols(
        Vec3::ZERO,
        Vec3::X,
        Vec3::X + Vec3::Y,
    );

    let theta = PI / 4.0;

    /* light */
    let disk_towards = Vec3::new(-0.046, -0.192, -1.314);
    let disk_dir = Vec3::new(-0.9132592442858147, 0.38194206881899756, -0.05342207528313975);

    // left, w.r.t camera
    scene.add_light(Disk::new(
	disk_towards + 0.3 * disk_dir,
	-disk_dir,
	0.05,
	Material::Light(Texture::Solid(30.0 * Color::WHITE))
    ));

    let right_disk_origin = Vec3::new(0.6, 0.15, -1.6);
    scene.add_light(Disk::new(
	right_disk_origin,
	disk_towards + 0.28 * Vec3::X - right_disk_origin,
	0.08,
	Material::Light(Texture::Solid(20.0 * Color::WHITE))
    ));

    scene.add_light(Rectangle::new(
        xy_rect,
        Material::Light(Texture::Solid(2.0 * Color::WHITE)))
                    .scale(0.4, 0.4, 1.0)
                    .rotate_y(-theta)
		    .rotate_axis(Vec3::new(theta.cos(), 0.0, theta.sin()), PI / 8.0)
                    .translate(0.6, 0.2, -1.7)
    );

    // behind
    scene.add_light(Rectangle::new(
        xy_rect,
        Material::Light(Texture::Solid(Color::WHITE)))
                    .scale(0.3, 0.3, 1.0)
                    .rotate_x(2.0 * PI - 2.0 * theta)
                    .translate(-0.15, 0.5, -0.8)
    );
    scene.add_light(Rectangle::new(
	xy_rect,
	Material::Light(Texture::Solid(2.0 * Color::WHITE)))
		    .scale(0.3, 0.3, 1.0)
		    .rotate_x(PI)
		    .translate(-0.1, 0.0, 0.0)
    );

    // above
    scene.add_light(Rectangle::new(
        xy_rect,
        Material::Light(Texture::Solid(2.0 * Color::WHITE)))
                    .scale(0.4, 0.4, 1.0)
                    .rotate_x(PI / 2.0)
                    .translate(-0.2, 0.5, -1.5)
    );

    let camera = if cfg!(debug_assertions) {
	Camera::perspective(
	    0.5 * Vec3::Z,
	    Vec3::NEG_Z,
	    Vec3::Y,
	    90.0,
            0.0,
            0.0,
	    1000,
	    1000
	)
    } else {
	Camera::orthographic(
            Vec3::new(0.12, -0.23, -1.205),
            Vec3::new(0.0, -0.26, -1.45),
            Vec3::Y,
            0.2,
            0.0,
            0.0,
            641,
            939,
	)
    };

    let renderer = Renderer::new(scene, camera);
    renderer.render().save("nefe.png")?;

    Ok(())
}
