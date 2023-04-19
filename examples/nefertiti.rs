use lumo::tracer::*;
use lumo::*;
use std::f64::consts::PI;

// By CosmoWenmann (https://www.thingiverse.com/thing:3974391) licensed under CC BY-NC-SA 4.0
const IMAGE_FILE: &str = "examples/aem_aem21300_3dsl01_mo08-03_p_img.png";
const NEFE_URL: &str = "https://cdn.thingiverse.com/assets/c7/e1/b6/f6/12/SPK_Nefertiti_Scan_FOIA_Response.zip";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let camera = Camera::default();
    let def_color = srgb_to_linear(242, 242, 242);
    let silver = srgb_to_linear(192, 192, 192);

    let mut scene = Scene::empty_box(
        def_color,
        Material::diffuse(Texture::Solid(srgb_to_linear(255, 0, 0))),
        Material::diffuse(Texture::Solid(srgb_to_linear(0, 255, 0))),
    );

    scene.add(
        Mesh::new(
            obj::obj_from_url(NEFE_URL)?,
            Material::diffuse(Texture::Image(Image::from_file(IMAGE_FILE)?)),
        )
        .to_unit_size()
        .to_origin()
        .scale(0.9, 0.9, 0.9)
        .rotate_x(-PI / 2.0)
        .translate(0.0, -0.07, -1.45),
    );

    scene.add(
        Cube::new(Material::specular(
            Texture::Marble(Perlin::new(), silver),
            0.04,
        ))
        .translate(-0.5, -0.5, -0.5)
        .scale(0.4, 1.2, 0.4)
        .translate(0.0, -1.1, -1.45),
    );

    let renderer = Renderer::new(scene, camera);
    renderer.render().save("nefe.png")?;

    Ok(())
}
