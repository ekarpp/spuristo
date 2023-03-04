extern crate glam;

mod image;
mod tracer;

const WIDTH: usize = 3840;
const HEIGHT: usize = 2160;

fn main() {
    let mut image = image::Image {
        buffer: Vec::with_capacity(WIDTH*HEIGHT),
        width: WIDTH,
        height: HEIGHT
    };

    let cam = tracer::camera::def();
    let scene = tracer::scene::def();
    for y in (0..HEIGHT).rev()
    {
        for x in 0..WIDTH
        {
            let u = x as f32 / (WIDTH-1) as f32;
            let v = y as f32 / (HEIGHT-1) as f32;
            let r = cam.ray_at(u, v);
            image.buffer.push(r.color(&scene));
        }
    }
    image.save();
}
