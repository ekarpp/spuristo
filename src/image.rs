use crate::DVec3;
use std::path::Path;
use std::fs::File;
use std::io::BufWriter;
use png::{Encoder, ColorType, BitDepth, EncodingError};

/// Contains the necessary data to write buffer to file.
pub struct Image {
    /// Image buffer storing RGB-channels in range \[0,1\].
    pub buffer: Vec<DVec3>,
    /// Width of rendered image.
    pub width: i32,
    /// Height of rendered image.
    pub height: i32,

}

impl Image {
    pub fn new(buffer: Vec<DVec3>, width: i32, height: i32) -> Self {
        Self {
            buffer,
            width,
            height,
        }
    }
    /// Translates the image buffer of RGB values in range \[0,1\]
    /// to discrete range \[0,255\]. Applies gamma correction.
    #[allow(clippy::identity_op)]
    fn rgb(&self) -> Vec<u8> {
        self.buffer.iter().flat_map(|px: &DVec3| {
            [self.to_srgb(px.x), self.to_srgb(px.y), self.to_srgb(px.z)]
        }).collect()
    }

    fn to_srgb(&self, c: f64) -> u8 {
        (c.powf(1.0 / 2.2) * 255.0) as u8
    }

    /// Creates the PNG file
    pub fn save(&self, fname: &str) -> Result<(), EncodingError> {
        let path = Path::new(fname);

        let mut binding = BufWriter::new(File::create(path)?);
        let mut encoder = Encoder::new(
            &mut binding,
            self.width as u32,
            self.height as u32);
        encoder.set_color(ColorType::Rgb);
        encoder.set_depth(BitDepth::Eight);

        encoder.write_header()?.write_image_data(&self.rgb())?;
        Ok(())
    }
}
