use super::*;

/// Abstraction of a vertex in the paths
pub struct Vertex<'a> {
    pub h: Hit<'a>,
    pub gathered: DVec3,
    pub pdf_fwd: f64,
    pub pdf_bck: f64,
    pub delta: bool,
}

impl<'a> Vertex<'a> {
    pub fn copy(&self) -> Self {
        let h = Hit::new(
            0.0,
            &Material::Blank,
            DVec3::NEG_X,
            DVec3::ZERO,
            DVec3::ZERO,
            DVec3::X,
            DVec3::X,
            DVec2::X,
        ).unwrap();
        Self {
            h,
            gathered: self.gathered,
            pdf_fwd: self.pdf_fwd,
            pdf_bck: self.pdf_bck,
            delta: self.delta,
        }
    }

    /// Camera vertex
    pub fn camera(xo: DVec3, gathered: DVec3) -> Self {
        let h = Hit::new(
            0.0,
            &Material::Blank,
            DVec3::NEG_X,
            xo,
            DVec3::ZERO,
            DVec3::X,
            DVec3::X,
            DVec2::X,
        ).unwrap();
        Self {
            h,
            gathered,
            pdf_bck: 1.0,
            pdf_fwd: 0.0,
            delta: true,
        }
    }

    /// Light vertex
    pub fn light(mut h: Hit<'a>, light: &'a dyn Sampleable, gathered: DVec3, pdf_bck: f64) -> Self {
        h.light = Some(light);
        Self {
            h,
            gathered,
            pdf_bck,
            // this might cause issues later on (if area is zero, point light?)...
            pdf_fwd: 1.0 / light.area(),
            // same as above ^
            delta: false,
        }
    }

    /// Surface vertex
    pub fn surface(
        h: Hit<'a>,
        gathered: DVec3,
        pdf_fwd: f64,
        prev: &Vertex,
    ) -> Self {
        let pdf_fwd = if h.material.is_delta() {
            0.0
        } else {
            let xo = prev.h.p;
            let xi = h.p;
            let wi = (xi - xo).normalize();
            let ng = h.ng;

            pdf_fwd * wi.dot(ng).abs() / xi.distance_squared(xo)
        };
        Self {
            h,
            gathered,
            pdf_fwd,
            pdf_bck: 0.0,
            delta: h.material.is_delta(),
        }
    }

    fn material(&self) -> &Material {
        self.h.material
    }

    /// Are we a surface/light vertex?
    pub fn is_surface(&self) -> bool {
        !matches!(self.material(), Material::Blank)
    }

    /// Are we on a light?
    pub fn is_light(&self) -> bool {
        self.h.is_light()
    }

    /// Are we on a surface with delta material?
    pub fn is_delta(&self) -> bool {
        self.material().is_delta()
    }

    /// Helper to get emittance at hit
    pub fn emittance(&self) -> DVec3 {
        self.material().emit(&self.h)
    }

    /// Helper to get shading cosine at hit
    pub fn shading_cosine(&self, wi: DVec3, ns: DVec3) -> f64 {
        self.material().shading_cosine(wi, ns)
    }

    /// Computes BSDF at hit of `self`
    pub fn bsdf(&self, prev: &Vertex, next: &Vertex, mode: Transport) -> DVec3 {
        // TODO (2)
        let wo = (self.h.p - prev.h.p).normalize();
        let wi = (next.h.p - self.h.p).normalize();

        self.material().bsdf_f(wo, wi, mode, &self.h)
    }

    /// Converts solid angle `pdf` to area PDF
    pub fn solid_angle_to_area(&self, pdf: f64, next: &Vertex) -> f64 {
        let xo = self.h.p;
        let xi = next.h.p;
        let wi = (xi - xo).normalize();

        if next.is_surface() {
            let ng = next.h.ng;
            pdf * wi.dot(ng).abs() / xi.distance_squared(xo)
        } else {
            pdf / xi.distance_squared(xo)
        }
    }

    /// Geometry term btwn `self` and `v` ... ...
    /// it is symmetric??
    pub fn g(&self, v: &Vertex) -> f64 {
        let xo = self.h.p;
        let xi = v.h.p;
        let no = self.h.ns;
        let ni = v.h.ns;

        let wi = (xi - xo).normalize();

        no.dot(wi).abs() * ni.dot(wi).abs() / xo.distance_squared(xi)
    }

    /// PDF to sample direction to `next` from `curr` w.r.t. surface area measure
    pub fn pdf_area(&self, prev: &Vertex, next: &Vertex, mode: Transport) -> f64 {
        let ho = &self.h;
        // prev
        let xo = prev.h.p;
        // curr
        let xi = ho.p;
        // prev -> curr
        let wo = xi - xo;
        let ro = Ray::new(xo, wo);
        // next
        let xii = next.h.p;
        // curr -> next
        let wi = xii - xi;
        let ri = Ray::new(xi, wi);
        // normalized
        let wi = ri.dir;
        let angle_pdf = match self.material().bsdf_pdf(ho, &ro) {
            None => 0.0,
            Some(pdf) => pdf.value_for(&ri, matches!(mode, Transport::Importance))
        };
        let ng = next.h.ng;

        angle_pdf * wi.dot(ng).abs() / xi.distance_squared(xii)
    }

    pub fn pdf_light_origin(&self) -> f64 {
        self.h.light.map_or(0.0, |light| 1.0 / light.area())
    }

    pub fn pdf_light_leaving(&self, next: &Vertex) -> f64 {
        if let Some(ref light) = self.h.light {
            let xo = self.h.p;
            let xi = next.h.p;
            let wi = xi - xo;
            let ri = Ray::new(xo, wi);
            // normalized
            let wi = ri.dir;
            let ng = self.h.ng;
            let (_, pdf_dir) = light.sample_leaving_pdf(&ri, ng);
            let ng = next.h.ng;
            // convert solid angle to area
            pdf_dir * wi.dot(ng).abs() / xo.distance_squared(xi)
        } else {
            0.0
        }
    }
}
