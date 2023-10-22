use crate::tracer::{ onb::Onb, Color };
use crate::{ Normal, Direction, Float, Vec2 };

/// Configurable parameters for a microsurface
#[derive(Copy, Clone)]
pub struct MicrofacetConfig {
    /// Roughness of the surface (α) [0,1]
    pub roughness: Float,
    /// Refraction index of the material >= 1.0
    pub eta: Float,
    /// Ratio of how metallic the material is [0,1]
    pub metallicity: Float,
}

impl MicrofacetConfig {
    pub fn new(
        roughness: Float,
        eta: Float,
        metallicity: Float,
    ) -> Self {
        assert!((0.0..=1.0).contains(&roughness));
        assert!((0.0..=1.0).contains(&metallicity));
        assert!(eta >= 1.0);

        Self {
            roughness: roughness.max(1e-5),
            eta,
            metallicity,
        }
    }
}

/// Defines a distribution of normals for a microfacet. `Float` parameter is the
/// roughness (α) of the surface.
#[derive(Copy, Clone)]
pub enum MfDistribution {
    /// Walter et al. 2007
    Ggx(MicrofacetConfig),
    /// Beckmann et al. 1987
    Beckmann(MicrofacetConfig),
}

impl MfDistribution {
    pub fn new(
        roughness: Float,
        eta: Float,
        metallicity: Float,
    ) -> Self {
        Self::Ggx(MicrofacetConfig::new(roughness, eta, metallicity))
    }

    /// might need tuning, send ratio that emittance is multiplied with?
    pub fn is_specular(&self) -> bool {
        self.get_config().roughness < 0.01
    }

    /// Does the material have delta scattering distribution?
    pub fn is_delta(&self) -> bool {
        self.get_config().roughness < 1e-2
    }

    pub fn eta(&self) -> Float {
        self.get_config().eta
    }

    /// Get roughness from config
    pub fn roughness(&self) -> Float {
        self.get_config().roughness
    }

    /// Getter, better way to do this?
    fn get_config(&self) -> &MicrofacetConfig {
        match self {
            Self::Ggx(cfg) | Self::Beckmann(cfg) => cfg,
        }
    }

    /// Disney diffuse (Burley 2012) with renormalization to conserve energy
    /// as done in Frostbite (Lagarde et al. 2014)
    pub fn disney_diffuse(
        &self,
        cos_theta_v: Float,
        cos_theta_wh: Float,
        cos_theta_wi: Float
    ) -> Float {
        let roughness2 = self.roughness().powi(2);
        let energy_bias = 0.5 * roughness2;
        let fd90 = energy_bias + 2.0 * cos_theta_wh.powi(2) * roughness2;

        let view_scatter = 1.0 + (fd90 - 1.0) * (1.0 - cos_theta_v).powi(5);
        let light_scatter = 1.0 + (fd90 - 1.0) * (1.0 - cos_theta_wi).powi(5);

        let energy_factor = 1.0 + roughness2 * (1.0 / 1.51 - 1.0);

        view_scatter * light_scatter * energy_factor
    }

    /// The microfacet distribution function.
    ///
    /// # Distributions
    /// * Beckmann - exp(-tan^2(θ) / α^2) / (π * α^2 * cos^4(θ))
    /// * GGX - α^2 / (π * (cos^4(θ) * (α^2 - 1.0) + 1.0)^2)
    ///
    /// # Arguments
    /// * `wh` - Microsurface normal in shading space
    pub fn d(&self, wh: Normal) -> Float {
        match self {
            Self::Ggx(cfg) => {
                let cos2_theta = wh.z.powi(2);

                if cos2_theta < crate::EPSILON {
                    0.0
                } else {
                    let roughness2 = cfg.roughness * cfg.roughness;

                    roughness2
                        / (crate::PI * (1.0 - cos2_theta * (1.0 - roughness2)).powi(2))
                }
            }
            Self::Beckmann(cfg) => {
                let cos2_theta = wh.z.powi(2);

                if cos2_theta < crate::EPSILON {
                    0.0
                } else {
                    let roughness2 = cfg.roughness * cfg.roughness;
                    let tan2_theta = (1.0 - cos2_theta) / cos2_theta;

                    (-tan2_theta / roughness2).exp()
                        / (crate::PI * roughness2 * cos2_theta.powi(2))
                }
            }
        }
    }

    /// Fresnel term with Schlick's approximation
    /// # Arguments
    /// * `v`      - Direction to viewer in shading space
    /// * `wh`     - Microsurface normal in shading space
    /// * `albedo` - Albedo at the point of impact
    pub fn f_reflection(&self, v: Direction, wh: Normal, albedo: Color) -> Color {
        let eta = self.eta();
        let cos_v = v.dot(wh).abs();

        let sin2_to = 1.0 - cos_v * cos_v;
        let sin2_ti = sin2_to * eta * eta;

        if sin2_ti >= 1.0 && v.z < 0.0 && false {
            // total internal reflection
            return Color::WHITE;
        }

        let metallicity = self.get_config().metallicity;

        let f0 = (eta - 1.0) / (eta + 1.0);
        let f0 = Color::splat(f0 * f0).lerp(albedo, metallicity);

        f0 + (Color::WHITE - f0) * (1.0 - cos_v).powi(5)
    }

    pub fn f_transmission(&self, v: Direction, wh: Normal) -> Float {
        let cos_o = v.dot(wh);
        let inside = cos_o < 0.0;
        let (eta_o, eta_i) = if inside {
            (self.eta(), 1.0)
        } else {
            (1.0, self.eta())
        };
        let cos_o = cos_o.abs();
        let sin_o = (1.0 - cos_o * cos_o).max(0.0).sqrt();
        let sin_i = eta_o / eta_i * sin_o;

        // total internal reflection
        if sin_i >= 1.0 {
            return 1.0;
        }

        let cos_i = (1.0 - sin_i * sin_i).max(0.0).sqrt();

        let r_par = ( (eta_i * cos_o) - (eta_o * cos_i) )
                  / ( (eta_i * cos_o) + (eta_o * cos_i) );
        let r_per = ( (eta_o * cos_o) - (eta_i * cos_i) )
                  / ( (eta_o * cos_o) + (eta_i * cos_i) );

        0.5 * (r_par * r_par + r_per * r_per)
    }


    /// Shadow-masking term. Used to make sure that only microfacets that are
    /// visible from `v` direction are considered. Uses the method described
    /// in Chapter 8.4.3 of PBR due to Heitz et al. 2013.
    ///
    /// # Arguments
    /// * `v`  - View direction in shading space
    /// * `wi` - Direction of ray away from the point of impact in shading space
    /// * `wh` - Microsurface normal in shading space
    pub fn g(&self, v: Direction, wi: Direction, wh: Normal) -> Float {
        // signum to fix refraction
        let cos_theta_wh = wh.z;
        let cos_theta_v = v.z;
        let chi = cos_theta_wh.signum() * v.dot(wh) / cos_theta_v;
        if chi < crate::EPSILON {
            0.0
        } else {
            1.0 / (1.0 + self.lambda(v) + self.lambda(wi))
        }
    }

    pub fn g1(&self, v: Direction, wh: Normal) -> Float {
        // signum to fix refraction
        let cos_theta_wh = wh.z;
        let cos_theta_v = v.z;
        let chi = cos_theta_wh.signum() * v.dot(wh) / cos_theta_v;
        if chi < crate::EPSILON {
            0.0
        } else {
            1.0 / (1.0 + self.lambda(v))
        }
    }

    /// Lambda function used in the definition of the shadow-masking term.
    /// Beckmann with polynomial approximation and GGX exactly. PBR Chapter 8.4.3
    ///
    /// # Arguments
    /// * `w` - Direction to consider in shading space
    fn lambda(&self, w: Direction) -> Float {
        match self {
            Self::Ggx(cfg) => {
                let cos2_theta = w.z.powi(2);
                if cos2_theta < crate::EPSILON {
                    0.0
                } else {
                    let tan2_theta = (1.0 - cos2_theta) / cos2_theta;
                    let roughness2 = cfg.roughness * cfg.roughness;

                    ((1.0 + roughness2 * tan2_theta).sqrt() - 1.0) / 2.0
                }
            }
            Self::Beckmann(cfg) => {
                let cos2_theta = w.z.powi(2);
                if cos2_theta < crate::EPSILON {
                    0.0
                } else {
                    let tan2_theta = ((1.0 - cos2_theta) / cos2_theta).abs();
                    let a = 1.0 / (cfg.roughness * tan2_theta);

                    if a >= 1.6 {
                        0.0
                    } else {
                        (1.0 - 1.259 * a + 0.396 * a * a)
                            / (3.535 * a + 2.181 * a * a)
                    }
                }
            }
        }
    }

    /// Probability to do importance sampling from NDF. Estimate based on
    /// the Fresnel term.
    pub fn probability_ndf_sample(&self, albedo: Color) -> Float {
        let cfg = self.get_config();

        let f0 = (cfg.eta - 1.0) / (cfg.eta + 1.0);
        let f0 = f0 * f0;

        (1.0 - cfg.metallicity) * f0 + cfg.metallicity * albedo.mean()
    }

    /// Probability that `wh` got sampled. `wh` and `v` in shading space.
    pub fn sample_normal_pdf(
        &self,
        wh: Normal,
        v: Direction,
    ) -> Float {
        let pdf = match self {
            Self::Beckmann(..) => {
                let cos_theta_wh = wh.z;
                self.d(wh) * cos_theta_wh
            }
            Self::Ggx(..) => {
                let wh_dot_v = wh.dot(v);
                let cos_theta_v = v.z;

                self.g1(v, wh) * self.d(wh) * wh_dot_v / cos_theta_v
            }
        };

        pdf.max(0.0)
    }

    /// Sampling microfacet normals per distribution for importance sampling.
    /// `v` in shading space.
    pub fn sample_normal(&self, v: Direction, rand_sq: Vec2) -> Normal {
        match self {
            Self::Ggx(cfg) => {
                // Heitz 2018 or
                // https://schuttejoe.github.io/post/ggximportancesamplingpart2/

                let roughness = cfg.roughness;
                // Map the GGX ellipsoid to a hemisphere
                let v_stretch = Direction::new(
                    v.x * roughness,
                    v.y * roughness,
                    v.z
                ).normalize();

                // ONB basis of the hemisphere configuration
                let hemi_basis = Onb::new(v_stretch);

                // compute a point on the disk
                let a = 1.0 / (1.0 + v_stretch.z);
                let r = rand_sq.x.sqrt();
                let phi = if rand_sq.y < a {
                    crate::PI * rand_sq.y / a
                } else {
                    crate::PI + crate::PI * (rand_sq.y - a) / (1.0 - a)
                };

                let x = r * phi.cos();
                let y = if rand_sq.y < a {
                    r * phi.sin()
                } else {
                    r * phi.sin() * v_stretch.z
                };

                // compute normal in hemisphere configuration
                let wm = Normal::new(
                    x,
                    y,
                    (1.0 - x*x - y*y).max(0.0).sqrt(),
                );
                let wm = hemi_basis.to_world(wm);

                // move back to ellipsoid
                Normal::new(
                    roughness * wm.x,
                    roughness * wm.y,
                    wm.z.max(0.0)
                ).normalize()
            }
            Self::Beckmann(cfg) => {
                let roughness2 = cfg.roughness * cfg.roughness;
                let theta = (-roughness2 * (1.0 - rand_sq.y).ln()).sqrt().atan();
                let phi = 2.0 * crate::PI * rand_sq.x;

                Normal::new(
                    theta.sin() * phi.cos(),
                    theta.sin() * phi.sin(),
                    theta.cos(),
                )
            }
        }
    }
}
