#![allow(warnings)]
use super::*;

// this could use the scoped assignment from PBRT...
/// Computes the MIS weight for the chosen sample strategy. PBRT, what orig paper
pub fn mis_weight(
    light_path: &[Vertex],
    s: usize,
    camera_path: &[Vertex],
    t: usize,
    sampled_vertex: Option<Vertex>,
) -> f64 {
    // assert!(t != 0)
    // if `sampled_vertex.is_some()` then t == 1 XOR s == 1

    if s + t == 2 {
        return 1.0;
    }

    let map0 = |pdf: f64| {
        if pdf == 0.0 { 1.0 } else { pdf }
    };

    // max(1) is lazy... if s == 0 we never call ls anyways. TODO fix later
    let mut ls: Option<Vertex> = None;
    let mut ls_m: Option<Vertex> = None;
    let mut ct: Option<Vertex> = None;
    let mut ct_m: Option<Vertex> = None;

    // if (s == 1) a1 = { qs, sampled };
    if s == 1 {
        ls = sampled_vertex;
    // else if (t == 1) a1 = { pt, sampled };
    } else if t == 1 {
        ct = sampled_vertex;
    }

    // *qs = s > 0 ? &lightVertices[s - 1]  : nullptr;
    if s > 1 {
        ls = Some(light_path[s - 1].copy());
    }
    // *pt = t > 0 ? &cameraVertices[t - 1] : nullptr;
    if t > 1 {
        ct = Some(camera_path[t - 1].copy());
    }

//    ls = ls.map(|v| v.delta = true; v);
//    ct = ct.map(|v| v.delta = true; v);
    // Update reverse density of vertex p_{t - 1}
    if let Some(ref mut ct) = ct {
        ct.delta = false;

        if s == 0 {
            ct.pdf_bck = ct.pdf_light_origin();
        } else if let Some(ref ls) = ls {
            let pdf = if s == 1 {
                ls.pdf_light_leaving(ct)
            } else {
                let ls_m = &light_path[s - 2];
                ls.pdf_area(ls_m, ct, Transport::Importance)
            };

            ct.pdf_bck = pdf;
        };
    }

    // update reverse density of vertex $p_{t-2}$
    if let Some(ref ct) = ct && t > 1 {
        // ct_m, but issues with shadowing
        let curr = &camera_path[t - 2];

        let pdf = if s == 0 {
            ct.pdf_light_leaving(curr)
        } else {
            ct.pdf_area(ls.as_ref().unwrap(), curr, Transport::Importance)
        };

        let mut curr = curr.copy();
        curr.pdf_bck = pdf;
        ct_m = Some(curr);
    }

    // update reverse density of vertex p_{s-1}
    if let Some(ref mut ls) = ls && let Some(ref ct) = ct {
        let pdf = if let Some(ref ct_m) = ct_m {
            ct.pdf_area(ct_m, ls, Transport::Radiance)
        } else {
            // should call camera.pdf here
            1.0
        };

        ls.pdf_bck = pdf;
        ls.delta = false;
    }

    // update reverse density of vertex p_{s-2}
    if let Some(ref ls) = ls && let Some(ref ct) = ct && s > 1 {
        let curr = &light_path[s - 2];
        let pdf = ls.pdf_area(ct, curr, Transport::Radiance);
        let mut curr = curr.copy();
        curr.pdf_bck = pdf;
        ls_m = Some(curr);
    }


    /*
    // max(1) is lazy... if s == 0 we never call ls anyways. TODO fix later
    let ls = if s != 1 { &light_path[s.max(1) - 1] } else { sampled_vertex.as_ref().unwrap() };
    let ct = if t != 1 { &camera_path[t - 1] } else { sampled_vertex.as_ref().unwrap() };
    */


    let mut sum_ri = 0.0;
    let mut ri = 1.0;

    // applies the updated PDF for camera_last of the connection
    if t > 0 {
        let pdf_prev = if s == 0 {
            // assert!(ct.h.light.is_some());
            // probability for the origin. uniformly sampled on light surface
            if ct.is_delta() {
                0.0
            } else {
                ct.h.light.map_or(0.0, |light| 1.0 / light.area())
            }
            // check camera and light path initializations. fix the mis weights to be like there. double check they are correct.
        } else if s == 1 {
            if let Some(light) = ls.h.light {
                let xo = ls.h.p;
                let xi = ct.h.p;
                let wi = xi - xo;
                let ri = Ray::new(xo, wi);
                // normalized
                let wi = ri.dir;
                let ng = ls.h.ng;
                let (_, pdf_dir) = light.sample_leaving_pdf(&ri, ng);
                let ng = ct.h.ng;
                // convert solid angle to area
                pdf_dir * wi.dot(ng).abs() / xo.distance_squared(xi)
            } else {
                unreachable!();
                0.0
            }
        } else {
            let ls_m = &light_path[s - 2];
            ls.pdf_area(ls_m, ct, Transport::Importance)
        };

        ri *= map0(pdf_prev) / map0(ct.pdf_fwd);
        let ls_is_delta = s == 0 || ls.is_delta();
        if !ct.is_delta() && !ls_is_delta {
            sum_ri += ri;
        }
    }

    // applies the updated PDF for camera t - 2 using the connection
    if t > 1 {
        let ct_m = &camera_path[t - 2];
        let pdf_prev = if s == 0 {
            if let Some(light) = ct.h.light {
                let xo = ct.h.p;
                let xi = ct_m.h.p;
                let wi = xi - xo;
                let ri = Ray::new(xo, wi);
                // normalized
                let wi = ri.dir;
                let ng = ct.h.ng;
                let (_, pdf_dir) = light.sample_leaving_pdf(&ri, ng);
                let ng = ct_m.h.ng;
                // convert solid angle to area
                pdf_dir * wi.dot(ng).abs() / xo.distance_squared(xi)
            } else {
                0.0
            }
        } else {
            ct.pdf_area(ls, ct_m, Transport::Importance)
        };
        ri *= map0(pdf_prev) / map0(ct_m.pdf_fwd);
        if !ct.is_delta() && !ct_m.is_delta() {
            sum_ri += ri;
        }
    }

    // vertices in camera path
    for i in (1..t.max(2) - 2).rev() {
        ri *= map0(camera_path[i].pdf_bck) / map0(camera_path[i].pdf_fwd);
        if !camera_path[i].is_delta() && !camera_path[i - 1].is_delta() {
            sum_ri += ri;
        }
    }

    let mut ri = 1.0;

    // applies the updated PDF at light_last using the connection
    if s > 0 {
        let pdf_prev = if t == 1 {
            // probability that the direction got sampled from camera.
            // should be 1.0? yes.
            1.0
        } else {
            let ct_m = &camera_path[t - 2];
            ct.pdf_area(ct_m, ls, Transport::Radiance)
        };
        ri *= map0(pdf_prev) / map0(ls.pdf_fwd);
        if !ls.is_delta() && !(t == 1 || camera_path[t - 1].is_delta()) {
            sum_ri += ri;
        }
    }

    // applies the updated PDF at light_last using the connection
    if s > 1 {
        let ls_m = &light_path[s - 2];
        let pdf_prev = ls.pdf_area(ct, ls_m, Transport::Radiance);

        ri *= map0(pdf_prev) / map0(ls_m.pdf_fwd);
        if !ls.is_delta() && !ls_m.is_delta() {
            sum_ri += ri;
        }
    }

    // vertices in light path
    for i in (0..s.max(2) - 2).rev() {
        ri *= map0(light_path[i].pdf_bck) / map0(light_path[i].pdf_fwd);
        if !light_path[i].is_delta() && !light_path[(i - 1).max(0)].is_delta() {
            sum_ri += ri;
        }
    }

    let weight = 1.0 / (1.0 + sum_ri);
    #[cfg(debug_assertions)]
    if weight < 0.0 {
        println!("negative weight in BDPT MIS");
    }
    weight
}
