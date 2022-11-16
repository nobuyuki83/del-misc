//! heat distance method


/// gradient on 3D triangle mesh
pub fn gradient_on_trimsh3(
    tri2gradxyz: &mut [f32],
    tri2vtx: &[usize],
    vtx2xyz: &[f32],
    vtx2heat: &[f32]) {
    assert_eq!(tri2gradxyz.len(), tri2vtx.len());
    for i_tri in 0..tri2vtx.len() / 3 {
        use del_geo::vec3::{sub, scale, cross, normalize};
        let i0 = tri2vtx[i_tri * 3 + 0];
        let i1 = tri2vtx[i_tri * 3 + 1];
        let i2 = tri2vtx[i_tri * 3 + 2];
        let p0 = &vtx2xyz[i0 * 3..i0 * 3 + 3];
        let p1 = &vtx2xyz[i1 * 3..i1 * 3 + 3];
        let p2 = &vtx2xyz[i2 * 3..i2 * 3 + 3];
        let (area, un) = del_geo::tri::area_and_unorm3(p0, p1, p2);
        let e01 = sub(p1, p0);
        let e12 = sub(p2, p1);
        let e20 = sub(p0, p2);
        let g0 = scale(&cross(&un, &e12), 0.5_f32 / area);
        let g1 = scale(&cross(&un, &e20), 0.5_f32 / area);
        let g2 = scale(&cross(&un, &e01), 0.5_f32 / area);
        tri2gradxyz[i_tri * 3 + 0] = g0[0] * vtx2heat[i0] + g1[0] * vtx2heat[i1] + g2[0] * vtx2heat[i2];
        tri2gradxyz[i_tri * 3 + 1] = g0[1] * vtx2heat[i0] + g1[1] * vtx2heat[i1] + g2[1] * vtx2heat[i2];
        tri2gradxyz[i_tri * 3 + 2] = g0[2] * vtx2heat[i0] + g1[2] * vtx2heat[i1] + g2[2] * vtx2heat[i2];
        normalize(&mut tri2gradxyz[i_tri * 3..i_tri * 3 + 3]);
    }
}

/// divergence on 3D triangle mesh
pub fn divergence_on_trimesh3(
    vtx2div: &mut [f32],
    tri2vtx: &[usize],
    vtx2xyz: &[f32],
    tri2gradxyz: &[f32]) {
    for i_tri in 0..tri2vtx.len() / 3 {
        use del_geo::vec3::{sub, dot};
        let i0 = tri2vtx[i_tri * 3 + 0];
        let i1 = tri2vtx[i_tri * 3 + 1];
        let i2 = tri2vtx[i_tri * 3 + 2];
        let p0 = &vtx2xyz[i0 * 3..i0 * 3 + 3];
        let p1 = &vtx2xyz[i1 * 3..i1 * 3 + 3];
        let p2 = &vtx2xyz[i2 * 3..i2 * 3 + 3];
        let cots = del_geo::tri::cot3(p0, p1, p2);
        let x = &tri2gradxyz[i_tri * 3..i_tri * 3 + 3];
        let t0 = cots[0] * dot(&sub(p2, p1), x);
        let t1 = cots[1] * dot(&sub(p0, p2), x);
        let t2 = cots[2] * dot(&sub(p1, p0), x);
        vtx2div[i0] += t2 - t1;
        vtx2div[i1] += t0 - t2;
        vtx2div[i2] += t1 - t0;
    }
}