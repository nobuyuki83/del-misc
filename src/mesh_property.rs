pub fn area_par_vertex_in_triangles_mesh(
    tri2vtx: &[usize],
    vtx2xyz: &[f32]) -> Vec<f32> {
    let num_vtx = vtx2xyz.len() / 3;
    let num_tri = tri2vtx.len() / 3;
    let mut areas = vec!(0_f32; num_vtx);
    for i_tri in 0..num_tri {
        let i0 = tri2vtx[i_tri * 3 + 0];
        let i1 = tri2vtx[i_tri * 3 + 1];
        let i2 = tri2vtx[i_tri * 3 + 2];
        let a0 = del_geo::tri3::area_(
            &vtx2xyz[i0 * 3..i0 * 3 + 3],
            &vtx2xyz[i1 * 3..i1 * 3 + 3],
            &vtx2xyz[i2 * 3..i2 * 3 + 3]) / 3_f32;
        areas[i0] += a0;
        areas[i1] += a0;
        areas[i2] += a0;
    }
    areas
}

pub fn mean_edge_length_triangles_mesh(
    tri2vtx: &[usize],
    vtx2xyz: &[f32]) -> f32 {
    let num_tri = tri2vtx.len() / 3;
    let mut sum = 0_f32;
    for i_tri in 0..num_tri {
        let i0 = tri2vtx[i_tri * 3 + 0];
        let i1 = tri2vtx[i_tri * 3 + 1];
        let i2 = tri2vtx[i_tri * 3 + 2];
        sum += del_geo::vec3::distance_(
            &vtx2xyz[i0 * 3..i0 * 3 + 3],
            &vtx2xyz[i1 * 3..i1 * 3 + 3]);
        sum += del_geo::vec3::distance_(
            &vtx2xyz[i1 * 3..i1 * 3 + 3],
            &vtx2xyz[i2 * 3..i2 * 3 + 3]);
        sum += del_geo::vec3::distance_(
            &vtx2xyz[i2 * 3..i2 * 3 + 3],
            &vtx2xyz[i0 * 3..i0 * 3 + 3]);
    }
    sum / (num_tri * 3) as f32
}
