
pub fn intersection_meshtri3(
    ray_org: &[f32],
    ray_dir: &[f32],
    vtx_xyz: &[f32],
    tri_vtx: &[usize]) -> Option<([f32; 3],usize)> {
    use crate::geo_tri;
    let mut hit_pos = Vec::<(f32, usize)>::new();
    for itri in 0..tri_vtx.len() / 3 {
        let i0 = tri_vtx[itri * 3 + 0];
        let i1 = tri_vtx[itri * 3 + 1];
        let i2 = tri_vtx[itri * 3 + 2];
        let res = geo_tri::ray_triangle_intersection(
            &ray_org, &ray_dir,
            &vtx_xyz[i0 * 3 + 0..i0 * 3 + 3],
            &vtx_xyz[i1 * 3 + 0..i1 * 3 + 3],
            &vtx_xyz[i2 * 3 + 0..i2 * 3 + 3]);
        match res {
            None => { continue; }
            Some(t) => {
                hit_pos.push((t, itri));
            }
        }
    }
    if hit_pos.is_empty() { return None; }
    hit_pos.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let t = hit_pos[0].0;
    let a = [
        t * ray_dir[0] + ray_org[0],
        t * ray_dir[1] + ray_org[1],
        t * ray_dir[2] + ray_org[2] ];
    return Some((a, hit_pos[0].1));
}