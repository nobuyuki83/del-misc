
pub fn intersection_meshtri3(
    ray_org: &[f32],
    ray_dir: &[f32],
    vtx_xyz: &[f32],
    tri_vtx: &[usize]) -> Option<([f32; 3],usize)> {
    use del_geo::tri;
    let mut hit_pos = Vec::<(f32, usize)>::new();
    for itri in 0..tri_vtx.len() / 3 {
        let i0 = tri_vtx[itri * 3 + 0];
        let i1 = tri_vtx[itri * 3 + 1];
        let i2 = tri_vtx[itri * 3 + 2];
        let res = tri::ray_triangle_intersection(
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

fn indexes_of_connected_triangle_in_sphere(
    pos: [f32; 3],
    rad: f32,
    itri0: usize,
    vtx_xyz: &[f32],
    tri_vtx: &[usize],
    tri_adjtri: &[usize]) -> Vec<usize> {
    use del_geo::{tri, vec3};
    let mut res = Vec::<usize>::new();
    let mut searched = std::collections::BTreeSet::<usize>::new();
    let mut next0 = Vec::<usize>::new();
    next0.push(itri0);
    while let Some(iel0) = next0.pop() {
        if searched.contains(&iel0) { continue; } // already studied
        searched.insert(iel0);
        let dist_min = {
            let i0 = tri_vtx[iel0 * 3 + 0];
            let i1 = tri_vtx[iel0 * 3 + 1];
            let i2 = tri_vtx[iel0 * 3 + 2];
            let (pn, _r0, _r1) = tri::nearest_triangle3_point3(
                &pos,
                &vtx_xyz[i0 * 3..i0 * 3 + 3],
                &vtx_xyz[i1 * 3..i1 * 3 + 3],
                &vtx_xyz[i2 * 3..i2 * 3 + 3]);
            vec3::distance(&pn, &pos)
        };
        if dist_min > rad { continue; }
        res.push(iel0);
        for ie in 0..3 {
            let iel1 = tri_adjtri[iel0 * 3 + ie];
            if iel1 == usize::MAX { continue; }
            next0.push(iel1);
        }
    }
    res
}

pub fn is_there_point_on_mesh_inside_sphere(
    smpli: &(usize, f32, f32),
    rad: f32,
    samples: &Vec<(usize, f32, f32)>,
    el_smpl: &std::collections::HashMap<usize, Vec<usize> >,
    vtx_xyz: &Vec<f32>,
    tri_vtx: &Vec<usize>,
    tri_adjtri: &Vec<usize>) -> bool
{
    use del_msh::sampling;
    use del_geo::vec3;
    let pos_i = sampling::position_on_mesh_tri3(
        smpli.0, smpli.1, smpli.2,
        &vtx_xyz, &tri_vtx);
    let indexes_tri = indexes_of_connected_triangle_in_sphere(
        pos_i, rad,
        smpli.0, vtx_xyz, tri_vtx, tri_adjtri);
    for idx_tri in indexes_tri.iter() {
        if !el_smpl.contains_key(idx_tri) {
            continue;
        }
        for j_smpl in el_smpl[idx_tri].iter() {
            let smpl_j = samples[*j_smpl];
            let pos_j = sampling::position_on_mesh_tri3(
                smpl_j.0, smpl_j.1, smpl_j.2,
                &vtx_xyz, &tri_vtx);
            let dist = vec3::distance(&pos_i, &pos_j);
            if dist < rad { return true; }
        }
    }
    return false;
}