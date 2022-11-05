

pub fn merge_trimesh3(
    sparse: &mut del_ls::sparse::BlockSparseMatrix<f32>,
    merge_buffer: &mut Vec<usize>,
    tri2vtx: &[usize],
    vtx2xyz: &[f32]) {
    for it in 0..tri2vtx.len() / 3 {
        let i0 = tri2vtx[it * 3 + 0];
        let i1 = tri2vtx[it * 3 + 1];
        let i2 = tri2vtx[it * 3 + 2];
        let cots = del_geo::tri::cot3(
            &vtx2xyz[(i0 * 3 + 0)..(i0 * 3 + 3)],
            &vtx2xyz[(i1 * 3 + 0)..(i1 * 3 + 3)],
            &vtx2xyz[(i2 * 3 + 0)..(i2 * 3 + 3)]);
        let emat: [f32; 9] = [
            cots[1] + cots[2], -cots[2], -cots[1],
            -cots[2], cots[2] + cots[0], -cots[0],
            -cots[1], -cots[0], cots[0] + cots[1]];
        sparse.merge(
            &[i0, i1, i2], &[i0, i1, i2], &emat,
            merge_buffer);
    }
}