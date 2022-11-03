pub fn centerize_normalize_boundingbox(vtx_xyz:Vec<f32>, ndim:usize) -> Vec<f32> {
    let mut a = nalgebra::DMatrix::from_vec(ndim, vtx_xyz.len()/ndim, vtx_xyz);
    let mut cnt = nalgebra::DVector::zeros(ndim);
    let mut len = nalgebra::DVector::zeros(ndim);
    for i in 0..ndim {
        let x0 = a.row(i).min();
        let x1 = a.row(i).max();
        cnt[i] = (x0+x1)*0.5_f32;
        len[i] = x1-x0;
    }
    a.column_iter_mut().for_each(|mut r| r -= cnt.clone());
    let size = len.max();
    a *= 1.0/size;
    a.as_slice().to_vec()
}