use num_traits::AsPrimitive;

pub fn area2<T>(p0: &[T], p1: &[T], p2: &[T]) -> T
    where T: std::ops::Sub<Output=T> + std::ops::Mul<Output=T> + 'static + Copy,
          f32: AsPrimitive<T>
{
    assert!(p0.len() == 2 && p1.len() == 2 && p2.len() == 2);
    0.5_f32.as_() * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))
}

pub fn dldx2(
    p0: &[f32],
    p1: &[f32],
    p2: &[f32]) -> [[f32; 3]; 3] {
    assert!(p0.len() == 2 && p1.len() == 2 && p2.len() == 2);
    let a0 = area2(p0, p1, p2);
    let tmp1 = 0.5 / a0;
    [
        [
            tmp1 * (p1[1] - p2[1]),
            tmp1 * (p2[1] - p0[1]),
            tmp1 * (p0[1] - p1[1]),
        ],
        [
            tmp1 * (p2[0] - p1[0]),
            tmp1 * (p0[0] - p2[0]),
            tmp1 * (p1[0] - p0[0]),
        ],
        [
            tmp1 * (p1[0] * p2[1] - p2[0] * p1[1]),
            tmp1 * (p2[0] * p0[1] - p0[0] * p2[1]),
            tmp1 * (p0[0] * p1[1] - p1[0] * p0[1]),
        ]
    ]
}

// above: 2D triangle
// --------------------
// below: 3D triangle

pub fn area3<T>(p0: &[T], p1: &[T], p2: &[T]) -> T
    where T: num_traits::real::Real + 'static,
          f32: num_traits::AsPrimitive<T>
{
    use crate::vec3;
    assert!(p0.len() == 3 && p1.len() == 3 && p2.len() == 3);
    let v1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let v2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let na = [
        v1[1] * v2[2] - v2[1] * v1[2],
        v1[2] * v2[0] - v2[2] * v1[0],
        v1[0] * v2[1] - v2[0] * v1[1]];
    return vec3::squared_norm(&na).sqrt() * 0.5_f32.as_();
}

pub fn normal3<T>(
    vnorm: &mut [T],
    v1: &[T],
    v2: &[T],
    v3: &[T])
    where T: std::ops::Sub<Output=T> + std::ops::Mul<Output=T> + std::ops::Sub + Copy
{
    vnorm[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1]);
    vnorm[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v2[0] - v1[0]) * (v3[2] - v1[2]);
    vnorm[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0]);
}

pub fn unit_normal3<T>(
    n: &mut [T],
    v1: &[T],
    v2: &[T],
    v3: &[T]) -> T
    where T: std::ops::Sub<Output=T> + std::ops::Mul<Output=T> + num_traits::Float + 'static + Copy + std::ops::MulAssign,
          f32: num_traits::AsPrimitive<T>
{
    use crate::vec3;
    normal3(
        n,
        v1, v2, v3);
    let a = vec3::norm(n) * 0.5_f32.as_();
    let invlen: T = 0.5_f32.as_() / a;
    n[0] *= invlen;
    n[1] *= invlen;
    n[2] *= invlen;
    a
}

pub fn cot3<T>(
    p0: &[T],
    p1: &[T],
    p2: &[T]) -> [T; 3]
    where T: num_traits::real::Real + 'static,
          f32: num_traits::AsPrimitive<T>
{
    use crate::vec3;
    assert!(p0.len() == 3 && p1.len() == 3 && p2.len() == 3);
    let v0 = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]];
    let v1 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let v2 = [p0[0] - p1[0], p0[1] - p1[1], p0[2] - p1[2]];
    let na = [
        v1[1] * v2[2] - v2[1] * v1[2],
        v1[2] * v2[0] - v2[2] * v1[0],
        v1[0] * v2[1] - v2[0] * v1[1]];
    let area: T = vec3::squared_norm(&na).sqrt() * 0.5_f32.as_();
    let tmp: T = 0.25_f32.as_() / area;
    let l0 = vec3::squared_norm(&v0);
    let l1 = vec3::squared_norm(&v1);
    let l2 = vec3::squared_norm(&v2);
    [
        (l1 + l2 - l0) * tmp,
        (l2 + l0 - l1) * tmp,
        (l0 + l1 - l2) * tmp
    ]
}

// ----------

/// Möller–Trumbore ray-triangle intersection algorithm
/// https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
/// https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection
pub fn ray_triangle_intersection(
    ray_org: &[f32],
    ray_dir: &[f32],
    p0: &[f32],
    p1: &[f32],
    p2: &[f32]) -> Option<f32> {
    use crate::vec3;
    // const EPSILON: f32 = 1.0e-4;
    let edge1 = vec3::sub(p1,p0);
    let edge2 = vec3::sub(p2, p0);
    let pvec = vec3::cross(&ray_dir, &edge2);
    let det = vec3::dot(&edge1, &pvec);
    // if det > -EPSILON && det < EPSILON { return None; }
    let invdet = 1.0 / det;
    let tvec = vec3::sub(ray_org, p0);
    let u = invdet * vec3::dot(&tvec, &pvec);
    if u < 0.0 || u > 1.0 { return None; }
    let qvec = vec3::cross(&tvec, &edge1);
    let v = invdet * vec3::dot(ray_dir, &qvec);
    if v < 0.0 || u + v > 1.0 { return None; }
    // At this stage we can compute t to find out where the intersection point is on the line.
    let t = invdet * vec3::dot(&edge2, &qvec);
    return Some(t);
    // ray intersection
    // if t > EPSILON { return Some(t); }
    // None
}