//! implementation of the "Adaptively Sampled Distance Field"
//!
//!　```Sarah F. Frisken, Ronald N. Perry, Alyn P. Rockwood, and Thouis R. Jones. 2000.
//! Adaptively sampled distance fields: a general representation of shape for computer graphics.
//! In Proceedings of the 27th annual conference on Computer graphics and interactive techniques (SIGGRAPH '00).```

#[derive(Clone)]
struct Node {
    // center
    c: [f64; 3],
    // radius
    r: f64,
    // index of children
    i: [usize; 8],
    // distance at the corners
    d: [f64; 8]
}

fn corner_distances(
    sdf: fn([f64; 3]) -> f64,
    c: [f64; 3],
    r: f64) -> [f64; 8] {
    [
        sdf([c[0] - r, c[1] - r, c[2] - r]),
        sdf([c[0] + r, c[1] - r, c[2] - r]),
        sdf([c[0] + r, c[1] + r, c[2] - r]),
        sdf([c[0] - r, c[1] + r, c[2] - r]),
        sdf([c[0] - r, c[1] - r, c[2] + r]),
        sdf([c[0] + r, c[1] - r, c[2] + r]),
        sdf([c[0] + r, c[1] + r, c[2] + r]),
        sdf([c[0] - r, c[1] + r, c[2] + r])
    ]
}

fn make_child_tree(
    nodes: &mut Vec<Node>,
    sdf: fn([f64; 3]) -> f64,
    c0: [f64; 3],
    r0: f64,
    d0: [f64; 8],
    min_r: f64,
    max_r: f64) -> Node
{
    if r0 * 0.5 < min_r {
        return Node { c: [0.; 3], d: [0.; 8], r: 0.0, i: [usize::MAX; 8] };
    }
    // Edges
    let v100 = sdf([c0[0], c0[1] - r0, c0[2] - r0]);
    let v210 = sdf([c0[0] + r0, c0[1], c0[2] - r0]);
    let v120 = sdf([c0[0], c0[1] + r0, c0[2] - r0]);
    let v010 = sdf([c0[0] - r0, c0[1], c0[2] - r0]);

    let v001 = sdf([c0[0] - r0, c0[1] - r0, c0[2]]);
    let v201 = sdf([c0[0] + r0, c0[1] - r0, c0[2]]);
    let v221 = sdf([c0[0] + r0, c0[1] + r0, c0[2]]);
    let v021 = sdf([c0[0] - r0, c0[1] + r0, c0[2]]);

    let v102 = sdf([c0[0], c0[1] - r0, c0[2] + r0]);
    let v212 = sdf([c0[0] + r0, c0[1], c0[2] + r0]);
    let v122 = sdf([c0[0], c0[1] + r0, c0[2] + r0]);
    let v012 = sdf([c0[0] - r0, c0[1], c0[2] + r0]);

    // Faces
    let v101 = sdf([c0[0], c0[1] - r0, c0[2]]);
    let v211 = sdf([c0[0] + r0, c0[1], c0[2]]);
    let v121 = sdf([c0[0], c0[1] + r0, c0[2]]);
    let v011 = sdf([c0[0] - r0, c0[1], c0[2]]);
    let v110 = sdf([c0[0], c0[1], c0[2] - r0]);
    let v112 = sdf([c0[0], c0[1], c0[2] + r0]);

    // Center
    let v111 = sdf([c0[0], c0[1], c0[2]]);

    let mut is_subdivide: Option<bool> = None;
    {
        if r0 * 0.5 > max_r { is_subdivide = Some(true) }
        if is_subdivide.is_none() {
            let min_dist = [
                v100, v210, v120, v010,
                v001, v201, v221, v021,
                v102, v212, v122, v012,
                v101, v211, v121, v011, v110, v112,
                v111].iter().fold(0. / 0., |m, v| v.min(m));
            if min_dist > r0 * 1.8 { // there is no object inside. (1.8 is larger than the diagonal)
                is_subdivide = Some(false);
            }
            if min_dist < min_r {
                is_subdivide = Some(true);
            }
        }
        if is_subdivide.is_none() { // edges
            let w100 = (d0[0] + d0[1]) * 0.5;
            let w210 = (d0[1] + d0[2]) * 0.5;
            let w120 = (d0[2] + d0[3]) * 0.5;
            let w010 = (d0[3] + d0[0]) * 0.5;
            let w102 = (d0[4] + d0[5]) * 0.5;
            let w212 = (d0[5] + d0[6]) * 0.5;
            let w122 = (d0[6] + d0[7]) * 0.5;
            let w012 = (d0[7] + d0[4]) * 0.5;
            let w001 = (d0[0] + d0[4]) * 0.5;
            let w201 = (d0[1] + d0[5]) * 0.5;
            let w221 = (d0[2] + d0[6]) * 0.5;
            let w021 = (d0[3] + d0[7]) * 0.5;
            let t = min_r * 0.8;
            if (v100 - w100).abs() > t { is_subdivide = Some(true); }
            if (v210 - w210).abs() > t { is_subdivide = Some(true); }
            if (v120 - w120).abs() > t { is_subdivide = Some(true); }
            if (v010 - w010).abs() > t { is_subdivide = Some(true); }
            if (v102 - w102).abs() > t { is_subdivide = Some(true); }
            if (v212 - w212).abs() > t { is_subdivide = Some(true); }
            if (v122 - w122).abs() > t { is_subdivide = Some(true); }
            if (v012 - w012).abs() > t { is_subdivide = Some(true); }
            if (v001 - w001).abs() > t { is_subdivide = Some(true); }
            if (v201 - w201).abs() > t { is_subdivide = Some(true); }
            if (v221 - w221).abs() > t { is_subdivide = Some(true); }
            if (v021 - w021).abs() > t { is_subdivide = Some(true); }
        }
        if is_subdivide.is_none() { // faces, center
            let w101 = (d0[0] + d0[1] + d0[4] + d0[5]) * 0.25;
            let w211 = (d0[1] + d0[2] + d0[5] + d0[6]) * 0.25;
            let w121 = (d0[2] + d0[3] + d0[6] + d0[7]) * 0.25;
            let w011 = (d0[3] + d0[0] + d0[7] + d0[4]) * 0.25;
            let w110 = (d0[0] + d0[1] + d0[2] + d0[3]) * 0.25;
            let w112 = (d0[4] + d0[5] + d0[6] + d0[7]) * 0.25;
            let w111 = (d0[0] + d0[1] + d0[2] + d0[3] + d0[4] + d0[5] + d0[6] + d0[7]) * 0.125;
            let t = min_r * 0.8;
            if (v101 - w101).abs() > t { is_subdivide = Some(true); }
            if (v211 - w211).abs() > t { is_subdivide = Some(true); }
            if (v121 - w121).abs() > t { is_subdivide = Some(true); }
            if (v011 - w011).abs() > t { is_subdivide = Some(true); }
            if (v110 - w110).abs() > t { is_subdivide = Some(true); }
            if (v112 - w112).abs() > t { is_subdivide = Some(true); }
            if (v111 - w111).abs() > t { is_subdivide = Some(true); }
        }
    }

    if is_subdivide != Some(true) {
        return Node { c: [0.; 3], d: [0.; 8], r: 0.0, i: [usize::MAX; 8] };
    }

    let nchild0 = nodes.len();
    nodes.resize(nodes.len() + 8, Node { c: [0.; 3], d: [0.; 8], r: 0.0, i: [usize::MAX; 8] });
    // left-bottom
    nodes[nchild0 + 0] = make_child_tree(
        nodes, sdf,
        [c0[0] - r0 * 0.5, c0[1] - r0 * 0.5, c0[2] - r0 * 0.5],
        r0 * 0.5,
        [d0[0], v100, v110, v010, v001, v101, v111, v011],
        min_r, max_r);
    nodes[nchild0 + 1] = make_child_tree(
        nodes, sdf,
        [c0[0] + r0 * 0.5, c0[1] - r0 * 0.5, c0[2] - r0 * 0.5],
        r0 * 0.5,
        [v100, d0[1], v210, v110, v101, v201, v211, v111],
        min_r, max_r);
    nodes[nchild0 + 2] = make_child_tree(
        nodes, sdf,
        [c0[0] + r0 * 0.5, c0[1] + r0 * 0.5, c0[2] - r0 * 0.5],
        r0 * 0.5,
        [v110, v210, d0[2], v120, v111, v211, v221, v121],
        min_r, max_r);
    nodes[nchild0 + 3] = make_child_tree(
        nodes, sdf,
        [c0[0] - r0 * 0.5, c0[1] + r0 * 0.5, c0[2] - r0 * 0.5],
        r0 * 0.5,
        [v010, v110, v120, d0[3], v011, v111, v121, v021],
        min_r, max_r);
    nodes[nchild0 + 4] = make_child_tree(
        nodes, sdf,
        [c0[0] - r0 * 0.5, c0[1] - r0 * 0.5, c0[2] + r0 * 0.5],
        r0 * 0.5,
        [v001, v101, v111, v011, d0[4], v102, v112, v012],
        min_r, max_r);
    nodes[nchild0 + 5] = make_child_tree(
        nodes, sdf,
        [c0[0] + r0 * 0.5, c0[1] - r0 * 0.5, c0[2] + r0 * 0.5],
        r0 * 0.5,
        [v101, v201, v211, v111, v102, d0[5], v212, v112],
        min_r, max_r);
    nodes[nchild0 + 6] = make_child_tree(
        nodes, sdf,
        [c0[0] + r0 * 0.5, c0[1] + r0 * 0.5, c0[2] + r0 * 0.5],
        r0 * 0.5,
        [v111, v211, v221, v121, v112, v212, d0[6], v122, ],
        min_r, max_r);
    nodes[nchild0 + 7] = make_child_tree(
        nodes, sdf,
        [c0[0] - r0 * 0.5, c0[1] + r0 * 0.5, c0[2] + r0 * 0.5],
        r0 * 0.5,
        [v011, v111, v121, v021, v012, v112, v122, d0[7]],
        min_r, max_r);
    return Node {
        c: c0,
        d: d0,
        r: r0,
        i: [nchild0 + 0, nchild0 + 1, nchild0 + 2, nchild0 + 3,
            nchild0 + 4, nchild0 + 5, nchild0 + 6, nchild0 + 7],
    };
}


pub struct AdaptiveDistanceField {
    nodes: Vec<Node>,
}

impl AdaptiveDistanceField {
    pub fn new() -> Self {
        AdaptiveDistanceField { nodes: vec!() }
    }
    pub fn build(&mut self, sdf: fn([f64; 3]) -> f64, center: [f64; 3], radius: f64) {
        self.nodes.reserve(1024 * 64);
        self.nodes.resize(1, Node { c: [0.;3], r: 0., d : [0.;8], i:[usize::MAX;8] });
        self.nodes[0] = make_child_tree(
            &mut self.nodes,
            sdf,
            center, radius, corner_distances(sdf, center, radius),
            radius * (0.99 / 64.0), radius * (1.01 / 4.0));
    }
    pub fn edge_to_xyz(&self) -> Vec<f32> {
        let mut edge2xyz = vec!();
        for v in self.nodes.iter() {
            let ps = [
                [v.c[0]-v.r,v.c[1]-v.r, v.c[2]-v.r],
                [v.c[0]-v.r,v.c[1]-v.r, v.c[2]+v.r],
                [v.c[0]-v.r,v.c[1]+v.r, v.c[2]-v.r],
                [v.c[0]-v.r,v.c[1]+v.r, v.c[2]+v.r],
                [v.c[0]+v.r,v.c[1]-v.r, v.c[2]-v.r],
                [v.c[0]+v.r,v.c[1]-v.r, v.c[2]+v.r],
                [v.c[0]+v.r,v.c[1]+v.r, v.c[2]-v.r],
                [v.c[0]+v.r,v.c[1]+v.r, v.c[2]+v.r] ];
            let es = [
                [0,1], [2,3], [4,5], [6,7],
                [0,2], [2,6], [6,4], [4,0],
                [1,3], [3,7], [7,5], [5,1]
            ];
            for e in es.iter() {
                let i0 = e[0];
                let i1 = e[1];
                edge2xyz.push(ps[i0][0] as f32);
                edge2xyz.push(ps[i0][1] as f32);
                edge2xyz.push(ps[i0][2] as f32);
                edge2xyz.push(ps[i1][0] as f32);
                edge2xyz.push(ps[i1][1] as f32);
                edge2xyz.push(ps[i1][2] as f32);
            }
        }
        edge2xyz
    }
}