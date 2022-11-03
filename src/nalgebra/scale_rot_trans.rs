pub struct ScaleRotTrans {
    pub s: f32,
    pub quaternion: nalgebra::UnitQuaternion<f32>,
    pub translation: nalgebra::geometry::Translation3<f32>,
}

impl ScaleRotTrans {
    pub fn new() -> Self {
        ScaleRotTrans {
            s: 1_f32,
            quaternion: nalgebra::UnitQuaternion::<f32>::identity(),
            translation: nalgebra::geometry::Translation3::<f32>::new(0., 0., 0.),
        }
    }

    pub fn to_homogenous(&self) -> nalgebra::Matrix4<f32> {
        let ms = nalgebra::geometry::Scale3::new(
            self.s, self.s, self.s).to_homogeneous();
        let mt = nalgebra::geometry::Translation3::new(
            self.translation.x, self.translation.y, self.translation.z).to_homogeneous();
        let mr = self.quaternion.to_homogeneous();
        mt * mr * ms
    }
}