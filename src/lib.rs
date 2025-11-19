pub mod marching_cubes;
pub mod region;
pub mod solid;
pub mod stl;
// 「領域」の式ツリー、AIに読みやすいように定義
pub type Vec2 = glam::DVec2;
pub type Vec3 = glam::DVec3;
pub type Float = f64;
pub type Triangle = [Vec3; 3];