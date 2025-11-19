use crate::region::Region;
use crate::{Float, Vec2, Vec3};

#[derive(Clone, Copy, Debug)]
pub enum Axis {
	X,
	Y,
	Z,
}

#[derive(Clone, Debug)]
pub enum Solid {
	// プリミティブ
	Box {
		half_extent: Vec3, // 原点中心の直方体: [-h.x, h.x] × [-h.y, h.y] × [-h.z, h.z]
	},
	Sphere {
		center: Vec3,
		radius: Float,
	},
	// CSG
	Union(Box<Solid>, Box<Solid>),
	Intersect(Box<Solid>, Box<Solid>),
	Subtract(Box<Solid>, Box<Solid>),
	// 変換
	Translate {
		solid: Box<Solid>,
		offset: Vec3,//InAxesがあっても残した方が良い
	},
	Scale {
		solid: Box<Solid>,
		factor: Float,
	},
	InAxes {
		solid: Box<Solid>,
		origin: Vec3,
		z_axis: Vec3, // 法線（Normal）。最も重要な「面の向き」
		x_hint: Vec3, // 面上のX方向（Zと直交化される）
	},
	// 2D Region の押し出し
	Extrude {
		region: Box<Region>, // xy平面上の 2D Region
		axis: Axis,
		range: [Float; 2],
	},
}

impl Solid {
	/// 3D SDF: inside < 0, outside > 0
	pub fn signed_distance(&self, p: Vec3) -> Float {
		match self {
			_=>todo!()
		}
	}
}
