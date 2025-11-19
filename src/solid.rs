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
		offset: Vec3, //InAxesがあっても残した方が良い
	},
	Rotate {
		solid: Box<Solid>,
		axis: Axis,
		radian: Float,
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
			// ------------------
			// Primitive: Box
			// ------------------
			Solid::Box { half_extent } => {
				// iq の標準 box SDF
				let q = Vec3::new(p.x.abs(), p.y.abs(), p.z.abs()) - *half_extent;
				let outside = Vec3::new(q.x.max(0.0), q.y.max(0.0), q.z.max(0.0)).length();
				let inside = q.x.max(q.y).max(q.z).min(0.0);
				outside + inside
			}

			// ------------------
			// Primitive: Sphere
			// ------------------
			Solid::Sphere { center, radius } => (p - *center).length() - *radius,

			// ------------------
			// CSG
			// ------------------
			Solid::Union(a, b) => {
				let da = a.signed_distance(p);
				let db = b.signed_distance(p);
				da.min(db)
			}
			Solid::Intersect(a, b) => {
				let da = a.signed_distance(p);
				let db = b.signed_distance(p);
				da.max(db)
			}
			Solid::Subtract(a, b) => {
				let da = a.signed_distance(p);
				let db = b.signed_distance(p);
				da.max(-db)
			}

			// ------------------
			// Translate（Move）
			// ------------------
			Solid::Translate { solid, offset } => solid.signed_distance(p - *offset),

			// ------------------
			// Rotate
			// ------------------
			Solid::Rotate {
				solid,
				axis,
				radian,
			} => {
				// 「solid を radian だけ回転させた」と解釈するなら
				// SDF的には点 p を「逆回転」させてから評価する
				let angle = -*radian; // ← ここを +*radian にすると回転の向きが逆になる

				let c = angle.cos();
				let s = angle.sin();

				let q = match axis {
					Axis::X => {
						// X軸まわり回転
						Vec3::new(p.x, c * p.y - s * p.z, s * p.y + c * p.z)
					}
					Axis::Y => {
						// Y軸まわり回転
						Vec3::new(c * p.x + s * p.z, p.y, -s * p.x + c * p.z)
					}
					Axis::Z => {
						// Z軸まわり回転
						Vec3::new(c * p.x - s * p.y, s * p.x + c * p.y, p.z)
					}
				};

				// Box<Solid> なので &Box<Solid> -> &Solid へ自動 Deref されてそのまま呼べる
				solid.signed_distance(q)
			}

			// ------------------
			// Scale（Uniform）
			// ------------------
			Solid::Scale { solid, factor } => {
				let inv = 1.0 / *factor;
				solid.signed_distance(p * inv) * *factor
			}

			// ------------------
			// InAxes（ローカル座標系で評価）
			// ------------------
			Solid::InAxes {
				solid,
				origin,
				z_axis,
				x_hint,
			} => {
				// ローカル軸を作る
				let z = z_axis.normalize();

				// x_hint を z と直交化（Gram–Schmidt）
				let mut x = *x_hint - z * x_hint.dot(z);
				if x.length_squared() < 1e-12 {
					// フォールバック方向
					let mut h = Vec3::new(1.0, 0.0, 0.0);
					if z.dot(h).abs() > 0.9 {
						h = Vec3::new(0.0, 1.0, 0.0);
					}
					x = h - z * h.dot(z);
				}
				x = x.normalize();

				// 右手系を構築
				let y = z.cross(x);

				// 世界座標 → ローカル座標
				let v = p - *origin;
				let local = Vec3::new(v.dot(x), v.dot(y), v.dot(z));

				solid.signed_distance(local)
			}

			// ------------------
			// Extrude (Region × axis-range)
			// ------------------
			Solid::Extrude {
				region,
				axis,
				range,
			} => {
				// normalize range
				let (min, max) = if range[0] <= range[1] {
					(range[0], range[1])
				} else {
					(range[1], range[0])
				};

				// p → (u, v, w)
				let (u, v, w) = match axis {
					Axis::Z => (p.x, p.y, p.z),
					Axis::X => (p.y, p.z, p.x), // YZ面 → X軸に押し出し
					Axis::Y => (p.z, p.x, p.y), // ZX面 → Y軸に押し出し
				};

				// 2D SDF
				let d2 = region.signed_distance(Vec2::new(u, v));

				// w のスラブ領域 SDF
				let c = 0.5 * (min + max);
				let h = 0.5 * (max - min);
				let dw = (w - c).abs() - h;

				// 押し出し SDF（有限高さ）
				d2.max(dw)
			}
		}
	}
	/// 軸平行バウンディングボックス [min, max]
	pub fn bounding_box(&self) -> [Vec3; 2] {
		fn vmin(a: Vec3, b: Vec3) -> Vec3 {
			Vec3::new(a.x.min(b.x), a.y.min(b.y), a.z.min(b.z))
		}
		fn vmax(a: Vec3, b: Vec3) -> Vec3 {
			Vec3::new(a.x.max(b.x), a.y.max(b.y), a.z.max(b.z))
		}

		match self {
			// ------------------
			// Primitive: Box
			// ------------------
			Solid::Box { half_extent } => {
				let min = -(*half_extent);
				let max = *half_extent;
				[min, max]
			}

			// ------------------
			// Primitive: Sphere
			// ------------------
			Solid::Sphere { center, radius } => {
				let r = Vec3::splat(*radius);
				[*center - r, *center + r]
			}

			// ------------------
			// CSG
			// ------------------
			Solid::Union(a, b) => {
				let [amin, amax] = a.bounding_box();
				let [bmin, bmax] = b.bounding_box();
				[vmin(amin, bmin), vmax(amax, bmax)]
			}
			Solid::Intersect(a, b) => {
				// 幾何的には「AABB同士の交差」を取る
				let [amin, amax] = a.bounding_box();
				let [bmin, bmax] = b.bounding_box();
				let min = vmax(amin, bmin);
				let max = vmin(amax, bmax);
				[min, max]
			}
			Solid::Subtract(a, _b) => {
				// A \ B のバウンディングは A の AABB 以上には広がらない
				a.bounding_box()
			}

			// ------------------
			// Translate
			// ------------------
			Solid::Translate { solid, offset } => {
				let [min, max] = solid.bounding_box();
				[min + *offset, max + *offset]
			}

			// ------------------
			// Rotate
			// ------------------
			Solid::Rotate {
				solid,
				axis,
				radian,
			} => {
				// まず「回転前」の AABB を取る
				let [bmin, bmax] = solid.bounding_box();

				// 物体を +radian だけ回転させたと解釈するので、AABB の頂点を +radian で回す
				let angle = *radian;
				let c = angle.cos();
				let s = angle.sin();

				// 軸ごとの回転
				fn rotate_vec(p: Vec3, axis: &Axis, c: Float, s: Float) -> Vec3 {
					match axis {
						Axis::X => Vec3::new(p.x, c * p.y - s * p.z, s * p.y + c * p.z),
						Axis::Y => Vec3::new(c * p.x + s * p.z, p.y, -s * p.x + c * p.z),
						Axis::Z => Vec3::new(c * p.x - s * p.y, s * p.x + c * p.y, p.z),
					}
				}

				let mut wmin = Vec3::new(Float::INFINITY, Float::INFINITY, Float::INFINITY);
				let mut wmax = Vec3::new(
					Float::NEG_INFINITY,
					Float::NEG_INFINITY,
					Float::NEG_INFINITY,
				);

				// 8頂点を全部回転して、その AABB を取る
				for &x in &[bmin.x, bmax.x] {
					for &y in &[bmin.y, bmax.y] {
						for &z in &[bmin.z, bmax.z] {
							let local = Vec3::new(x, y, z);
							let world = rotate_vec(local, axis, c, s);
							wmin = vmin(wmin, world);
							wmax = vmax(wmax, world);
						}
					}
				}

				[wmin, wmax]
			}

			// ------------------
			// Scale（Uniform）
			// ------------------
			Solid::Scale { solid, factor } => {
				let [min, max] = solid.bounding_box();
				// factor > 0 を想定。負の場合も一応対応
				if *factor >= 0.0 {
					[min * *factor, max * *factor]
				} else {
					let f = *factor;
					[max * f, min * f]
				}
			}

			// ------------------
			// InAxes（任意姿勢のローカル形状）
			// ------------------
			Solid::InAxes {
				solid,
				origin,
				z_axis,
				x_hint,
			} => {
				// まずローカル空間での AABB を取得
				let [lmin, lmax] = solid.bounding_box();

				// ローカル座標系を構築（signed_distance と同じ）
				let z = z_axis.normalize();
				let mut x = *x_hint - z * x_hint.dot(z);
				if x.length_squared() < 1e-12 {
					// フォールバック
					let mut h = Vec3::new(1.0, 0.0, 0.0);
					if z.dot(h).abs() > 0.9 {
						h = Vec3::new(0.0, 1.0, 0.0);
					}
					x = h - z * h.dot(z);
				}
				x = x.normalize();
				let y = z.cross(x);

				// ローカルの8頂点を世界座標に変換して AABB を取る
				let mut wmin = Vec3::new(Float::INFINITY, Float::INFINITY, Float::INFINITY);
				let mut wmax = Vec3::new(
					Float::NEG_INFINITY,
					Float::NEG_INFINITY,
					Float::NEG_INFINITY,
				);

				for &lx in &[lmin.x, lmax.x] {
					for &ly in &[lmin.y, lmax.y] {
						for &lz in &[lmin.z, lmax.z] {
							let local = Vec3::new(lx, ly, lz);
							let world = *origin + x * local.x + y * local.y + z * local.z;
							wmin = vmin(wmin, world);
							wmax = vmax(wmax, world);
						}
					}
				}

				[wmin, wmax]
			}

			// ------------------
			// Extrude (Region × axis-range)
			// ------------------
			Solid::Extrude {
				region,
				axis,
				range,
			} => {
				// range を正規化
				let (minr, maxr) = if range[0] <= range[1] {
					(range[0], range[1])
				} else {
					(range[1], range[0])
				};

				// 2D Region の AABB
				let [rmin, rmax] = region.bounding_box(); // [Vec2; 2]

				match axis {
					Axis::Z => {
						let min = Vec3::new(rmin.x, rmin.y, minr);
						let max = Vec3::new(rmax.x, rmax.y, maxr);
						[min, max]
					}
					Axis::X => {
						// YZ 面に region, X に range
						let min = Vec3::new(minr, rmin.x, rmin.y);
						let max = Vec3::new(maxr, rmax.x, rmax.y);
						[min, max]
					}
					Axis::Y => {
						// ZX 面に region, Y に range
						// region.x -> Z, region.y -> X という対応（signed_distance に合わせる）
						let min = Vec3::new(rmin.y, minr, rmin.x);
						let max = Vec3::new(rmax.y, maxr, rmax.x);
						[min, max]
					}
				}
			}
		}
	}
}
