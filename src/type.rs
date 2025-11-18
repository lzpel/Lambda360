// 「領域」の式ツリー、AIに読みやすいように定義
pub type Vec2 = glam::DVec2;
pub type Float = f64;
#[derive(Clone, Debug)]
pub enum Region {
	// プリミティブ
	Circle { center: Vec2, radius: Float }, //例えばCircleやRectはTranslateとScaleを用いれば変数なしにできるけどそこまでプリミティブにしなくていいという判断？
	Rect { min: Vec2, max: Vec2 },
	Polygon { vertices: Vec<Vec2> },

	// CSG
	Union(Box<Self>, Box<Self>),
	Intersect(Box<Self>, Box<Self>),
	Subtract(Box<Self>, Box<Self>),
	Minkowski(Box<Self>, Box<Self>),

	// 変換
	Translate { region: Box<Self>, offset: Vec2 },
	Scale { region: Box<Self>, factor: Float },
	Rotate { region: Box<Self>, radian: Float },
}

impl Region {
	fn contains(&self, p: Vec2) -> bool {
		// match self { ... } で再帰的に評価
		match self {
			// 中心との距離が半径以下なら内側
			Region::Circle { center, radius } => {
				let d2 = (p - center).length_squared();
				d2 <= radius * radius
			}
			// 軸平行矩形：min/max の間にあれば内側
			Region::Rect { min, max } => {
				p.x >= min.x && p.y >= min.y && p.x <= max.x && p.y <= max.y
			}
			// 多角形：典型的な「レイキャスト法」で内外判定
			Region::Polygon { vertices } => {
				let n = vertices.len();
				if n < 3 {
					return false;
				}

				let mut inside = false;
				let px = p.x;
				let py = p.y;

				// (i, j) 辺を順番に見ていく
				for i in 0..n {
					let vi = vertices[i];
					let vj = vertices[if i == 0 { n - 1 } else { i - 1 }]; //前

					// 辺 (vj -> vi) が p.y をまたいでいるか？
					let intersect = ((vi.y > py) != (vj.y > py))
						&& (px < (vj.x - vi.x) * (py - vi.y) / (vj.y - vi.y) + vi.x);

					if intersect {
						inside = !inside;
					}
				}

				inside
			}
			Region::Union(a, b) => a.contains(p) || b.contains(p),
			Region::Intersect(a, b) => a.contains(p) && b.contains(p),
			Region::Subtract(a, b) => a.contains(p) && !b.contains(p),
			Region::Minkowski(_a, _b) => {
				// 本気でやるなら「∃x s.t. x∈A ∧ p-x∈B」を探す必要があるので、
				// contains で厳密実装するのはかなり重い・難しい。
				// ここでは未実装にしておく。
				unimplemented!("contains() for Minkowski is not implemented")
			}
			Region::Translate { region, offset } => {
				// 逆平行移動 (p - offset)
				region.contains(p - *offset)
			}
			Region::Scale { region, factor } => {
				if *factor == 0.0 {
					// scale=0 → 1 点に潰れるので、代表点として原点を使う
					region.contains(Vec2::ZERO)
				} else {
					// 逆スケール (p / factor)
					region.contains(p / *factor)
				}
			}
			Region::Rotate { region, radian } => {
				// 逆回転 (回転角 -θ)
				let c = radian.cos();
				let s = radian.sin();
				// R(-θ) * p
				let local = Vec2::new(c * p.x + s * p.y, -s * p.x + c * p.y);
				region.contains(local)
			}
		}
	}

	fn bounding_box(&self) -> [Vec2; 2] {
		match self {
			// 円：中心 ± (r, r)
			Region::Circle { center, radius } => {
				let r = *radius;
				let min = Vec2::new(center.x - r, center.y - r);
				let max = Vec2::new(center.x + r, center.y + r);
				[min, max]
			}

			// Rect はそのまま
			Region::Rect { min, max } => [*min, *max],

			// 多角形：全頂点の min/max
			Region::Polygon { vertices } => {
				assert!(!vertices.is_empty(), "Polygon with no vertices");
				let mut min = vertices[0];
				let mut max = vertices[0];
				for v in vertices.iter().skip(1) {
					min.x = min.x.min(v.x);
					min.y = min.y.min(v.y);
					max.x = max.x.max(v.x);
					max.y = max.y.max(v.y);
				}
				[min, max]
			}

			// CSG
			Region::Union(a, b) => {
				let [amin, amax] = a.bounding_box();
				let [bmin, bmax] = b.bounding_box();
				let min = Vec2::new(amin.x.min(bmin.x), amin.y.min(bmin.y));
				let max = Vec2::new(amax.x.max(bmax.x), amax.y.max(bmax.y));
				[min, max]
			}

			Region::Intersect(a, b) => {
				// 交差の AABB は両者 AABB の交差
				let [amin, amax] = a.bounding_box();
				let [bmin, bmax] = b.bounding_box();
				let min = Vec2::new(amin.x.max(bmin.x), amin.y.max(bmin.y));
				let max = Vec2::new(amax.x.min(bmax.x), amax.y.min(bmax.y));
				// min.x > max.x 等の場合は「空」だけど、そのまま返す（扱いは呼び出し側に任せる）
				[min, max]
			}

			Region::Subtract(a, _b) => {
				// A \ B の AABB は A の AABB 以上には広がらないので、
				// 過小でない近似として A の AABB をそのまま返す
				a.bounding_box()
			}

			Region::Minkowski(a, b) => {
				// Minkowski 和の AABB は、両者の AABB の和で OK
				// AABB(A ⊕ B) = [amin + bmin, amax + bmax]
				let [amin, amax] = a.bounding_box();
				let [bmin, bmax] = b.bounding_box();
				let min = amin + bmin;
				let max = amax + bmax;
				[min, max]
			}

			// 変換
			Region::Translate { region, offset } => {
				let [min, max] = region.bounding_box();
				[min + *offset, max + *offset]
			}

			Region::Scale { region, factor } => {
				let [min, max] = region.bounding_box();
				if *factor == 0.0 {
					// 原点に潰れる
					[Vec2::ZERO, Vec2::ZERO]
				} else {
					let s = *factor;
					let smin = min * s;
					let smax = max * s;
					// s が負のとき、min/max が反転しうるので取り直す
					let new_min = Vec2::new(smin.x.min(smax.x), smin.y.min(smax.y));
					let new_max = Vec2::new(smin.x.max(smax.x), smin.y.max(smax.y));
					[new_min, new_max]
				}
			}

			Region::Rotate { region, radian } => {
				// 元の AABB の 4頂点を θ 回転させて、その AABB を取る
				let [min, max] = region.bounding_box();
				let c = radian.cos();
				let s = radian.sin();

				let corners = [
					Vec2::new(min.x, min.y),
					Vec2::new(min.x, max.y),
					Vec2::new(max.x, min.y),
					Vec2::new(max.x, max.y),
				];

				let mut rmin = Vec2::new(f64::INFINITY, f64::INFINITY);
				let mut rmax = Vec2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);

				for v in &corners {
					// R(θ) * v
					let x = c * v.x - s * v.y;
					let y = s * v.x + c * v.y;
					rmin.x = rmin.x.min(x);
					rmin.y = rmin.y.min(y);
					rmax.x = rmax.x.max(x);
					rmax.y = rmax.y.max(y);
				}

				[rmin, rmax]
			}
		}
	}

	fn signed_distance(&self, p: Vec2) -> Float {
		match self {
			// 円: 距離-半径
			Region::Circle { center, radius } => (p - *center).length() - *radius,

			// 軸平行矩形: AABB の SDF
			Region::Rect { min, max } => {
				// IQ: box SDF
				let c = (*min + *max) * 0.5;
				let r = (*max - *min) * 0.5;
				let q = (p - c).abs() - r;

				let outside = Vec2::new(q.x.max(0.0), q.y.max(0.0)).length();
				let inside = q.x.max(q.y).min(0.0);
				outside + inside
			}

			// 多角形: エッジへの最近距離 + 内外で符号付け
			Region::Polygon { vertices } => {
				let n = vertices.len();
				assert!(n >= 3, "Polygon with fewer than 3 vertices");

				// 内外判定（contains と同じレイキャスト）
				let mut inside = false;
				let px = p.x;
				let py = p.y;

				for i in 0..n {
					let vi = vertices[i];
					let vj = vertices[if i == 0 { n - 1 } else { i - 1 }];

					let intersect = ((vi.y > py) != (vj.y > py))
						&& (px < (vj.x - vi.x) * (py - vi.y) / (vj.y - vi.y) + vi.x);

					if intersect {
						inside = !inside;
					}
				}

				// エッジへの最短距離
				let mut min_dist2 = f64::INFINITY;
				for i in 0..n {
					let v0 = vertices[i];
					let v1 = vertices[(i + 1) % n];
					let e = v1 - v0;
					let w = p - v0;
					let t = w.dot(e) / e.length_squared();
					let t_clamped = t.clamp(0.0, 1.0);
					let proj = v0 + e * t_clamped;
					let d2 = (p - proj).length_squared();
					if d2 < min_dist2 {
						min_dist2 = d2;
					}
				}

				let dist = min_dist2.sqrt();
				if inside { -dist } else { dist }
			}

			// CSG
			Region::Union(a, b) => {
				let da = a.signed_distance(p);
				let db = b.signed_distance(p);
				da.min(db)
			}

			Region::Intersect(a, b) => {
				let da = a.signed_distance(p);
				let db = b.signed_distance(p);
				da.max(db)
			}

			Region::Subtract(a, b) => {
				let da = a.signed_distance(p);
				let db = b.signed_distance(p);
				da.max(-db)
			}

			Region::Minkowski(_a, _b) => {
				// 厳密な Minkowski 和の SDF は
				//   d_{A⊕B}(p) = inf_x { d_A(x) + d_B(p - x) }
				// 的な形になり、汎用にはかなり重いのでここでは未実装にしておく。
				unimplemented!("signed_distance() for Minkowski is not implemented")
			}

			// 変換
			Region::Translate { region, offset } => {
				// 逆平行移動
				region.signed_distance(p - *offset)
			}

			Region::Scale { region, factor } => {
				if *factor == 0.0 {
					// 原点に潰れるので、とりあえず原点での距離を返す（近似）
					region.signed_distance(Vec2::ZERO)
				} else {
					// 一様スケール: d'(p) = d(p/s) * |s|
					let s = *factor;
					region.signed_distance(p / s) * s.abs()
				}
			}

			Region::Rotate { region, radian } => {
				// 逆回転してから元の距離関数を評価
				let c = radian.cos();
				let s = radian.sin();
				let local = Vec2::new(c * p.x + s * p.y, -s * p.x + c * p.y);
				region.signed_distance(local)
			}
		}
	}
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn contains_circle() {
		let c = Region::Circle {
			center: Vec2::new(0.0, 0.0),
			radius: 1.0,
		};

		assert!(c.contains(Vec2::new(0.0, 0.0)));
		assert!(c.contains(Vec2::new(0.5, 0.5)));
		assert!(!c.contains(Vec2::new(2.0, 0.0)));
	}

	#[test]
	fn contains_rect() {
		let r = Region::Rect {
			min: Vec2::new(-1.0, -1.0),
			max: Vec2::new(1.0, 1.0),
		};

		assert!(r.contains(Vec2::new(0.0, 0.0)));
		assert!(r.contains(Vec2::new(1.0, 1.0)));
		assert!(!r.contains(Vec2::new(2.0, 0.0)));
	}

	#[test]
	fn contains_polygon_triangle() {
		let poly = Region::Polygon {
			vertices: vec![
				Vec2::new(0.0, 0.0),
				Vec2::new(2.0, 0.0),
				Vec2::new(1.0, 2.0),
			],
		};

		assert!(poly.contains(Vec2::new(1.0, 0.5)));
		assert!(!poly.contains(Vec2::new(0.0, 3.0)));
	}

	#[test]
	fn contains_translate() {
		let base = Region::Circle {
			center: Vec2::new(0.0, 0.0),
			radius: 1.0,
		};
		let r = Region::Translate {
			region: Box::new(base),
			offset: Vec2::new(2.0, 0.0),
		};

		assert!(r.contains(Vec2::new(2.0, 0.0)));
		assert!(!r.contains(Vec2::new(0.0, 0.0)));
	}

	#[test]
	fn contains_scale() {
		let base = Region::Circle {
			center: Vec2::new(0.0, 0.0),
			radius: 1.0,
		};

		let r = Region::Scale {
			region: Box::new(base),
			factor: 2.0,
		};

		// 半径 2 の円になったとみなせる
		assert!(r.contains(Vec2::new(1.5, 0.0)));
		assert!(!r.contains(Vec2::new(3.0, 0.0)));
	}

	#[test]
	fn contains_rotate() {
		let rect = Region::Rect {
			min: Vec2::new(-1.0, -0.5),
			max: Vec2::new(1.0, 0.5),
		};
		let rotated = Region::Rotate {
			region: Box::new(rect),
			radian: std::f64::consts::FRAC_PI_2, // 90°
		};

		// (1,0) は元の rect の (0, -1) に対応 → outside
		assert!(!rotated.contains(Vec2::new(1.0, 0.0)));

		// (0,1) は元の rect の (1,0) に対応 → inside
		assert!(rotated.contains(Vec2::new(0.0, 1.0)));
	}

	// -------------------------
	// bounding_box tests
	// -------------------------
	#[test]
	fn bounding_box_circle() {
		let c = Region::Circle {
			center: Vec2::new(1.0, 2.0),
			radius: 1.5,
		};

		let [min, max] = c.bounding_box();
		assert_eq!(min, Vec2::new(-0.5, 0.5));
		assert_eq!(max, Vec2::new(2.5, 3.5));
	}

	#[test]
	fn bounding_box_union() {
		let a = Region::Rect {
			min: Vec2::new(0.0, 0.0),
			max: Vec2::new(1.0, 1.0),
		};
		let b = Region::Rect {
			min: Vec2::new(2.0, 2.0),
			max: Vec2::new(3.0, 3.0),
		};

		let u = Region::Union(Box::new(a), Box::new(b));
		let [min, max] = u.bounding_box();
		assert_eq!(min, Vec2::new(0.0, 0.0));
		assert_eq!(max, Vec2::new(3.0, 3.0));
	}

	#[test]
	fn bounding_box_translate() {
		let base = Region::Rect {
			min: Vec2::new(-1.0, -1.0),
			max: Vec2::new(1.0, 1.0),
		};

		let t = Region::Translate {
			region: Box::new(base),
			offset: Vec2::new(5.0, -2.0),
		};

		let [min, max] = t.bounding_box();
		assert_eq!(min, Vec2::new(4.0, -3.0));
		assert_eq!(max, Vec2::new(6.0, -1.0));
	}

	// -------------------------
	// signed_distance tests
	// -------------------------
	#[test]
	fn signed_distance_circle() {
		let c = Region::Circle {
			center: Vec2::new(0.0, 0.0),
			radius: 1.0,
		};

		assert!((c.signed_distance(Vec2::new(0.0, 0.0)) + 1.0).abs() < 1e-6);
		assert!((c.signed_distance(Vec2::new(2.0, 0.0)) - 1.0).abs() < 1e-6);
	}

	#[test]
	fn signed_distance_rect() {
		let r = Region::Rect {
			min: Vec2::new(-1.0, -1.0),
			max: Vec2::new(1.0, 1.0),
		};

		assert!(r.signed_distance(Vec2::new(0.0, 0.0)) < 0.0); // inside
		assert!((r.signed_distance(Vec2::new(2.0, 0.0)) - 1.0).abs() < 1e-6);
	}

	#[test]
	fn signed_distance_union() {
		let a = Region::Circle {
			center: Vec2::new(-1.0, 0.0),
			radius: 1.0,
		};
		let b = Region::Circle {
			center: Vec2::new(1.0, 0.0),
			radius: 1.0,
		};

		let u = Region::Union(Box::new(a), Box::new(b));

		assert!(u.signed_distance(Vec2::new(1.5, 0.0)) < 0.0);
		assert!(u.signed_distance(Vec2::new(3.0, 0.0)) > 0.0);
	}

	#[test]
	fn signed_distance_translate() {
		let c = Region::Circle {
			center: Vec2::new(0.0, 0.0),
			radius: 1.0,
		};
		let t = Region::Translate {
			region: Box::new(c),
			offset: Vec2::new(3.0, 0.0),
		};

		// (3,0) は中心 → d=-1
		assert!((t.signed_distance(Vec2::new(3.0, 0.0)) + 1.0).abs() < 1e-6);
		// (5,0) は外側 → d=1
		assert!((t.signed_distance(Vec2::new(5.0, 0.0)) - 1.0).abs() < 1e-6);
	}

	#[test]
	fn signed_distance_scale() {
		let c = Region::Circle {
			center: Vec2::new(0.0, 0.0),
			radius: 1.0,
		};
		let s = Region::Scale {
			region: Box::new(c),
			factor: 2.0,
		};

		// 半径2の円として動く
		assert!((s.signed_distance(Vec2::new(0.0, 0.0)) + 2.0).abs() < 1e-6);
		assert!((s.signed_distance(Vec2::new(3.0, 0.0)) - 1.0).abs() < 1e-6);
	}

	#[test]
	fn signed_distance_rotate() {
		let rect = Region::Rect {
			min: Vec2::new(-1.0, -0.5),
			max: Vec2::new(1.0, 0.5),
		};
		let r = Region::Rotate {
			region: Box::new(rect),
			radian: std::f64::consts::FRAC_PI_2, // 90度
		};

		// (0,1) は元の rect の (1,0) → inside
		assert!(r.signed_distance(Vec2::new(0.0, 0.9)) < 0.0);
		// (1,0) は元の rect の (0,-1) → outside
		assert!(r.signed_distance(Vec2::new(0.9, 0.0)) > 0.0);
	}
}
