use std::io::{self, Write};
use crate::{Float, Vec3};

/// バイナリ STL を書き出す。
/// - `tris`: 三角形配列
/// - `w`: 書き込み先（ファイル・メモリバッファなど）
/// - `header_text`: 先頭 80 バイトのヘッダに入れる文字列（None ならゼロ埋め）
///
/// 注意:
/// - STL の float は常に little-endian の f32
/// - normal は v0,v1,v2 から計算（ゼロ面積なら (0,0,0)）
pub fn write_stl_binary<W: Write>(
    tris: &[crate::Triangle],
    mut w: W,
    header_text: Option<&str>,
) -> io::Result<()> {
    // 80バイトヘッダ
    let mut header = [0u8; 80];
    if let Some(text) = header_text {
        let bytes = text.as_bytes();
        let n = bytes.len().min(80);
        header[..n].copy_from_slice(&bytes[..n]);
    }
    w.write_all(&header)?;

    // 三角形数（u32, little-endian）
    let count = tris.len() as u32;
    w.write_all(&count.to_le_bytes())?;

    for tri in tris {
        let [v0, v1, v2] = tri;

        // 法線を計算（右手系）
        let e1 = v1 - v0;
        let e2 = v2 - v0;
        let n = e1.cross(e2);
        let len2 = n.length_squared();

        let normal = if len2 > 0.0 {
            n / len2.sqrt()
        } else {
            // 退化三角形の場合は (0,0,0) にしておく
            Vec3::new(0.0, 0.0, 0.0)
        };

        // f64 (Float) → f32 に落として little-endian で書き込むヘルパ
        fn write_f32_le<W: Write>(w: &mut W, x: Float) -> io::Result<()> {
            let v = x as f32;
            w.write_all(&v.to_le_bytes())
        }

        // normal
        write_f32_le(&mut w, normal.x)?;
        write_f32_le(&mut w, normal.y)?;
        write_f32_le(&mut w, normal.z)?;

        // vertices
        write_f32_le(&mut w, v0.x)?;
        write_f32_le(&mut w, v0.y)?;
        write_f32_le(&mut w, v0.z)?;

        write_f32_le(&mut w, v1.x)?;
        write_f32_le(&mut w, v1.y)?;
        write_f32_le(&mut w, v1.z)?;

        write_f32_le(&mut w, v2.x)?;
        write_f32_le(&mut w, v2.y)?;
        write_f32_le(&mut w, v2.z)?;

        // attribute byte count（通常は 0）
        w.write_all(&0u16.to_le_bytes())?;
    }

    Ok(())
}