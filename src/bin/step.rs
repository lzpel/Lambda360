use std::fs::File;
use std::io::{self, Write};

fn step_cylinder_csg(radius_mm: f64, height_mm: f64, name: &str) -> String {
    let r = format!("{:.6}", radius_mm);
    let h = format!("{:.6}", height_mm);

    format!(
"ISO-10303-21;
HEADER;
  FILE_DESCRIPTION(('STEP AP214'),'2;1');
  FILE_NAME('{name}.step','2025-11-21T12:00:00',('rust-generator'),('rust-generator'),'','', '');
  FILE_SCHEMA(('AUTOMOTIVE_DESIGN_CC2'));
ENDSEC;
DATA;

// ----- context & units -----
#100 = APPLICATION_CONTEXT('mechanical design');
#101 = APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#100);

#110 = (LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));
#111 = (PLANE_ANGLE_UNIT() NAMED_UNIT(*) SI_UNIT($,.RADIAN.));
#112 = (SOLID_ANGLE_UNIT() NAMED_UNIT(*) SI_UNIT($,.STERADIAN.));
#113 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-06),#110,'distance_accuracy_value','confusion');
#114 = (GEOMETRIC_REPRESENTATION_CONTEXT(3)
         GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#113))
         GLOBAL_UNIT_ASSIGNED_CONTEXT((#110,#111,#112))
         REPRESENTATION_CONTEXT('',''));

// ----- placement -----
#1 = DIRECTION('',(0.,0.,1.));
#2 = DIRECTION('',(1.,0.,0.));
#3 = CARTESIAN_POINT('',(0.,0.,0.));
#4 = AXIS2_PLACEMENT_3D('',#3,#1,#2);

// ----- CSG primitive: right circular cylinder -----
#10 = RIGHT_CIRCULAR_CYLINDER('',#4,{r},{h});
#11 = CSG_SOLID(#10);

// ----- product structure -----
#20 = PRODUCT('{name}','{name}','',(#100));
#21 = PRODUCT_DEFINITION_FORMATION('','',#20);
#22 = DESIGN_CONTEXT('design',#100,'design');
#23 = PRODUCT_DEFINITION('','',#21,#22);

#30 = SHAPE_REPRESENTATION('',(#11),#114);
#31 = PRODUCT_DEFINITION_SHAPE('','',#23);
#32 = SHAPE_DEFINITION_REPRESENTATION(#31,#30);

ENDSEC;
END-ISO-10303-21;
",
        name = name, r = r, h = h
    )
}

fn write_step(path: &str, contents: &str) -> io::Result<()> {
    let mut f = File::create(path)?;
    f.write_all(contents.as_bytes())
}

fn main() -> io::Result<()> {
    let txt = step_cylinder_csg(10.0, 50.0, "cylinder_csg");
    write_step("cylinder_csg.step", &txt)?;
    println!("Wrote cylinder_csg.step");
    Ok(())
}