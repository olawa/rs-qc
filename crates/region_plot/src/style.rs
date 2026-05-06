use plotters::style::RGBColor;

#[derive(Clone, Debug)]
pub struct PlotStyle {
    pub background: RGBColor,
    pub text: RGBColor,
    pub axis: RGBColor,
    pub coverage: RGBColor,
    pub read_forward: RGBColor,
    pub read_reverse: RGBColor,
    pub deletion: RGBColor,
    pub insertion: RGBColor,
    pub gene: RGBColor,
    pub gene_reverse: RGBColor,
    pub exon: RGBColor,
    pub exon_reverse: RGBColor,
    pub haplotype_1: RGBColor,
    pub haplotype_2: RGBColor,
    pub modification: RGBColor,
}

impl Default for PlotStyle {
    fn default() -> Self {
        Self {
            background: RGBColor(255, 255, 255),
            text: RGBColor(35, 35, 35),
            axis: RGBColor(120, 120, 120),
            coverage: RGBColor(184, 184, 184),
            read_forward: RGBColor(180, 180, 180),
            read_reverse: RGBColor(180, 180, 180),
            deletion: RGBColor(0, 0, 0),
            insertion: RGBColor(0, 0, 255),
            gene: RGBColor(100, 150, 220),
            gene_reverse: RGBColor(228, 105, 135),
            exon: RGBColor(40, 120, 220),
            exon_reverse: RGBColor(220, 88, 122),
            haplotype_1: RGBColor(48, 125, 190),
            haplotype_2: RGBColor(204, 121, 47),
            modification: RGBColor(255, 165, 0),
        }
    }
}
