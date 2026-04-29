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
    pub exon: RGBColor,
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
            coverage: RGBColor(71, 123, 191),
            read_forward: RGBColor(86, 156, 109),
            read_reverse: RGBColor(191, 103, 92),
            deletion: RGBColor(40, 40, 40),
            insertion: RGBColor(119, 74, 157),
            gene: RGBColor(76, 88, 102),
            exon: RGBColor(43, 65, 98),
            haplotype_1: RGBColor(48, 125, 190),
            haplotype_2: RGBColor(204, 121, 47),
            modification: RGBColor(210, 167, 38),
        }
    }
}
