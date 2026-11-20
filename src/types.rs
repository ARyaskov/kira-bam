use std::path::Path;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OutputFormat {
    Sam,
    Bam,
    UncompressedBam,
    Cram,
}

impl OutputFormat {
    pub fn from_path_ext(path: &Path) -> Self {
        match path.extension().and_then(|e| e.to_str()) {
            Some("sam") => OutputFormat::Sam,
            Some("bam") => OutputFormat::Bam,
            Some("cram") => OutputFormat::Cram,
            _ => OutputFormat::Bam,
        }
    }
}
pub(crate) fn set_header_order(
    header: &mut noodles_sam::Header,
    sort_order: &[u8],
    group_order: Option<&[u8]>,
) {
    use noodles_sam::header::record::value::map::header::tag;

    let fields = header
        .header_mut()
        .get_or_insert_with(Default::default)
        .other_fields_mut();
    fields.insert(tag::SORT_ORDER, sort_order.into());
    fields.shift_remove(&tag::SUBSORT_ORDER);
    fields.shift_remove(&tag::GROUP_ORDER);
    if let Some(group_order) = group_order {
        fields.insert(tag::GROUP_ORDER, group_order.into());
    }
}
