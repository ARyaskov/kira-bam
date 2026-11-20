use clap::Parser;
use kira_bam::cli::SortArgs;

#[test]
fn library_sort_does_not_initialize_the_global_rayon_pool() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("input.sam");
    let bam = dir.path().join("sorted.bam");
    std::fs::write(&sam, "@HD\tVN:1.6\nr\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n").unwrap();
    let args = SortArgs::try_parse_from([
        "sort",
        "-@",
        "1",
        "-o",
        bam.to_str().unwrap(),
        sam.to_str().unwrap(),
    ])
    .unwrap();
    kira_bam::sort::run(args).unwrap();
    assert!(bam.is_file());
    rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build_global()
        .expect("a library sort must leave the host's global Rayon configuration available");
    assert_eq!(rayon::current_num_threads(), 2);
}
