use anyhow::Result;
use clap::{Parser, Subcommand};

use kira_bam::cli::{
    AddReplaceRgArgs, AmpliconclipArgs, AmpliconstatsArgs, BedcovArgs, CalmdArgs, CatArgs,
    CollateArgs, ConsensusArgs, CoverageArgs, DepthArgs, DictArgs, FaidxArgs, FastqArgs,
    FixmateArgs, FlagstatArgs, FqidxArgs, HeadArgs, IdxstatsArgs, ImportArgs, IndexArgs,
    MarkdupArgs, MergeArgs, MpileupArgs, QuickcheckArgs, ReheaderArgs, SamplesArgs, SortArgs,
    SplitArgs, StatsArgs, TviewArgs, ViewArgs, cmd_addreplacerg, cmd_ampliconclip,
    cmd_ampliconstats, cmd_bedcov, cmd_calmd, cmd_cat, cmd_collate, cmd_consensus, cmd_coverage,
    cmd_depth, cmd_dict, cmd_faidx, cmd_fastq, cmd_fixmate, cmd_flagstat, cmd_fqidx, cmd_head,
    cmd_idxstats, cmd_import, cmd_index, cmd_markdup, cmd_merge, cmd_mpileup, cmd_quickcheck,
    cmd_reheader, cmd_samples, cmd_sort, cmd_split, cmd_stats, cmd_tview, cmd_view,
};

#[derive(Parser)]
#[command(name = "kira-bam")]
#[command(
    about = "High-performance BAM/SAM/CRAM toolkit (samtools-compatible)",
    version,
    author
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    View(ViewArgs),
    Flagstat(FlagstatArgs),
    Sort(SortArgs),
    Index(IndexArgs),
    Merge(MergeArgs),
    Markdup(MarkdupArgs),
    Fixmate(FixmateArgs),
    Idxstats(IdxstatsArgs),
    Quickcheck(QuickcheckArgs),
    Depth(DepthArgs),
    Coverage(CoverageArgs),
    Stats(StatsArgs),
    Fastq(FastqArgs),
    Calmd(CalmdArgs),
    Head(HeadArgs),
    Cat(CatArgs),
    Reheader(ReheaderArgs),
    Dict(DictArgs),
    Split(SplitArgs),
    Addreplacerg(AddReplaceRgArgs),
    Collate(CollateArgs),
    Bedcov(BedcovArgs),
    Samples(SamplesArgs),
    Consensus(ConsensusArgs),
    Mpileup(MpileupArgs),
    Faidx(FaidxArgs),
    Fqidx(FqidxArgs),
    Ampliconclip(AmpliconclipArgs),
    Ampliconstats(AmpliconstatsArgs),
    Import(ImportArgs),
    Tview(TviewArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::View(a) => cmd_view(a),
        Commands::Flagstat(a) => cmd_flagstat(a),
        Commands::Sort(a) => cmd_sort(a),
        Commands::Index(a) => cmd_index(a),
        Commands::Merge(a) => cmd_merge(a),
        Commands::Markdup(a) => cmd_markdup(a),
        Commands::Fixmate(a) => cmd_fixmate(a),
        Commands::Idxstats(a) => cmd_idxstats(a),
        Commands::Quickcheck(a) => cmd_quickcheck(a),
        Commands::Depth(a) => cmd_depth(a),
        Commands::Coverage(a) => cmd_coverage(a),
        Commands::Stats(a) => cmd_stats(a),
        Commands::Fastq(a) => cmd_fastq(a),
        Commands::Calmd(a) => cmd_calmd(a),
        Commands::Head(a) => cmd_head(a),
        Commands::Cat(a) => cmd_cat(a),
        Commands::Reheader(a) => cmd_reheader(a),
        Commands::Dict(a) => cmd_dict(a),
        Commands::Split(a) => cmd_split(a),
        Commands::Addreplacerg(a) => cmd_addreplacerg(a),
        Commands::Collate(a) => cmd_collate(a),
        Commands::Bedcov(a) => cmd_bedcov(a),
        Commands::Samples(a) => cmd_samples(a),
        Commands::Consensus(a) => cmd_consensus(a),
        Commands::Mpileup(a) => cmd_mpileup(a),
        Commands::Faidx(a) => cmd_faidx(a),
        Commands::Fqidx(a) => cmd_fqidx(a),
        Commands::Ampliconclip(a) => cmd_ampliconclip(a),
        Commands::Ampliconstats(a) => cmd_ampliconstats(a),
        Commands::Import(a) => cmd_import(a),
        Commands::Tview(a) => cmd_tview(a),
    }
}
