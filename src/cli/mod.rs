pub mod args;
pub mod commands;

pub use args::{
    AddReplaceRgArgs, AmpliconclipArgs, AmpliconstatsArgs, BedcovArgs, CalmdArgs, CatArgs,
    CollateArgs, ConsensusArgs, CoverageArgs, DepthArgs, DictArgs, FaidxArgs, FastqArgs,
    FixmateArgs, FlagstatArgs, FqidxArgs, HeadArgs, IdxstatsArgs, ImportArgs, IndexArgs,
    MarkdupArgs, MergeArgs, MpileupArgs, QuickcheckArgs, ReheaderArgs, SamplesArgs, SortArgs,
    SplitArgs, StatsArgs, TviewArgs, ViewArgs,
};
pub use commands::{
    cmd_addreplacerg, cmd_ampliconclip, cmd_ampliconstats, cmd_bedcov, cmd_calmd, cmd_cat,
    cmd_collate, cmd_consensus, cmd_coverage, cmd_depth, cmd_dict, cmd_faidx, cmd_fastq,
    cmd_fixmate, cmd_flagstat, cmd_fqidx, cmd_head, cmd_idxstats, cmd_import, cmd_index,
    cmd_markdup, cmd_merge, cmd_mpileup, cmd_quickcheck, cmd_reheader, cmd_samples, cmd_sort,
    cmd_split, cmd_stats, cmd_tview, cmd_view,
};
