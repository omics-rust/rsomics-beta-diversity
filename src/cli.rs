use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_beta_diversity::{CountTable, DistanceMatrix, Metric, Report};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-beta-diversity", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Count table (feature-by-sample TSV); reads stdin when "-" or omitted.
    #[arg(default_value = "-")]
    input: PathBuf,

    /// Distance metric: braycurtis, jaccard, euclidean, canberra, cityblock.
    #[arg(short = 'm', long, default_value = "braycurtis")]
    metric: String,

    /// Parse the input as comma-separated instead of tab-separated.
    #[arg(long, default_value_t = false)]
    csv: bool,

    /// Output path; writes stdout when "-".
    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    fn execute(&self) -> Result<Report> {
        let metric = Metric::parse(self.metric.trim())?;
        let delim = if self.csv { ',' } else { '\t' };

        let reader: Box<dyn std::io::BufRead> = if self.input.as_os_str() == "-" {
            Box::new(BufReader::new(std::io::stdin().lock()))
        } else {
            Box::new(BufReader::new(File::open(&self.input).map_err(|e| {
                RsomicsError::InvalidInput(format!("{}: {e}", self.input.display()))
            })?))
        };
        let table = CountTable::parse(reader, delim)?;
        let dm = DistanceMatrix::compute(&table, metric);

        if !self.common.json {
            let mut out: Box<dyn Write> = if self.output == "-" {
                Box::new(BufWriter::new(std::io::stdout().lock()))
            } else {
                Box::new(BufWriter::new(
                    File::create(&self.output).map_err(RsomicsError::Io)?,
                ))
            };
            dm.write_tsv(&mut out)?;
            out.flush().map_err(RsomicsError::Io)?;
        }

        Ok(dm.report(metric))
    }
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        Cli::execute(&self)?;
        Ok(())
    }

    // The default `run` discards the body's value, so `--json` would emit a
    // `result: null` envelope after the text matrix. Override to carry the
    // populated Report into the single envelope while leaving the non-json
    // path (the TSV matrix on stdout/-o) byte-for-byte intact.
    fn run(self) -> std::process::ExitCode {
        let common = self.common().clone();
        rsomics_common::run(&common, Self::meta(), move || Cli::execute(&self))
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Pairwise between-sample beta-diversity distance matrix from a feature count table.",
    origin: Some(Origin {
        upstream: "scikit-bio skbio.diversity.beta_diversity",
        upstream_license: "BSD-3-Clause",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["[table.tsv] [-m braycurtis] [-o dm.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('m'),
                long: "metric",
                aliases: &[],
                value: Some("<name>"),
                type_hint: Some("String"),
                required: false,
                default: Some("braycurtis"),
                description: "braycurtis | jaccard | euclidean | canberra | cityblock(manhattan).",
                why_default: Some("Bray-Curtis is the microbiome default"),
            },
            FlagSpec {
                short: None,
                long: "csv",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: Some("false"),
                description: "Parse the table as comma-separated.",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "Output path (- for stdout).",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Bray-Curtis distance matrix from a TSV count table",
            command: "rsomics-beta-diversity counts.tsv",
        },
        Example {
            description: "Jaccard (presence/absence) matrix from stdin, 8 threads",
            command: "cat counts.tsv | rsomics-beta-diversity -m jaccard -t 8",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }

    #[test]
    fn metric_aliases() {
        assert_eq!(Metric::parse("manhattan").unwrap(), Metric::Cityblock);
        assert_eq!(Metric::parse("bray").unwrap(), Metric::BrayCurtis);
        assert!(Metric::parse("unifrac").is_err());
    }
}
