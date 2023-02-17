#![allow(clippy::too_many_arguments)]

use clap::Parser;
use problem_1::Problem1;

pub mod integration;
pub mod problem_1;
pub mod problem_3;
pub mod utils;

#[derive(Parser)]
#[clap(author, version, about, long_about = None, rename_all = "snake-case")]
pub enum Cli {
    Problem1(Problem1),
}

impl Cli {
    fn exec(self) {
        match self {
            Self::Problem1(problem1) => problem1.exec(),
        }
    }
}

fn main() {
    Cli::parse().exec()
}
