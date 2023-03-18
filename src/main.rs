#![allow(clippy::too_many_arguments)]

use clap::{Parser, ValueEnum};
use problem_1::Problem1;
use problem_2::Problem2;
use problem_4::Problem4;

mod integration;
mod matrices;
mod polynomials;
mod problem_1;
mod problem_2;
mod problem_3;
mod problem_4;
mod utils;

#[derive(Parser)]
#[clap(author, version, about, long_about = None, rename_all = "snake-case")]
pub enum Cli {
    Problem1(Problem1),
    Problem2(Problem2),
    Problem4(Problem4),
}

#[derive(Parser, Clone, ValueEnum)]
#[clap(rename_all = "snake-case")]
enum FunctionType {
    U,
    V,
    SigmaX,
    SigmaY,
}

#[derive(Parser, Clone, ValueEnum)]
#[clap(rename_all = "snake-case")]
enum LoadFunction {
    // function: (x - 2.5) * (x - 2.5)
    Type1,
    // function: (x - 2.5) * (x - 2.5) * (x - 2.5)
    Type2,
}

impl LoadFunction {
    fn call(&self, x: f64) -> f64 {
        match self {
            Self::Type1 => (x - 2.5) * (x - 2.5),
            Self::Type2 => (x - 2.5) * (x - 2.5) * (x - 2.5),
        }
    }
}

impl Cli {
    fn exec(self) {
        match self {
            Self::Problem1(problem1) => problem1.exec(),
            Self::Problem2(problem2) => problem2.exec(),
            Self::Problem4(problem4) => problem4.exec(),
        }
    }
}

fn main() {
    Cli::parse().exec()
}
