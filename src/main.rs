#![allow(clippy::too_many_arguments)]

use clap::{Parser, ValueEnum};
use problem_1::Problem1;
use problem_3::Problem3;

pub mod integration;
pub mod problem_1;
pub mod problem_3;
pub mod utils;

#[derive(Parser)]
#[clap(author, version, about, long_about = None, rename_all = "snake-case")]
pub enum Cli {
    Problem1(Problem1),
    Problem3(Problem3),
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
            Self::Problem3(problem3) => problem3.exec(),
        }
    }
}

fn main() {
    Cli::parse().exec()
}
