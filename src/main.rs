#![allow(clippy::too_many_arguments)]

use std::fmt::Display;

use clap::{Parser, ValueEnum};
use problem_1::Problem1;
use problem_1_2::Problem1_2;
use problem_2::Problem2;
use problem_3::Problem3;

mod integration;
mod matrices;
mod polynomials;
mod problem_1;
mod problem_1_2;
mod problem_2;
mod problem_3;
mod utils;

#[derive(Parser)]
#[clap(author, version, about, long_about = None, rename_all = "snake-case")]
pub enum Cli {
    Problem1(Problem1),
    Problem1_2(Problem1_2),
    Problem2(Problem2),
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

impl Display for FunctionType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::U => write!(f, "u"),
            Self::V => write!(f, "v"),
            Self::SigmaX => write!(f, "sigma_x"),
            Self::SigmaY => write!(f, "sigma_y"),
        }
    }
}

#[derive(Parser, Clone, ValueEnum)]
#[clap(rename_all = "snake-case")]
enum LoadFunction {
    // function: (x - 2.5) * (x - 2.5)
    Type1,
    // function: (x - 2.5) * (x - 2.5) * (x - 2.5)
    Type2,
}

impl Display for LoadFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Type1 => write!(f, "(x - 2.5) * (x - 2.5)"),
            Self::Type2 => write!(f, "(x - 2.5) * (x - 2.5) * (x - 2.5)"),
        }
    }
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
            Self::Problem1_2(problem2) => problem2.exec(),
            Self::Problem2(problem3) => problem3.exec(),
            Self::Problem3(problem4) => problem4.exec(),
        }
    }
}

fn main() {
    Cli::parse().exec()
}
