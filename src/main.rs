#![feature(get_mut_unchecked)]

use std::{fs::File, io::BufReader};

use clap::Parser;
use input::Input;
use scanner::Scanner;
use solver::try_solve;

mod huffman;
mod input;
mod scanner;
mod solver;
mod tree;

#[derive(Parser)]
struct App {
    #[arg(short, long)]
    /// The other format for input
    occurences: bool,
    #[arg(short, long, default_value_t = 0.0)]
    /// The factor on which each probability will be multiplied (in some cases it can lead to better results)
    ///
    /// If the value is 0, then instead of probabilities occurences will be used.
    mul: f64,
    #[arg(short = 'u', long, default_value_t = false)]
    /// Whether the Huffman-Algorithm will be used.
    huffman: bool,
    #[arg()]
    file_name: String,
}

fn main() {
    let app = App::parse();

    let input = Input::read(
        !app.occurences,
        app.mul,
        &mut Scanner::new(BufReader::new(File::open(app.file_name).unwrap())),
    );

    try_solve(
        &input.letters,
        &input.probabilities,
        &input.symbols,
        if app.mul == 0.0 { input.sum } else { app.mul },
        app.huffman,
    );
}
