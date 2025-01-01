use std::{
    cmp::Reverse,
    collections::{BinaryHeap, HashMap, HashSet},
};

use length::{ESolver, ILPSolver, LengthSolver};
use scanner::Scanner;
use tree::{FirstToFindTreeSolver, TreeIndex, TreeSolver};

pub mod length;
pub mod scanner;
pub mod tree;

fn gen_tree(
    max: u32,
    letters: &[u32],
    exclude: &HashSet<TreeIndex>,
    mut write: impl std::io::Write,
) {
    let mut heap = BinaryHeap::new();

    write.write(b"DiGraph {").unwrap();

    heap.push((Reverse(0u32), 1));

    while let Some((Reverse(cost), mut i)) = heap.pop() {
        if exclude.contains(&i) {
            write
                .write(format!(r#""{:?}" [color=red];"#, (i, cost)).as_bytes())
                .unwrap();
            continue;
        }
        if cost < max {
            let i_c = i;
            i *= letters.len() as TreeIndex;
            for &letter in letters {
                heap.push((Reverse(cost + letter), i));
                write
                    .write(
                        format!(r#""{:?}" -> "{:?}";"#, (i_c, cost), (i, cost + letter)).as_bytes(),
                    )
                    .unwrap();
                i += 1;
            }
        }
    }

    write.write(b"}").unwrap();
}

#[derive(Clone, Copy)]
enum ChosenInputFormat {
    Message,
    Occurences,
}

#[derive(Clone, Copy)]
enum ChosenLengthSolver {
    ILP,
    E,
}

#[derive(Clone, Copy)]
enum ChosenTreeSolver {
    FirstToFind,
}

fn read_input(
    format: ChosenInputFormat,
    scanner: &mut Scanner<impl std::io::BufRead>,
) -> (Vec<u32>, Vec<(char, u32)>) {
    let n = scanner.read::<usize>();
    let mut letters = Vec::with_capacity(n);
    for _ in 0..n {
        letters.push(scanner.read());
    }
    let occurences = match format {
        ChosenInputFormat::Message => {
            let message = scanner.read_line().to_string();
            let mut occurences = HashMap::new();
            for c in message.chars() {
                *occurences.entry(c).or_insert(0) += 1;
            }
            let occurences = occurences.into_iter().collect::<Vec<_>>();
            occurences
        }
        ChosenInputFormat::Occurences => {
            let n = scanner.read::<usize>();
            let mut occurences = Vec::with_capacity(n);
            for _ in 0..n {
                occurences.push((scanner.read(), scanner.read()));
            }
            occurences
        }
    };
    (letters, occurences)
}

fn solve<L: LengthSolver, T: TreeSolver>(
    letters: &[u32],
    occurences: &[u32],
) -> Option<(f64, Vec<u32>, Vec<u128>)> {
    let m = 30; // max possible letter cost

    let mut length_solver = L::new(&occurences, &letters, m);
    let mut tree_solver = T::new(letters);

    let mut max_cost = length_solver.theoretical_max();
    let mut lengths = vec![0; occurences.len()];
    let mut indices = vec![0; occurences.len()];

    let mut best_lengths = vec![0; occurences.len()];
    let mut best_indices = vec![0; occurences.len()];

    loop {
        let Some(mut cost) = length_solver.solve(max_cost, &mut lengths) else {
            break;
        };

        let non_adjusted_cost = cost;

        println!("{cost} {lengths:?}");

        match tree_solver.solve(&mut lengths, &mut indices) {
            Some(true) => {
                cost = 0.0;
                for i in 0..occurences.len() {
                    cost += occurences[i] as f64 * lengths[i] as f64;
                }
            }
            Some(false) => {
                // Nothing
            }
            None => {
                continue;
            }
        }

        if cost < max_cost {
            best_lengths.clone_from(&lengths);
            best_indices.clone_from(&indices);
        }

        max_cost = cost;

        if L::BEST_SOLVER && (non_adjusted_cost == cost || non_adjusted_cost >= max_cost) {
            break;
        }
    }

    Some((max_cost, best_lengths, best_indices))
}

fn main() {
    let mut input_format = ChosenInputFormat::Message;
    let mut length_solver = ChosenLengthSolver::ILP;
    let mut tree_solver = ChosenTreeSolver::FirstToFind;

    let mut output_graph = false;
    let mut output_lengths = true;
    let mut output_indices = true;
    let mut output_codes = true;

    for arg in std::env::args().skip(1) {
        if arg == "--ilp" || arg == "--lp" {
            length_solver = ChosenLengthSolver::ILP;
        } else if arg == "--e" {
            length_solver = ChosenLengthSolver::E;
        } else if arg == "--msg" {
            input_format = ChosenInputFormat::Message;
        } else if arg == "--occ" {
            input_format = ChosenInputFormat::Occurences;
        } else if arg == "--ftf" {
            tree_solver = ChosenTreeSolver::FirstToFind;
        } else if arg == "--graph" {
            output_graph = !output_graph;
        } else if arg == "--lengths" {
            output_lengths = !output_lengths;
        } else if arg == "--indices" || arg == "--indexes" {
            output_indices = !output_indices;
        } else if arg == "--codes" {
            output_codes = !output_codes;
        } else {
            let file = std::fs::File::open(&arg).unwrap();
            let mut scanner = Scanner::new(std::io::BufReader::new(file));

            let (letters, mut occurences) = read_input(input_format, &mut scanner);

            occurences.sort_unstable_by_key(|v| Reverse(v.1));

            let occurences_count = occurences.iter().map(|v| v.1).collect::<Vec<_>>();

            let start = std::time::Instant::now();

            let solve_fn = match (length_solver, tree_solver) {
                (ChosenLengthSolver::ILP, ChosenTreeSolver::FirstToFind) => {
                    solve::<ILPSolver, FirstToFindTreeSolver>
                }
                (ChosenLengthSolver::E, ChosenTreeSolver::FirstToFind) => {
                    solve::<ESolver, FirstToFindTreeSolver>
                }
            };

            match solve_fn(&letters, &occurences_count) {
                Some((cost, lengths, indices)) => {
                    println!("cost: {}", cost);

                    if output_graph {
                        let mut exclude = HashSet::new();
                        for i in 0..indices.len() {
                            exclude.insert(indices[i]);
                        }
                        gen_tree(
                            *lengths.iter().max().unwrap(),
                            &letters,
                            &exclude,
                            std::io::stdout().lock(),
                        );
                        println!();
                    }

                    if output_lengths {
                        println!(
                            "lengths: {}",
                            lengths
                                .iter()
                                .map(|v| v.to_string())
                                .collect::<Vec<_>>()
                                .join(" ")
                        );
                    }

                    if output_indices {
                        println!(
                            "indices: {}",
                            indices
                                .iter()
                                .map(|v| v.to_string())
                                .collect::<Vec<_>>()
                                .join(" ")
                        );
                    }

                    if output_codes {
                        for i in 0..indices.len() {
                            print!("{}: ", occurences[i].0);
                            let mut i = indices[i];
                            loop {
                                let j = i % letters.len() as u128;
                                i /= letters.len() as u128;
                                print!("{} ", j);
                                if i == 1 {
                                    break;
                                }
                            }
                            println!();
                        }
                    }
                }
                None => {
                    println!("found nothing");
                }
            }

            println!("time elapsed: {:?}", start.elapsed());
        }
    }
}
