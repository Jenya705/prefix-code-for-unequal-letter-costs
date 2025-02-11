use std::{
    cmp::Reverse,
    collections::{BinaryHeap, HashMap, HashSet},
    io::Read,
    sync::{
        atomic::{AtomicU8, Ordering},
        Arc,
    },
};

use clap::{Parser, ValueEnum};
use length::{ESolver, GivenSolver, ILPSolver, LengthSolver, RelaxedLPSolver};
use scanner::Scanner;
use tree::{FirstToFindTreeSolver, TreeIndex, TreeOptimizer, TreeSolver};

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

    heap.push((Reverse(0u32), 0));

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
                i += 1;
                heap.push((Reverse(cost + letter), i));
                write
                    .write(
                        format!(r#""{:?}" -> "{:?}";"#, (i_c, cost), (i, cost + letter)).as_bytes(),
                    )
                    .unwrap();
            }
        }
    }

    write.write(b"}").unwrap();
}

fn read_input(
    format: AppInputFormat,
    scanner: &mut Scanner<impl std::io::BufRead>,
) -> (Vec<u32>, Vec<(char, u32)>) {
    let n = scanner.read::<usize>();
    let mut letters = Vec::with_capacity(n);
    for _ in 0..n {
        letters.push(scanner.read());
    }
    let occurences = match format {
        AppInputFormat::Message => {
            let message = scanner.read_line().to_string();
            let mut occurences = HashMap::new();
            for mut c in message.chars() {
                if c == '\n' {
                    c = ' ';
                }
                *occurences.entry(c).or_insert(0) += 1;
            }
            let occurences = occurences.into_iter().collect::<Vec<_>>();
            occurences
        }
        AppInputFormat::Occurences => {
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

type Attempt = (f64, Vec<u32>, Vec<u128>);

fn solve<L: LengthSolver, T: TreeSolver>(
    letters: &[u32],
    occurences: &[u32],
    time_limit: Option<u32>,
    mut attempts_history: Option<&mut Vec<Attempt>>,
) -> Option<Attempt> {
    let mut length_solver = L::new(&occurences, &letters);
    let mut tree_solver = T::new(letters);

    let mut max_cost = length_solver.theoretical_max();
    let mut lengths = vec![0; occurences.len()];
    let mut indices = vec![0; occurences.len()];

    let mut best_lengths = vec![0; occurences.len()];
    let mut best_indices = vec![0; occurences.len()];

    let mut tree_optimizer = TreeOptimizer::new(occurences.len());

    let start = std::time::Instant::now();

    let counter = Arc::new(AtomicU8::new(0));

    {
        let counter = Arc::clone(&counter);
        std::thread::spawn(move || loop {
            let mut buf = [0; 10];
            if let Ok(l) = std::io::stdin().lock().read(&mut buf) {
                const WORD: &[u8] = b"break";
                if l > WORD.len() && &buf[..WORD.len()] == WORD {
                    counter.store(1, Ordering::Relaxed);
                    return;
                }
                if counter.load(Ordering::Relaxed) == 2 {
                    return;
                }
            }
        });
    }

    loop {
        if matches!(time_limit, Some(limit) if start.elapsed().as_secs() > limit as _) {
            break;
        }

        if counter.load(Ordering::Relaxed) == 1 {
            break;
        }

        let Some(mut cost) = length_solver.solve(max_cost, &mut lengths) else {
            break;
        };

        let non_adjusted_cost = cost;

        if non_adjusted_cost > max_cost {
            if L::BEST_SOLVER {
                // we will find nothing better anyway
                break;
            } else {
                continue;
            }
        }

        let changed = match tree_solver.solve(&mut lengths, &mut indices) {
            Some(changed) => changed,
            None => {
                continue;
            }
        };

        let changed = tree_optimizer.optimize(&letters, &mut lengths, &mut indices) || changed;

        if changed {
            cost = 0.0;
            for i in 0..occurences.len() {
                cost += occurences[i] as f64 * lengths[i] as f64;
            }
        }

        if cost < max_cost {
            best_lengths.clone_from(&lengths);
            best_indices.clone_from(&indices);
        }

        if let Some(ref mut attempts_story) = attempts_history {
            attempts_story.push((cost, lengths.clone(), indices.clone()));
        }

        max_cost = cost.min(max_cost);

        if L::BEST_SOLVER && non_adjusted_cost == cost {
            break;
        }
    }

    counter.store(2, Ordering::Relaxed);

    if best_lengths[0] == 0 {
        None
    } else {
        Some((max_cost, best_lengths, best_indices))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum AppInputFormat {
    Message,
    Occurences,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum AppLengthSolver {
    ILP,
    RLP,
    E,
    Given,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum AppTreeSolver {
    FirstToFind,
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct App {
    #[arg(short, value_enum, default_value_t = AppLengthSolver::ILP)]
    length_solver: AppLengthSolver,

    #[arg(short, value_enum, default_value_t = AppInputFormat::Message)]
    input_format: AppInputFormat,

    #[arg(short, value_enum, default_value_t = AppTreeSolver::FirstToFind)]
    tree_solver: AppTreeSolver,

    #[arg(long = "og", default_value_t = false)]
    output_graph: bool,

    #[arg(long = "ol", default_value_t = true)]
    output_lengths: bool,

    #[arg(long = "oi", default_value_t = false)]
    output_indices: bool,

    #[arg(long = "ot", default_value_t = false)]
    output_table: bool,

    #[arg(short = 's', long, default_value_t = false)]
    history: bool,

    #[arg(long)]
    time_limit: Option<u32>,

    #[arg()]
    input_file: String,
}

fn output_attempt(attempt: &Attempt, letters: &[u32], occurences: &[(char, u32)], app: &App) {
    let (cost, lengths, indices) = attempt;

    println!(
        "cost: {}, avg. cost pro symbol: {:.3}",
        cost,
        cost / occurences.iter().map(|v| v.1).sum::<u32>() as f64
    );

    if app.output_graph {
        let mut exclude = HashSet::new();
        for i in 0..indices.len() {
            if indices[i] == 0 {
                break;
            }
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

    if app.output_lengths {
        println!(
            "lengths: {}",
            lengths
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(" ")
        );
    }

    if app.output_indices {
        println!(
            "indices: {}",
            indices
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(" ")
        );
    }

    if app.output_table {
        let mut to_output = vec![];
        for i in 0..indices.len() {
            print!("{}: ", occurences[i].0);
            let mut i = indices[i];
            loop {
                i -= 1;
                let j = i % letters.len() as u128;
                i /= letters.len() as u128;
                to_output.push(j as u32);
                if i == 0 {
                    break;
                }
            }
            for to_output in to_output.drain(..).rev() {
                print!("{to_output} ");
            }
            println!();
        }
    }
}

fn main() {
    let app = App::parse();

    let file = std::fs::File::open(&app.input_file).unwrap();
    let mut scanner = Scanner::new(std::io::BufReader::new(file));

    let (letters, mut occurences) = read_input(app.input_format, &mut scanner);
    occurences.sort_unstable_by_key(|v| Reverse(v.1));

    let occurences_count = occurences.iter().map(|v| v.1).collect::<Vec<_>>();

    let solve_fn = match (app.length_solver, app.tree_solver) {
        (AppLengthSolver::ILP, AppTreeSolver::FirstToFind) => {
            solve::<ILPSolver, FirstToFindTreeSolver>
        }
        (AppLengthSolver::E, AppTreeSolver::FirstToFind) => solve::<ESolver, FirstToFindTreeSolver>,
        (AppLengthSolver::Given, AppTreeSolver::FirstToFind) => {
            solve::<GivenSolver, FirstToFindTreeSolver>
        }
        (AppLengthSolver::RLP, AppTreeSolver::FirstToFind) => {
            solve::<RelaxedLPSolver, FirstToFindTreeSolver>
        }
    };

    let mut attempts_history = app.history.then(|| vec![]);

    let start = std::time::Instant::now();

    match solve_fn(
        &letters,
        &occurences_count,
        app.time_limit,
        attempts_history.as_mut(),
    ) {
        Some(attempt) => {
            output_attempt(&attempt, &letters, &occurences, &app);
        }
        None => {
            println!("found nothing");
        }
    }

    println!("time elapsed: {:?}", start.elapsed());

    if let Some(attempts_history) = attempts_history {
        let mut scanner = Scanner::new(std::io::BufReader::new(std::io::stdin().lock()));

        for i in 0..attempts_history.len() {
            println!(
                "{i}: {} {:?}",
                attempts_history[i].0, &attempts_history[i].1
            );
        }

        loop {
            let i = scanner.read::<isize>();
            if i == -1 {
                break;
            }
            let i = i as usize;
            output_attempt(&attempts_history[i], &letters, &occurences, &app);
        }
    }
}
