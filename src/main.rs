use std::{
    cmp::Reverse,
    collections::{BinaryHeap, HashMap, HashSet},
    io::Write,
};

use good_lp::{
    default_solver, variable, variables, Constraint, Expression, Solution, SolverModel, Variable,
};
use scanner::Scanner;

pub mod length;
pub mod scanner;

#[derive(Debug)]
pub struct InputData {
    pub balls: Vec<u32>,
    pub message: String,
}

impl InputData {
    pub fn read(scanner: &mut Scanner<impl std::io::BufRead>) -> Self {
        let n = scanner.read::<usize>();
        let mut balls = Vec::with_capacity(n);
        for _ in 0..n {
            balls.push(scanner.read());
        }
        let message = scanner.read_line().to_string();
        Self { balls, message }
    }
}

#[derive(Debug)]
struct Symbol {
    count: u32,
    symbol: char,
}

fn ilp_solve(occurences: &[u32], v_cacher: &mut VCacher, n: usize, m: usize) {
    variables! {vars: };

    let mut vars_mat = vec![];

    for _ in 0..n {
        let mut v = vec![];
        for _ in 0..m {
            v.push(vars.add(variable().binary()));
        }
        vars_mat.push(v);
    }

    let mut constraints = vec![];

    for j in 0..m {
        let mut expr = Expression::default();
        for p in 0..j {
            for i in 0..n {
                expr = expr + vars_mat[i][p] * v_cacher.get(j as i32 - p as i32);
            }
        }
        constraints.push(expr.leq(v_cacher.get(j as i32)));
    }

    for i in 0..n {
        let mut expr = Expression::default();
        for j in 0..m {
            expr = expr + vars_mat[i][j];
        }
        constraints.push(expr.eq(1));
    }

    let mut minimise_expr = Expression::default();

    for i in 0..n {
        for j in 0..m {
            minimise_expr = minimise_expr + vars_mat[i][j] * (occurences[i] * (j + 1) as u32);
        }
    }

    let mut problem = vars.minimise(minimise_expr).using(default_solver);

    for constraint in constraints {
        problem = problem.with(constraint);
    }

    let solution = problem.solve().unwrap();

    for i in 0..n {
        for j in 0..m {
            if solution.value(vars_mat[i][j]) > 0.5 {
                print!("{} ", j);
            }
        }
    }
    println!();
}

fn solve(data: InputData) {
    let mut occurences = HashMap::new();

    for message in data.message.chars() {
        *occurences.entry(message).or_insert(0) += 1;
    }

    println!("{:?}", occurences);

    let mut symbols = occurences
        .into_iter()
        .map(|(symbol, occurences)| Symbol {
            count: occurences,
            symbol,
        })
        .collect::<Vec<_>>();

    // bigger occurences number = bigger probability
    symbols.sort_unstable_by_key(|symbol| Reverse(symbol.count));

    println!("{:?}", symbols);

    let t = calculate_t(&data);

    println!("{:?}", t);

    let each_symbols_count = symbols.iter().map(|v| v.count as f64).collect::<Vec<_>>();
    let symbols_count = each_symbols_count.iter().sum::<f64>();
    let n = symbols.len();
    let m = 100;
    let goal = (n as f64 - 1.0) / (1.0 - t);

    let mut sorted = Vec::with_capacity(n * m);
    for i in 0..n {
        for j in 0..m {
            sorted.push((i, j, each_symbols_count[i] * t.powi(-(j as i32 + 1))));
        }
    }
    sorted.sort_unstable_by(|v1, v2| {
        v1.2.partial_cmp(&v2.2)
            .unwrap()
            .then(v2.0.cmp(&v1.0))
            .then(v2.1.cmp(&v1.1))
    });
    let sorted = sorted.into_iter().map(|v| (v.0, v.1)).collect::<Vec<_>>();

    let mut t_powers = vec![1.0; m];
    for i in 1..m {
        t_powers[i] = (t_powers[i - 1]) * t;
    }

    let mut rows = vec![0; n];

    let sub = symbols_count.log(t);

    let e_min = each_symbols_count
        .iter()
        .map(|&v| v * (v.log(t) - sub))
        .sum::<f64>()
        / 0.90;

    println!("{e_min}");

    let mut v_cacher = VCacher {
        cache: vec![1],
        letters: &data.balls,
    };

    let mut best_e = e_min;
    let mut best_result = vec![];
    let mut best_rows = vec![];

    for j0 in 0..m {
        'e_calc: for i0 in 0..n {
            let subject0 = calculate_subject0(&t_powers, i0, j0, n);

            // if subject0 > goal {
            //     continue;
            // }

            // e represent cost of representing the MESSAGE itself and not a symbol on average
            let Some(e) = calculate_e(
                &each_symbols_count,
                t,
                &t_powers,
                i0,
                j0,
                n,
                m,
                goal,
                &sorted,
                &mut rows,
                subject0,
                best_e,
            ) else {
                continue 'e_calc;
            };

            if e > e_min {
                continue 'e_calc;
            }

            // rows array is sorted
            // let mut pr = rows[0];
            // let mut count = 1;
            // for i in 1..rows.len() {
            //     if pr == rows[i] {
            //         count += 1;
            //     } else {
            //         if v_cacher.get(pr as i32) < count {
            //             println!("Skip due not having enough ways of representing letters with length {}", pr);
            //             continue 'e_calc;
            //         }
            //         count = 0;
            //         pr = rows[i];
            //     }
            // }

            // if v_cacher.get(pr as i32) < count {
            //     println!(
            //         "Skip due not having enough ways of representing letters with length {}",
            //         pr
            //     );
            //     continue 'e_calc;
            // }

            // let mut res = vec![0; rows.len()];
            // let mut dfs = DFS {
            //     rows: &rows,
            //     letters: &data.balls,
            //     set: &mut HashSet::with_capacity(rows.len()),
            //     res: &mut res,
            // };

            // let res = dfs.dfs(0, 0, 0);
            // println!("{i0} {j0} {e} {res} {:?}", dfs.res);

            let res = tree_search(&mut rows, &data.balls);

            let new_e = rows
                .iter()
                .zip(&each_symbols_count)
                .map(|(&row, &count)| row as f64 * count)
                .sum::<f64>();

            if e != new_e {
                println!("{e} {new_e} adjusted");
            }

            if new_e < best_e {
                best_e = e;
                best_result.clone_from(&res);
                best_rows.clone_from(&rows);
            }

            // let mut exclude = HashSet::new();
            // for e in res.iter() {
            //     exclude.insert(*e);
            // }
            // gen_tree(
            //     *rows.last().unwrap(),
            //     &data.balls,
            //     &exclude,
            //     std::io::stdout().lock(),
            // );
            // println!();

            // exclude.clear();
            // for e in dfs.res.iter() {
            //     exclude.insert(*e as u32);
            // }
            // gen_tree(
            //     *rows.last().unwrap(),
            //     &data.balls,
            //     &exclude,
            //     std::io::stdout().lock(),
            // );
            // println!();

            // if res.len() == rows.len() {
            //     println!("success: {e} {:?} {:?}", rows, res);
            // } else {
            //     println!("failure: {e} {:?} {:?}", rows, res);
            // }
        }
    }

    let mut exclude = HashSet::new();
    for e in best_result.iter() {
        exclude.insert(*e);
    }
    gen_tree(
        *best_rows.last().unwrap(),
        &data.balls,
        &exclude,
        std::io::stdout().lock(),
    );
    println!();
    println!("success: {best_e} {:?} {:?}", best_rows, best_result);

    // let v = calculate_e(&each_symbols_count, t, &t_powers, 1, 2, n, m, goal, &sorted);
    // println!("{:?}", (v.0 / symbols_count, v.1),);
}

struct VCacher<'a> {
    cache: Vec<u32>,
    letters: &'a [u32],
}

impl<'a> VCacher<'a> {
    pub fn get(&mut self, h: i32) -> u32 {
        if h < 0 {
            0
        } else if self.cache.len() as i32 > h {
            self.cache[h as usize]
        } else {
            for h in self.cache.len() as i32..h {
                let _ = self.get(h);
            }
            let mut v = 0;
            for i in 0..self.letters.len() {
                v += self.get(h - self.letters[i] as i32);
            }
            self.cache.push(v);
            v
        }
    }
}

struct DFS<'a> {
    rows: &'a [u32],
    letters: &'a [u32],
    set: &'a mut HashSet<usize>,
    res: &'a mut [usize],
}

impl<'a> DFS<'a> {
    fn dfs(&mut self, row_i: usize, mut set_index: usize, current_cost: u32) -> bool {
        set_index *= self.letters.len() + 1;
        for i in 0..self.letters.len() {
            set_index += 1;

            // println!("{row_i} {set_index} {current_cost}");

            if self.set.contains(&set_index) {
                continue;
            }

            let new_cost = current_cost + self.letters[i];

            if new_cost > self.rows[row_i] {
                // TODO: sort costs
                break;
            } else if new_cost == self.rows[row_i] {
                self.set.insert(set_index);
                if row_i == self.rows.len() - 1 || self.dfs(row_i + 1, 0, 0) {
                    self.res[row_i] = set_index;
                    return true;
                }
                self.set.remove(&set_index);
            } else {
                if self.dfs(row_i, set_index, new_cost) {
                    return true;
                }
            }
        }

        false
    }
}

fn gen_tree(max: u32, letters: &[u32], exclude: &HashSet<u64>, mut write: impl std::io::Write) {
    let mut heap = BinaryHeap::new();

    write.write(b"DiGraph {").unwrap();

    heap.push((Reverse(0u32), 0u64));

    while let Some((Reverse(cost), mut i)) = heap.pop() {
        if exclude.contains(&i) {
            write
                .write(format!(r#""{:?}" [color=red];"#, (i, cost)).as_bytes())
                .unwrap();
            continue;
        }
        if cost < max {
            let i_c = i;
            i *= letters.len() as u64 + 1;
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

/// dijkstra like algorithm
fn tree_search(mut rows: &mut [u32], letters: &[u32]) -> Vec<u64> {
    let mut result = vec![];
    let mut heap = BinaryHeap::new();
    heap.push((Reverse(0u32), 0u64));

    while let Some((Reverse(cost), mut i)) = heap.pop() {
        if cost < rows[0] {
            i *= letters.len() as u64 + 1;
            for &letter in letters {
                i += 1;
                heap.push((Reverse(cost + letter), i));
            }
        } else {
            result.push(i);
            rows[0] = cost;
            rows = &mut rows[1..];
            if rows.len() == 0 {
                break;
            }
        }
    }

    result
}

fn calculate_subject0(t_powers: &[f64], i0: usize, j0: usize, n: usize) -> f64 {
    let mut subject = 0.0;
    for j in 0..=j0 {
        subject += t_powers[j];
    }
    subject *= (n - i0) as f64;
    subject
}

fn calculate_e(
    symbols: &[f64],
    t: f64,
    t_powers: &[f64],
    i0: usize,
    j0: usize,
    n: usize, // max i
    m: usize, // max j
    goal: f64,
    sorted: &[(usize, usize)],
    rows: &mut [u32],
    mut subject: f64,
    res_to_beat: f64,
) -> Option<f64> {
    rows.fill(0);

    let mut res = 0.0;

    for i in i0..n {
        rows[i] = j0 as u32 + 1;
        res += symbols[i] * rows[i] as f64;
    }

    // println!("{subject}");

    // println!("{:?}", rows);

    for (si, &(i, j)) in sorted
        .iter()
        .filter(|&&(i, j)| !(i <= i0 && j > j0) && !(i >= i0 && j <= j0))
        .enumerate()
    {
        subject += t_powers[j];
        rows[i] += 1;
        res += symbols[i];

        if res > res_to_beat {
            break;
        }

        // let mut rec_subject = 0.0;
        // for i in 0..n {
        //     rec_subject += rows[i].0 * rows[i].1;
        // }
        // if rec_subject != subject {
        //     println!("{si}. {rec_subject} {subject}");
        // }

        if subject >= goal {
            // for i in 0..n {
            //     println!("{:?}", mat[i]);
            // }
            // println!("{:?}", rows);
            return Some(res);
        }
    }

    None
}

fn calculate_t(data: &InputData) -> f64 {
    // using newton's method
    let mut t0 = 1.0f64;

    for _ in 0..10 {
        let mut f_t0 = -1.0;
        let mut fdt_t0 = 0.0;

        for &ball in data.balls.iter() {
            f_t0 += t0.powi(ball as i32);
            // derivative of the power function
            fdt_t0 += ball as f64 * t0.powi(ball as i32 - 1);
        }

        // step
        t0 -= f_t0 / fdt_t0;
    }

    t0
}

fn main() {
    // let mut message = String::new();

    // let mut occurences = vec![];

    // for (char, count) in [
    //     (' ', 2000),
    //     ('E', 1050),
    //     ('T', 720),
    //     ('O', 654),
    //     ('A', 630),
    //     ('N', 590),
    //     ('I', 550),
    //     ('S', 540),
    //     ('H', 520),
    //     ('R', 470),
    //     ('D', 350),
    //     ('L', 290),
    //     ('C', 230),
    //     ('U', 225),
    //     ('F', 225),
    //     ('M', 210),
    //     ('P', 175),
    //     ('W', 120),
    //     ('Y', 120),
    //     ('G', 110),
    //     ('B', 105),
    //     ('V', 80),
    //     ('K', 30),
    //     ('X', 20),
    //     ('J', 10),
    //     ('Q', 10),
    //     ('Z', 10),
    // ] {
    //     occurences.push(count as u32);
    //     // message.push_str(&char.to_string().repeat(count));
    // }

    // ilp_solve(
    //     &occurences,
    //     &mut VCacher {
    //         cache: vec![1],
    //         letters: &[2, 3, 3],
    //     },
    //     occurences.len(),
    //     30,
    // );

    // solve(InputData {
    //     balls: vec![2, 3, 3],
    //     message,
    // });

    enum InputFormat {
        Message,
        Occurences,
    }

    enum LengthSolver {
        ILP,
        E,
    }

    let mut input_format = InputFormat::Message;
    let mut length_solver = LengthSolver::ILP;

    for arg in std::env::args().skip(1) {
        if arg == "--ilp" || arg == "--lp" {
            length_solver = LengthSolver::ILP;
        } else if arg == "--e" {
            length_solver = LengthSolver::E;
        } else if arg == "--msg" {
            input_format = InputFormat::Message;
        } else if arg == "--occ" {
            input_format = InputFormat::Occurences;
        } else {
            let file = std::fs::File::open(&arg).unwrap();
            let mut scanner = Scanner::new(std::io::BufReader::new(file));

            let data = InputData::read(&mut scanner);

            let start = std::time::Instant::now();

            solve(data);

            println!("time elapsed: {:?}", start.elapsed());
        }
    }
}
