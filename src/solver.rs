use std::time::Instant;

use good_lp::{
    Expression, ProblemVariables, ResolutionError, Solution, SolverModel, Variable,
    VariableDefinition,
};

use crate::{
    huffman::{huffman_optimal, huffman_solve},
    tree::{adjust_costs, build_prefixes},
};

/// der in der Doku beschriebene Wert L
const LIMIT: u32 = u32::MAX;

pub struct VCacher {
    values: Vec<u32>,
}

impl VCacher {
    pub fn new(letters: &[u32], m: u32, n: u32) -> Self {
        let mut values = vec![1];
        'o: for i in 1usize..=(m as usize) {
            let mut v = 0u32;
            for &letter in letters {
                if let Some(&to_add) = values.get(i - letter as usize) {
                    match v.checked_add(to_add) {
                        Some(nv) => v = nv,
                        None => break 'o,
                    }
                }
            }
            if v > LIMIT || values[i - 1] * n > LIMIT {
                break;
            }
            values.push(v);
        }
        Self { values }
    }

    pub fn m(&self) -> u32 {
        self.values.len() as u32
    }

    pub fn get(&self, i: usize) -> u32 {
        self.values[i]
    }
}

pub fn calculate_max_lengths(cost: f64, probabilities: &[f64], m: u32, lengths: &mut [u32]) {
    let mut f = 0.0;
    let mut s = probabilities.iter().sum::<f64>();
    for i in 0..lengths.len() {
        let mut l = 1;
        let mut r = m;

        // wenn i = 0, dann ist i - 1 = -1
        if i != 0 {
            f += probabilities[i - 1];
            s -= probabilities[i - 1];
        }
        
        // f = p_1 + ... + p_{i-1}
        // s = p_i + ... + p_n

        // Binäresuche in (de facto) einer geordneten Liste, wobei die Elemente gleich
        // den Ausdruck aus der Sektion über Optimierung vom ILP-Problem sind.
        while r - l > 1 {
            let m = (l + r) / 2;
            let c = f + s * (m as f64);
            if c >= cost {
                r = m;
            } else {
                l = m;
            }
        }

        // die obere Grenze ist gebraucht (also r)
        // (l ist die untere Grenze)
        lengths[i] = r;
    }
}

pub fn try_solve(
    letters: &[u32],
    probabilities: &[f64],
    symbols: &[(char, u64)],
    mul: f64,
    use_huffman: bool,
) {
    let timer = Instant::now();

    let mut lengths = vec![0; probabilities.len()];

    if use_huffman {
        println!("huffman optimal: {}", huffman_optimal(letters));
        huffman_solve(letters, probabilities, &mut lengths);
        // kann den PC noch optimieren
        lengths.sort_unstable();
    } else {
        let mut adjusted_costs = vec![0; probabilities.len()];
        adjust_costs(letters, &mut adjusted_costs);

        let m = adjusted_costs.last().cloned().unwrap();
        let v_cacher = VCacher::new(letters, m, probabilities.len() as u32);
        let m = m.min(v_cacher.m());

        // Kalkuliert Kosten des PC aus den gegebenen Längen der Präfixe
        let calc_costs = |costs: &[u32]| -> f64 {
            costs
                .iter()
                .zip(probabilities.iter())
                .map(|(&f, &s)| f as f64 * s)
                .sum::<f64>()
        };

        // Es kann theoretisch sein, dass der von der Bestimmung des m PC besser als vom Huffman ist.
        let mut best_cost = calc_costs(&adjusted_costs);
        let mut best_adjusted_costs = adjusted_costs.clone();

        let mut costs = vec![0; probabilities.len()];
        huffman_solve(letters, probabilities, &mut costs);
        costs.sort_unstable();
        let cost = calc_costs(&costs);
        // Wählen der beste PC von den beiden
        if cost < best_cost {
            best_cost = cost;
            best_adjusted_costs = costs;
        }

        println!("best pre-ilp found cost: {}", best_cost / mul);

        let mut theoretical_max_lengths = vec![0; probabilities.len()];

        calculate_max_lengths(best_cost, probabilities, m, &mut theoretical_max_lengths);

        let result = ilp_solve(
            &v_cacher,
            letters,
            m,
            probabilities,
            &theoretical_max_lengths,
            &mut lengths,
        );

        match result {
            // ResolutionError::Infeasible meint dass es keine Lösung gefunden wurde, weil
            // die Lösung (im Fall) anhand der Optimierung des Problems nicht zu finden ist.
            Ok(_) | Err(ResolutionError::Infeasible) => {}
            Err(ref err) => println!("ILP Solver error: {:?}", err),
        }

        // Entweder ein Fehler oder der PC, der durch Huffman gefunden wurde, ist besser
        if result.is_err() || calc_costs(&lengths) > best_cost {
            lengths = best_adjusted_costs;
            println!("using already calculated result as the ilp solver didn't find any better possibilities");
        }
    }

    let mut prefixes = vec![vec![]; probabilities.len()];
    println!(
        "avg cost: {}, cost: {}",
        probabilities
            .iter()
            .zip(lengths.iter())
            .map(|(&f, &s)| f / mul * s as f64)
            .sum::<f64>(),
        symbols
            .iter()
            .zip(lengths.iter())
            .map(|(f, &s)| f.1 * s as u64)
            .sum::<u64>(),
    );
    println!("lengths: {lengths:?}");
    build_prefixes(letters, &lengths, &mut prefixes);
    for i in 0..probabilities.len() {
        println!("'{}': {:?}", symbols[i].0, prefixes[i]);
    }

    println!("time elapsed: {:?}", timer.elapsed());
}

pub fn ilp_solve(
    v: &VCacher,
    _letters: &[u32],
    m: u32,
    probabilities: &[f64],
    max_lengths: &[u32],
    lengths: &mut [u32],
) -> Result<(), good_lp::ResolutionError> {
    let n = probabilities.len() as u32;

    // Enthält alle Variablen und ihre Konfigurationen
    let mut problem = ProblemVariables::new();
    
    // Enthält Referenzen zu den Variablen
    let mut variables = vec![vec![Option::<Variable>::None; m as usize]; n as usize];

    for i in 0..n {
        for j in 0..(max_lengths[i as usize].min(m)) {
            variables[i as usize][j as usize] =
                Some(problem.add(VariableDefinition::new().binary()));
        }
    }

    // Der Ausdruck, der minimiert sein muss
    let mut minimise = Expression::default();
    for i in 0..n {
        for j in 0..m {
            if let Some(var) = variables[i as usize][j as usize] {
                minimise = minimise + probabilities[i as usize] * ((j + 1) as f64) * var;
            }
        }
    }

    let mut problem = problem
        .minimise(minimise)
        .using(good_lp::default_solver);

    // Die erste Einschänkungen
    for i in 1..m {
        let mut expr = Expression::default();
        for p in 1..=i {
            for j in 1..=n {
                if let Some(var) = variables[j as usize - 1][p as usize - 1] {
                    expr = expr + var * v.get((i - p) as usize);
                }
            }
        }
        problem = problem.with(expr.leq(v.get(i as usize)));
    }

    // Die zweite Einschränkungen
    for i in 0..n {
        let mut expr = Expression::default();
        for j in 0..m {
            if let Some(var) = variables[i as usize][j as usize] {
                expr = expr + var;
            }
        }
        problem = problem.with(expr.eq(1));
    }

    let res = problem.solve()?;

    for i in 0..n {
        for j in 0..m {
            // >= 0.5, weil die Vergleichoperationen bei Floats nicht sehr sicher sind.
            if res.value(variables[i as usize][j as usize].unwrap()) >= 0.5 {
                lengths[i as usize] = j + 1;
                break;
            }
        }
    }

    lengths.sort_unstable();

    Ok(())
}
