use good_lp::{
    default_solver, solvers::highs::HighsParallelType, variable, variables, Expression, Solution,
    SolverModel, VariableDefinition,
};

pub trait LengthSolver {
    const BEST_SOLVER: bool;

    fn new(occurences: &[u32], letters: &[u32]) -> Self;

    fn m(&self) -> usize;

    fn theoretical_max(&self) -> f64;

    fn solve(&mut self, max_cost: f64, lengths: &mut [u32]) -> Option<f64>;
}

pub struct LPSolver {
    occurences: Vec<f64>,
    letters: Vec<f64>,
    m: usize,
    v_cacher: VCacher,
    previous_attempts: Vec<Vec<u32>>,
}

impl LPSolver {
    fn new(occurences: &[u32], letters: &[u32]) -> Self {
        let mut v_cacher = VCacher::new();

        let letters = letters.iter().map(|v| *v as f64).collect::<Vec<_>>();

        Self {
            m: (v_cacher.determine_max_h(100, &letters) as usize).min(100),
            occurences: occurences.iter().map(|v| *v as f64).collect(),
            letters,
            v_cacher,
            previous_attempts: vec![],
        }
    }

    fn solve<F>(&mut self, lengths: &mut [u32], mut configurate_variable: F) -> Option<f64>
    where
        F: FnMut(VariableDefinition) -> VariableDefinition,
    {
        loop {
            let mut variables = variables!();

            let mut vars_mat = vec![];

            for _ in 0..self.occurences.len() {
                let mut v = vec![];
                for _ in 0..self.m {
                    v.push(variables.add((configurate_variable)(variable())));
                }
                vars_mat.push(v);
            }

            let mut minimise = Expression::default();

            for i in 0..self.occurences.len() {
                for j in 0..self.m {
                    minimise = minimise + vars_mat[i][j] * (self.occurences[i] * (j + 1) as f64);
                }
            }

            let mut problem = variables.minimise(minimise.clone()).using(default_solver);

            for j in 0..self.m {
                let mut expr = Expression::default();
                for p in 0..=j {
                    for i in 0..self.occurences.len() {
                        expr = expr
                            + vars_mat[i][p]
                                * self.v_cacher.get_f64(j as i32 - p as i32, &self.letters);
                    }
                }
                problem =
                    problem.with(expr.leq(self.v_cacher.get_f64(j as i32 + 1, &self.letters)));
            }

            for i in 0..self.occurences.len() {
                let mut expr = Expression::default();
                for j in 0..self.m {
                    expr = expr + vars_mat[i][j];
                }
                problem = problem.with(expr.eq(1));
            }

            for prev in self.previous_attempts.iter() {
                let mut expr = Expression::default();
                for i in 0..self.occurences.len() {
                    expr = expr + vars_mat[i][prev[i] as usize - 1];
                }
                problem = problem.with(expr.leq(self.occurences.len() as f64 - 1.0));
            }

            // making the lengths sorted

            if !self.previous_attempts.is_empty() {
                let gen_expr = |i: usize| -> Expression {
                    let mut expr = Expression::default();
                    for j in 0..self.m {
                        expr = expr + j as f64 * vars_mat[i][j];
                    }
                    expr
                };

                for i in 1..self.occurences.len() {
                    problem = problem.with(gen_expr(i - 1).leq(gen_expr(i)));
                }
            }

            problem = problem.set_parallel(HighsParallelType::On).set_threads(8);

            // I cannot prove that generation of max m is right, thus
            // I need to handle the case when it is not and try to adjust the value
            // in order to be able to solve the problem in the first place
            // but the library's author didn't make the error handling for cases, when the problem
            // itself is defined wrongly, so I have to do it myself :P
            let solution = match std::panic::catch_unwind(|| problem.solve()) {
                Ok(sol) => sol.ok()?,
                Err(_) => {
                    self.m -= 1;
                    continue;
                }
            };

            let calculate_l = |i: usize| -> u32 {
                for j in 0..self.m {
                    if solution.value(vars_mat[i][j]) > 0.5 {
                        return j as u32 + 1;
                    }
                }
                unreachable!()
            };

            for i in 0..self.occurences.len() {
                lengths[i] = calculate_l(i);
            }

            if self.previous_attempts.is_empty() {
                lengths.sort_unstable();
            }
            
            self.previous_attempts.push(lengths.to_vec());

            println!("{:?}", lengths);

            println!("{}", solution.eval(minimise));

            let cost = lengths
                .iter()
                .zip(self.occurences.iter())
                .map(|(&l, &o)| o * l as f64)
                .sum();

            return Some(dbg!(cost));
        }
    }
}

pub struct ILPSolver {
    solver: LPSolver,
}

pub struct VCacher {
    cache: Vec<u64>,
}

impl VCacher {
    pub fn new() -> Self {
        Self { cache: vec![1] }
    }

    // TODO: make it work for ilp
    pub fn determine_max_h(&mut self, max: i32, letters: &[f64]) -> i32 {
        let mut sum = 0u64;

        for h in 0..max {
            if self.try_get(h, letters).is_none() {
                return h;
            }

            match sum
                .checked_add(self.get(h, letters))
                .filter(|&v| v < u32::MAX as u64)
            {
                Some(v) => sum = v,
                None => return h,
            }
        }
        max
    }

    pub fn try_get(&mut self, h: i32, letters: &[f64]) -> Option<u64> {
        let res = if h < 0 {
            0
        } else if self.cache.len() as i32 > h {
            self.cache[h as usize]
        } else {
            for h in self.cache.len() as i32..h {
                let _ = self.try_get(h, letters)?;
            }
            let mut v = 0u64;
            for i in 0..letters.len() {
                v = v.checked_add(self.try_get(h - letters[i] as i32, letters)?)?;
            }
            self.cache.push(v);
            v
        };

        Some(res)
    }

    pub fn get(&mut self, h: i32, letters: &[f64]) -> u64 {
        self.try_get(h, letters).unwrap()
    }

    pub fn get_f64(&mut self, h: i32, letters: &[f64]) -> f64 {
        let res = self.get(h, letters);
        debug_assert!(res < f64::MAX as u64);
        res as f64
    }
}

impl LengthSolver for ILPSolver {
    const BEST_SOLVER: bool = true;

    fn new(occurences: &[u32], letters: &[u32]) -> Self {
        Self {
            solver: LPSolver::new(occurences, letters),
        }
    }

    fn m(&self) -> usize {
        self.solver.m
    }

    fn solve(&mut self, _max_cost: f64, lengths: &mut [u32]) -> Option<f64> {
        self.solver.solve(lengths, |variable| variable.binary())
    }

    fn theoretical_max(&self) -> f64 {
        // we do not need it here, because it always returns the BEST result possible and doesn't use
        // the max cost anyway
        f64::MAX
    }
}

pub struct RelaxedLPSolver {
    solver: LPSolver,
}

impl LengthSolver for RelaxedLPSolver {
    const BEST_SOLVER: bool = false;

    fn new(occurences: &[u32], letters: &[u32]) -> Self {
        Self {
            solver: LPSolver::new(occurences, letters),
        }
    }

    fn m(&self) -> usize {
        self.solver.m
    }

    fn theoretical_max(&self) -> f64 {
        f64::MAX
    }

    fn solve(&mut self, _max_cost: f64, lengths: &mut [u32]) -> Option<f64> {
        self.solver
            .solve(lengths, |variable| variable.bounds(0..=1))
    }
}

pub struct ESolver {
    occurences: Vec<f64>,
    occurences_prefix_sums: Vec<f64>,
    letters: Vec<f64>,
    t_powers: Vec<f64>,
    t_powers_prefix_sums: Vec<f64>,
    i0: usize,
    j0: usize,
    m: usize,
    goal: f64,
    order: Vec<(usize, usize)>,
}

fn find_t(letters: &[f64]) -> f64 {
    // using newton's method
    let mut t0 = 1.0f64;

    for _ in 0..10 {
        let mut f_t0 = -1.0;
        let mut fdt_t0 = 0.0;

        for &letter in letters {
            f_t0 += t0.powi(letter as i32);
            // derivative of the power function
            fdt_t0 += letter as f64 * t0.powi(letter as i32 - 1);
        }

        // step
        t0 -= f_t0 / fdt_t0;
    }

    t0
}

impl ESolver {
    fn progress(&mut self) -> bool {
        self.i0 += 1;
        if self.i0 == self.occurences.len() {
            self.j0 += 1;
            self.i0 = 0;
            if self.j0 == self.m {
                return false;
            }
        }
        true
    }
}

impl LengthSolver for ESolver {
    const BEST_SOLVER: bool = false;

    fn new(occurences: &[u32], letters: &[u32]) -> Self {
        let occurences = occurences.iter().map(|v| *v as f64).collect::<Vec<_>>();
        let letters = letters.iter().map(|v| *v as f64).collect::<Vec<_>>();

        let m = VCacher::new().determine_max_h(100, &letters) as usize;

        let t = find_t(&letters);
        let mut t_powers = vec![1.0; m];
        let mut t_powers_prefix_sums = vec![1.0; m];
        for i in 1..m {
            t_powers[i] = t_powers[i - 1] * t;
            t_powers_prefix_sums[i] = t_powers_prefix_sums[i - 1] + t_powers[i];
        }

        let mut occurences_prefix_sums = vec![occurences[0]; occurences.len()];
        for i in 1..occurences.len() {
            occurences_prefix_sums[i] = occurences_prefix_sums[i - 1] + occurences[i];
        }

        let goal = (occurences.len() as f64 - 1.0) / (1.0 - t);

        let mut sorted = Vec::with_capacity(m * occurences.len());
        for i in 0..occurences.len() {
            for j in 0..m {
                sorted.push((i, j, occurences[i] * t.powi(-(j as i32 + 1))));
            }
        }
        sorted.sort_unstable_by(|v1, v2| {
            v1.2.partial_cmp(&v2.2)
                .unwrap()
                .then(v2.0.cmp(&v1.0))
                .then(v2.1.cmp(&v1.1))
        });
        let sorted = sorted.into_iter().map(|v| (v.0, v.1)).collect::<Vec<_>>();

        Self {
            occurences,
            occurences_prefix_sums,
            letters,
            t_powers,
            t_powers_prefix_sums,
            i0: 0,
            j0: 0,
            m,
            goal,
            order: sorted,
        }
    }

    fn m(&self) -> usize {
        self.m
    }

    fn theoretical_max(&self) -> f64 {
        const DIV: f64 = 0.90;

        let sub = self
            .occurences_prefix_sums
            .last()
            .unwrap()
            .log(self.t_powers[1]);

        self.occurences
            .iter()
            .map(|&v| v * (v.log(self.t_powers[1]) - sub))
            .sum::<f64>()
            / DIV
    }

    fn solve(&mut self, max_cost: f64, lengths: &mut [u32]) -> Option<f64> {
        if self.j0 == self.m {
            // exhausted
            return None;
        }

        loop {
            lengths.fill(0);

            let mut res = 0.0;

            // setting initial lengths
            for i in self.i0..self.occurences.len() {
                lengths[i] = self.j0 as u32 + 1;
                // TODO: occurences sums
                // res += self.occurences[i];
            }
            res += self.occurences_prefix_sums[self.occurences.len() - 1];
            res *= (self.j0 + 1) as f64;

            // calculating initial subject
            let mut subject = 0.0;
            /*
            for j in 0..=self.j0 {
                // TODO: t_powers sums
                subject += self.t_powers[j];
            }
            */
            subject += self.t_powers_prefix_sums[self.j0];
            subject *= (self.occurences.len() - self.i0) as f64;

            for &(i, j) in self.order.iter().filter(|&&(i, j)| {
                !(i <= self.i0 && j > self.j0) && !(i >= self.i0 && j <= self.j0)
            }) {
                subject += self.t_powers[j];
                lengths[i] += 1;
                res += self.occurences[i];

                if res > max_cost {
                    break;
                }

                if subject >= self.goal {
                    self.progress();
                    return Some(res);
                }
            }

            if !self.progress() {
                return None;
            }
        }
    }
}

// A hacky way to generate trees for a given input
pub struct GivenSolver {
    result: Option<Vec<u32>>,
}

impl LengthSolver for GivenSolver {
    const BEST_SOLVER: bool = true;

    fn new(occurences: &[u32], _letters: &[u32]) -> Self {
        Self {
            result: Some(occurences.to_vec()),
        }
    }

    fn m(&self) -> usize {
        self.result
            .as_ref()
            .and_then(|v| v.last().cloned())
            .unwrap_or(0) as usize
    }

    fn theoretical_max(&self) -> f64 {
        f64::MAX
    }

    fn solve(&mut self, _max_cost: f64, lengths: &mut [u32]) -> Option<f64> {
        let result = self.result.take()?;

        lengths.clone_from_slice(&result);

        // cost doesn't matter at this point
        Some(1.0)
    }
}
