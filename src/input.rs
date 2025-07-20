use std::collections::HashMap;

use crate::scanner::Scanner;

pub struct Input {
    pub letters: Vec<u32>,
    pub probabilities: Vec<f64>,
    pub symbols: Vec<(char, u64)>,
    pub sum: f64,
}

impl Input {
    pub fn read(message: bool, mut mul: f64, scanner: &mut Scanner<impl std::io::BufRead>) -> Self {
        let r = scanner.read::<usize>();
        let mut letters = vec![];
        for _ in 0..r {
            letters.push(scanner.read());
        }

        letters.sort_unstable();

        let mut symbols = vec![];
        let mut probabilities = vec![];

        // Es gibt zwei Möglichkeiten vom Dateiformat:
        // - eine reine Nachricht
        // - eine Liste von Wahrscheinlichkeiten
        if message {
            let msg = scanner.read_line();
            let mut map = HashMap::new();
            for c in msg.chars() {
                *map.entry(c).or_insert(0) += 1;
            }
            symbols = map.into_iter().collect::<Vec<_>>();
        } else {
            for _ in 0..scanner.read::<usize>() {
                let c = scanner.read::<char>();
                let p = scanner.read::<u64>();
                symbols.push((c, p));
            }
        }
        let sum = symbols.iter().map(|&(_, v)| v).sum::<u64>() as f64;
        // sortieren in absteigender Reihenfolge (größere Werte zuerst)
        symbols.sort_unstable_by_key(|&(_, v)| u64::MAX - v);
        // wenn mul gleich sum ist, dann werden alle Wahrscheinlichkeiten gleich 
        // die Anzahl der Symbole in der Nachricht (bzw. die von Datei gegebene Wahrscheinlichkeiten)
        if mul == 0.0 {
            mul = sum;
        }
        for &(_, v) in symbols.iter() {
            probabilities.push(v as f64 * (mul / sum));
        }

        Self {
            letters,
            probabilities,
            symbols,
            sum,
        }
    }
}
