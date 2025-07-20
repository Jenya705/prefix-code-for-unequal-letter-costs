use std::{
    cmp::Reverse,
    collections::{BinaryHeap, HashMap},
};

/// In der Sektion "Finden eines zufälligen Codes" beschriebene Funktion.
/// adjusted muss q Werten enthalten.
/// Die Funktion schreibt dann in die Liste die Ergebniswerte.
pub fn adjust_costs(letters: &[u32], adjusted: &mut [u32]) {
    let mut heap = BinaryHeap::new();
    heap.push(Reverse(0u32));
    for i in 0..adjusted.len() {
        loop {
            let node = heap.pop().unwrap().0;
            if node < adjusted[i] || (heap.is_empty() && i != adjusted.len() - 1) {
                for k in 0..letters.len() {
                    heap.push(Reverse(letters[k] + node));
                }
            } else {
                adjusted[i] = node;
                break;
            }
        }
    }
}

pub fn build_prefixes(letters: &[u32], costs: &[u32], prefixes: &mut [Vec<usize>]) {
    /// Implementiert Vergleichfunktionen,
    /// so dass die Liste ignoriert wird und der binäre Heap Min Heap war
    struct InHeap(u32, Vec<usize>);

    impl PartialEq for InHeap {
        fn eq(&self, other: &Self) -> bool {
            self.0 == other.0
        }
    }

    impl PartialOrd for InHeap {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            other.0.partial_cmp(&self.0)
        }
    }

    impl Eq for InHeap {}

    impl Ord for InHeap {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            self.partial_cmp(other).unwrap()
        }
    }

    // Ein Max-Heap
    let mut heap = BinaryHeap::new();
    heap.push(InHeap(0u32, vec![0usize; 0]));
    for i in 0..costs.len() {
        loop {
            let InHeap(node, prefix) = heap.pop().unwrap();
            // nicht notwendig
            assert!(node <= costs[i]);
            if node == costs[i] {
                prefixes[i] = prefix;
                break;
            } else {
                // der Knoten gehört keinem Symbol, also muss er erweitert werden
                for i in 0..letters.len() {
                    let mut prefix = prefix.clone();
                    prefix.push(i);
                    heap.push(InHeap(letters[i] + node, prefix));
                }
            }
        }
    }
}

/// Die Funktion optimiert die Präfixe, damit sie die verkürzt, wenn es möglich ist
#[allow(unused)]
pub fn optimize_prefixes(prefixes: &mut [Vec<usize>]) {
    let mut map = HashMap::new();

    for prefix in prefixes.iter() {
        for i in 1..prefix.len() {
            let key = &prefix[..i];
            if map.contains_key(key) {
                *map.get_mut(key).unwrap() += 1;
            } else {
                map.insert(key.to_vec(), 1);
            }
        }
    }

    for prefix in prefixes {
        loop {
            let l = prefix.len() - 1;
            if matches!(map.get(&prefix[0..l]), Some(&1)) {
                prefix.remove(l);
            } else {
                break;
            }
        }
    }
}

/// Kalkuliert Kosten des PC mit den gegebenen Präfixen
#[allow(unused)]
pub fn calculate_costs(letters: &[u32], prefixes: &[Vec<usize>], costs: &mut [u32]) {
    for i in 0..prefixes.len() {
        costs[i] = 0;
        for &letter in &prefixes[i] {
            costs[i] += letters[letter];
        }
    }
}
