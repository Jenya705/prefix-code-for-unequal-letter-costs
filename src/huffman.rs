use std::{cmp::Ordering, collections::BinaryHeap, rc::Rc};

/// Gibt zurück, ob der Huffman-Algorithmus für den Fall optimal wird.
pub fn huffman_optimal(letters: &[u32]) -> bool {
    letters.iter().all(|&v| v == letters[0])
}

/// Eine Knotenstructur
#[derive(Debug)]
struct Node {
    probability: f64,
    parent: Option<Rc<Node>>,
    /// Länge (bzw. Kosten) zwischen dem Elternknoten und diesem Knoten  
    length: u32,
}

impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.probability == other.probability
    }
}

impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        other.probability.partial_cmp(&self.probability)
    }
}

impl Eq for Node {}

impl Ord for Node {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

pub fn huffman_solve(letters: &[u32], probabilities: &[f64], lengths: &mut [u32]) {
    let n = probabilities.len() % (letters.len() - 1);
    let n = if n == 0 { letters.len() - 1 } else { n };

    let mut heap = BinaryHeap::new();
    // Alle Blätter
    let mut original = Vec::with_capacity(probabilities.len());

    for &probability in probabilities.iter() {
        let node = Rc::new(Node {
            probability,
            parent: None,
            length: 0,
        });

        original.push(Rc::clone(&node));
        heap.push(node);
    }

    // l ist die Anzahl der Kinder
    let mut create_node = |l: usize| -> bool {
        // root node
        // am Ende bleibt der Kernknoten immer
        if heap.len() == 1 {
            return false;
        }

        let mut parent = Rc::new(Node {
            probability: 0.0,
            parent: None,
            length: 0,
        });

        for i in (0..l).rev() {
            let mut node = heap.pop().unwrap();
            // SAFETY:
            // - no lifetimes
            // - no raw references to any rc exist
            unsafe { Rc::get_mut_unchecked(&mut parent) }.probability += node.probability;
            unsafe { Rc::get_mut_unchecked(&mut node) }.length = letters[i];
            unsafe { Rc::get_mut_unchecked(&mut node) }.parent = Some(Rc::clone(&parent));
        }

        heap.push(parent);

        true
    };

    // Zuerst ein Knoten mit n Kindern
    if n != 1 {
        create_node(n);
    }
    // Danach mit r Kindern
    while create_node(letters.len()) {}

    // Kalkulieren der a_i Werte
    for (i, node) in original.iter().enumerate() {
        let mut current = node;
        // Gehen vom Blatt zum Kernknoten
        while let Some(ref parent) = current.parent {
            lengths[i] += current.length;
            current = parent;
        }
    }
}
