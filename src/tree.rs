use std::{
    cmp::Reverse,
    collections::{BinaryHeap, HashMap},
};

pub type TreeIndex = u128;

pub trait TreeSolver {
    fn new(letters: &[u32]) -> Self;

    fn solve(&mut self, lengths: &mut [u32], indices: &mut [TreeIndex]) -> Option<bool>;
}

pub struct FirstToFindTreeSolver {
    letters: Vec<u32>,
    heap: BinaryHeap<(Reverse<u32>, u128)>,
}

impl TreeSolver for FirstToFindTreeSolver {
    fn new(letters: &[u32]) -> Self {
        Self {
            letters: letters.to_vec(),
            heap: BinaryHeap::new(),
        }
    }

    fn solve(&mut self, mut lengths: &mut [u32], mut indices: &mut [TreeIndex]) -> Option<bool> {
        let mut length_changed = false;

        self.heap.clear();
        self.heap.push((Reverse(0), 0));

        while let Some((Reverse(cost), mut i)) = self.heap.pop() {
            if cost < lengths[0] {
                i *= self.letters.len() as TreeIndex;
                for &letter in self.letters.iter() {
                    i += 1;
                    self.heap.push((Reverse(cost + letter), i));
                }
            } else {
                length_changed = length_changed || lengths[0] != cost;
                indices[0] = i;
                lengths[0] = cost;
                lengths = &mut lengths[1..];
                indices = &mut indices[1..];
                if lengths.len() == 0 {
                    break;
                }
            }
        }

        (lengths.len() == 0).then_some(length_changed)
    }
}

pub struct TreeOptimizer {
    map: HashMap<TreeIndex, u32>,
    vec: Vec<usize>,
    l_copy: Vec<u32>,
    i_copy: Vec<TreeIndex>,
}

impl TreeOptimizer {
    pub fn new(occurences: usize) -> Self {
        Self {
            map: HashMap::new(),
            vec: vec![0; occurences],
            l_copy: vec![0; occurences],
            i_copy: vec![0; occurences],
        }
    }

    pub fn optimize(
        &mut self,
        letters: &[u32],
        lengths: &mut [u32],
        indices: &mut [TreeIndex],
    ) -> bool {
        self.map.clear();

        let mut optimized = false;

        for i in indices.iter() {
            let mut i = *i;
            loop {
                i -= 1;
                i /= letters.len() as TreeIndex;
                if i == 0 {
                    break;
                }
                *self.map.entry(i).or_insert(0) += 1;
            }
        }

        for ii in 0..indices.len() {
            let mut i = indices[ii];
            loop {
                let mut ni = i;
                ni -= 1;
                ni /= letters.len() as TreeIndex;
                if ni == 0 {
                    break;
                }
                if self.map.get(&ni) != Some(&1) {
                    break;
                }
                i = ni;
            }
            if i != indices[ii] {
                indices[ii] = i;
                let mut new_length = 0u32;
                loop {
                    i -= 1;
                    new_length += letters[(i % letters.len() as TreeIndex) as usize];
                    i /= letters.len() as TreeIndex;
                    if i == 0 {
                        break;
                    }
                }
                lengths[ii] = new_length;
                optimized = true;
            }
        }

        if optimized {
            for i in 0..self.vec.len() {
                self.vec[i] = i;
            }

            self.vec.sort_unstable_by_key(|&v| lengths[v]);

            self.l_copy.copy_from_slice(&lengths);
            self.i_copy.copy_from_slice(&indices);

            for i in 0..self.vec.len() {
                lengths[i] = self.l_copy[self.vec[i]];
                indices[i] = self.i_copy[self.vec[i]];
            }
        }

        optimized
    }
}
