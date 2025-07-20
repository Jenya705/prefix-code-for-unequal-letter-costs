An implementation of an algorithm proposed by Richard Karp that generates an optimal huffman tree for unequal letter costs. The program also includes the original optimal for equal letter costs algorithm by Huffman, which is unfortunately not optimal for unequal letter costs.
# Help page
Usage: bwinf4321 [OPTIONS] <FILE_NAME>

Arguments:
  <FILE_NAME>
          

Options:
  -o, --occurences
          The other format for input

  -m, --mul <MUL>
          The factor on which each probability will be multiplied (in some cases it can lead to better results)
          
          If the value is 0, then instead of probabilities occurences will be used.
          
          [default: 0]

  -u, --huffman
          Whether the Huffman-Algorithm will be used

  -h, --help
          Print help (see a summary with '-h')
