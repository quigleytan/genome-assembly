# M4OEP-tequigle: Gene Reconstruction and Assembly Toolkit CURRENT
Author: Tanner Quigley
### About
This project is a C++ implementation of a bioinformatics toolkit to analyze DNA sequences, reconstruct linear and 
circular genomes, and visualize genome assembly.
### Features
- DNA sequence representation and k-mer analysis.
- Performance optimization for large datasets.
- Implementation of custom data structures.
- Eulerian path computation and visualization.
- K-mer, contig, and scaffold construction and visualization.

## Brief Overview
***
### Module 1: DNA Sequence Representation and K-mer Analysis
Implementation of data structures to hold, analyze, and conduct basic analysis on genomic sequences.
- **Classes**:
- `DNASequence`: Represents a DNA sequence with methods for validation and manipulation.
- `KmerEncoder`: Encodes k-mers for efficient storage and retrieval.

### Module 2: DNA Sequence Handling and File I/O
Allow proper file reading and input of data from FASTA files.
- `SequenceReader`: Reads DNA sequences from FASTA files.
- `OpenAddressingTable`: Linear probing hash table implementation.
- `KmerTable`: Child class of `OpenAddressingTable` efficient k-mer storage and retrieval.
- `DeBruijnGraph`: Writes DNA sequences to FASTA files, has an `OpenAddressingTable` for node storage. Does not compress
edges, as this is crucial from Eulerian walk assembly.
- `InitializationTests`, `ProcessingTests`, `AssemblyTests`: Tests proper logical data manipulation, file I/O, and
implementation of data structures.

### Module 3: Genome Reconstruction and Assembly
Eulerian walk implementation for genome assembly.

### Module 4:
Visualization of DNA sequences and assembly results using C++ graphics.

## Full Overview and Report

***

### Module 1
This module focused on creating classes and laying out the groundwork for future scaling of my gene sequencer. The main
classes created were DNASequence, KmerCounter, and KmerEncoder. Each class has its own header and implementation files,
with proper documentation for each method and member variable. These files allow for storage of DNA sequences and its
important information, which will be important when simulating assembly. At the time of submission, I have no
significant bugs to report from my testing of the file. 

I think that for this module, I earned approximately 50-60 points, as I fully implemented several classes split into
header and cpp files. Additionally, while I only included one example of enumerated types and overloaded operators, I
included them in the places where I felt they were necessary only. I plan to scale this project and use it for future
modules, and I didn't want to overcomplicate the code with unnecessary features. My program's main function is for this
stage is complete, and it serves as the middle ground between reading file information and producing a final output of
an assembled genome.

***

### Module 2

This module was focused around reading in FASTA files and implementation of a De Bruijn graph. File I/O will allow for
testing with much larger datasets and better program flow, as before I would prompt for a short sequence from the user.
The De Bruijn graph implementation lays out the path for a future Eulerian walk algorithm. This algorithm is the actual
assembly step, as it will trace a non-repeating path through the graph. 

Regarding the changes I made, I realized that I needed a data structure to store the nodes of the graph. To do this, I
made a parent class from `KmerTable` called `OpenAddressingTable` that stores items in a hash table. This was 
challenging as`KmerTable` was very specialized, as it was designed to store k-mer counts without having to manually
access data from outside the insert() function. To solve this, I generalized the insert() method to return the item, 
which allows me to update node information for `DeBruijnGraph`. Now that I couldn't just increment value, I added
virtual functions in my insert() method to allow for different behavior on duplicate and unique insertion cases when
called from `KmerTable`.

I think from this module I have earned approximately 100 points, as I implemented both examples of IS A and HAS A in
my classes for hash table functions. Additionally, I incorporated file I/O to read in sequences from a FASTA
(bioinformatics formatted text) file. Finally, I built three testing files to ensure that my encoding, information
storage, file reading, and graph logic all work as intended. They are split up into the stages of genome assembly and
are intended to be run in the order they are listed in the Module 2 quick overview, as each stage is reliant on the ones
before it. For more information regarding my testing files and outputs, see the testing section of the `README.md`.

### Testing Summary
I wanted to explain my testing cases, especially for `AssemblyTests.cpp`, as the values being tested for can seem a bit
arbitrary. Starting with `InitializationTests.cpp`, I am just checking to make sure that the information being read in
from the FASTA file is correct, and that basic DNA sequence information is correct. Next, `ProcessingTests.cpp` is
slightly more involved, as I am using 2-bit encoding logic to reduce the space complexity of my project. Aditionally,
it is making sure that the hooks I have built into my child class, KmerTable, are working and ensuring that the k-mer
counts are being incremented properly. Finally, `AssemblyTests.cpp` is the most complex, as we are testing proper node
and transition logic for a De Bruijn graph. Firstly, I am ensuring that the proper number of nodes and transitions are
being created. Since our sequence is "AGTGCGTCAGT" and we use a k value of 3, we find the k-mers:

| 3-mer | Prefix | Suffix |
|------:|:------:|:------:|
| AGT   | AG     | GT     |
| GTG   | GT     | TG     |
| TGC   | TG     | GC     |
| GCG   | GC     | CG     |
| CGT   | CG     | GT     |
| GTC   | GT     | TC     |
| TCA   | TC     | CA     |
| CAG   | CA     | AG     |
| AGT   | AG     | GT     |

Note: Having 9 k-mers makes sense thanks to the formula: |sequence| − k + 1 = 11 − 3 + 1 = 9

This gives us 7 unique 2-mers, giving us a node count of 7. Additionally, if we track the number of transitions between
these nodes, we find that there are 9, which is our edge count. Additionally, using AG as an example, we see that there
are two prefix occurrences and one suffix occurrence. This gives us our indegree and outdegree values. In the future, I
will include more robust testing for checking expected neighbors and transitions. As of now, the tests I have run
through `main.cpp` align with expected outputs regarding possible neighbors, but I still need to test for edge cases.

### Module 3

Whilst I did not submit anything for this module, I still wanted to write this section to give a bit of background and
explanation to the processes and file seen in the jump from module 2 to 4. The freedom of disregarding module
requirements allowed me to make significant progress on the assembly algorithms/pipelines. I implemented a simple
Eulerian walk algorithm that traces a non-repeating path through the graph. This was able to reconstruct simple genomes 
given that the size of repeats in the sequence was below a certain k threshold. I expanded my handling of k-values up
to 63 by using `__uint128_t`. I was able to reconstruct the genome of Escherichia phage phiX174, which consists of
roughly 5000 bases. However, as predicted, this algorithm was not able to reconstruct the genome of Escherichia coli K12
due to the sheer size of the repeats in the sequence. Whilst Eulerian walks are able to process repeats, it requires
the k value used to exceed the length of said repeats. In a genome with a very large number of repeats such as 
Escherichia coli, the algorithm constantly takes the wrong path, resulting in an assembled but incorrect genome. 
`EulerianTraversal.cpp` effectively served as a proof of concept for sequence reconstruction, and gave me an outline for
`ContigTraversal.cpp` and `ContigScaffolder.cpp`.

By using a new pipeline using multi-stage traversal which included different construction strategies that help the
algorithm pick which path is most likely to be correct, greatly increase the chance of a correct assembly. The contigs
are generated from the k-mers in the graph by tracing through looking for overlaps and k-mer adjacency, similar to the
original algorithm. Once the contigs are generated, they are ordered and assembled into scaffolds, which are then
assembled into the final genome.

### Module 4




Visualization piepline
Step 1 — CMakeLists.txt changes (FetchContent + Visualizer executable)
Step 2 — VisData.h (all shared structs, no dependencies)
Step 3 — Recorder.h (header-only step recorder, optional ptr in traversal classes)
Step 4 — DataExporter.h/.cpp (writes .visdata from VisSession)
Step 5 — DataLoader.h/.cpp (reads .visdata back into VisSession)
Step 6 — Modify ContigTraversal.cpp (~10 lines, non-breaking)
Step 7 — Modify ContigScaffolder.cpp (~10 lines, non-breaking)
Step 8 — ContigView.h/.cpp (bar renderer + animation)
Step 9 — VisualizerApp.cpp (window, loop, top-level state)





### References

| Module / Area       | Topic                          | Link                                                                 |
|--------------------|--------------------------------|----------------------------------------------------------------------|
| General            | Usage of static                | https://www.geeksforgeeks.org/cpp/static-keyword-cpp/               |
| General            | Documentation guidelines       | https://developer.lsst.io/cpp/api-docs.html                         |
| General            | Usage of auto                  | https://www.geeksforgeeks.org/cpp/type-inference-in-c-auto-and-decltype/ |
| General            | Exceptions                     | https://www.geeksforgeeks.org/cpp/how-to-throw-custom-exception-in-cpp/ |
| DNASequence        | Helper functions               | https://www.w3tutorials.net/blog/what-are-helper-functions-in-c/     |
| DNASequence        | size_t                         | https://www.geeksforgeeks.org/cpp/difference-between-int-and-size_t-in-cpp/ |
| DNASequence        | Switch and cases               | https://www.w3schools.com/cpp/cpp_switch.asp                        |
| KmerEncoder        | Bitshift and masking           | https://www.geeksforgeeks.org/cpp/left-shift-right-shift-operators-c-cpp/ |
| OpenAddressingTable| Iterator implementation        | https://stackoverflow.com/questions/46431762/how-to-implement-standard-iterators-in-class |
| DeBruijnGraph      | Usage of auto                  | Reference source from Module 1                                       |
| DeBruijnGraph      | Use of [[nodiscard]]           | https://stackoverflow.com/questions/76489630/explanation-of-nodiscard-in-c17 |
| DeBruijnGraph      | Use of explicit                | https://www.geeksforgeeks.org/cpp/use-of-explicit-keyword-in-cpp/    |