# M4OEP-tequigle: Gene Reconstruction and Assembly Toolkit
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

## Brief Overview and Key Files
***
### Module 1: DNA Sequence Representation and K-mer Analysis
Implementation of data structures to hold, analyze, and conduct basic analysis on genomic sequences.
- **Classes**:
- `dna_sequence`: Represents a DNA sequence with methods for validation and manipulation.
- `kmer_encoding`: Encodes k-mers for efficient storage and retrieval.

### Module 2: DNA Sequence Handling and File I/O
Allow proper file reading and input of data from FASTA files.
- `sequence_reader`: Reads DNA sequences from FASTA files.
- `open_addressing_table`: Linear probing hash table implementation.
- `kmer_table`: Child class of `open_addressing_table` efficient k-mer storage and retrieval.
- `de_bruijn_graph`: Writes DNA sequences to FASTA files, has an `open_addressing_table` for node storage. Does not compress
edges, as this is crucial from Eulerian walk assembly.
- `InitializationTests`, `ProcessingTests`, `AssemblyTests`: Tests proper logical data manipulation, file I/O, and
implementation of data structures.

### Module 3: Genome Reconstruction and Assembly
Implementation of main assembly algorithms and pipelines.
- `eulerian_traversal`: Implements a simple Eulerian walk algorithm to reconstruct genomes from De Bruijn graphs.
- `contig_traversal`: Implements a multi-stage traversal strategy to generate contigs from De Bruijn graphs.
- `contig_scaffolder`: Orders and assembles contigs into scaffolds, and scaffolds into the final genome assembly.

### Module 4:
Visualization of DNA sequences and assembly results using C++ graphics.
- `contig_view`: Contains rendering logic for contigs and scaffolds.
- `vis_loader`, `vis_exporter`: Handles loading and exporting of visualization data.
- `recorder`: Allows for recording of assembly results to a visdata file.
- `gui`: Main visualization program.

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
which allows me to update node information for `de_bruijn_graph`. Now that I couldn't just increment value, I added
virtual functions in my insert() method to allow for different behavior on duplicate and unique insertion cases when
called from `KmerTable`.

I think from this module I have earned approximately 100 points, as I implemented both examples of IS A and HAS A in
my classes for hash table functions. Additionally, I incorporated file I/O to read in sequences from a FASTA
(bioinformatics formatted text) file. Finally, I built three testing files to ensure that my encoding, information
storage, file reading, and graph logic all work as intended. They are split up into the stages of genome assembly and
are intended to be run in the order they are listed in the Module 2 quick overview, as each stage is reliant on the ones
before it. For more information regarding my testing files and outputs, see the testing section of the `README.md`.

### Testing Summary (Module 2)
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

NOTE: Any files mentioned in the modules prior have likely been renamed and reorganized, explanations below.

For this module, I focused on creating a simple yet independent visualization of the assembly process. Up until this
point, I had been printing out stdout-based visualizations of the assembly process built into the main pipeline. 
However, I want to separate the visualization from the assembly process, as it is not necessary for the assembly to run,
and it allows for quicker testing of the algorithmic logic and overall reconstruction accuracy. This is especially
important if the pipeline was to open up a graphics window on every run. From this, I currently have two main files:
`scaffold_assembly.cpp` and `gui.cpp`. `scaffold_assembly.cpp` is the main file that runs the assembly pipeline and
creates informative files regarding the results of the assembly. This also creates the visdata files that my graphics
pipeline reads in. `gui.cpp` is the file that creates the graphical user interface and currently contains two tabs once
the user provides the proper visdata file. The first tab shows each scaffold and its corresponding contigs, while
the second tab shows the final assembly and its corresponding scaffolds. I have yet to implement gap size estimation or
resolution, so the sequences seen have grey spaces that are 10 bases wide, representing unknown regions between known
scaffolds. I used ImGui to create interactive elements such as the scrollable sections, and I used relatively simple
GLFW/Open GL to represent the contigs and scaffolds. Finally, this module also included a complete project structure
overhaul, as it moved files into `src` and `include` directories to help make file paths more intuitive. Prior to this
change, I was constantly confusing myself with where the file I was working on was located relative to the graphics
modules or project-wide files. 

Notable name changes:
- All files were changed from PascalCase to snake_case.
- DataLoader/Exporter were renamed to vis_exporter/loader to help clear up what type of data it was used for.
- All the assembly pipelines were prefixed with the strategy and filetype they used.

As far as current bugs go, I have an issue where when my sequence wraps on the second tab
of the GUI, only the upper sequence is responding to cursor placement. At the time of writing this, I have no idea as
to why this is happening, and I need to find my programming duck to fix this issue. I plan to solve this before Sunday's
submission, but I am writing this here in case it is not fixed.

As far as future work goes, I would like to implement an easier file selection system for the user, as the current
filepath system is not intuitive and can be very frustrating to use. Additionally, I plan on implementing gap size
estimation and resolution, as this will allow for me to construct a "polished" finished genome. For this process, I 
currently would like to visualize the steps of the resolution process step-by-step using another header in my visdata 
files. Finally, I need to build in more support for scaffold tracing from sink nodes accurately, as that is the main
source of error in my reconstructed products. 

For this module, I think that I have earned approximately 100 points, as I have two functional main programs for each of
my main processes. The graphical user interface is interactive with the inclusion of a scrollable timeline and
information panels when scaffolds/contigs are clicked on. I understand if I lose points on user-friendliness for the
file input section, which I understand, but I think that the rest of the program is relatively easy to navigate.
Finally, I think that my code is organized into folders that make sense given the scope of a file's functionality. I
made sure to have graphical, data, and core algorithms separated, as well as separating cpp and header files.

### References

| Module / Area        | Topic                   | Link                                                                 |
|----------------------|-------------------------|----------------------------------------------------------------------|
| general              | Usage of static         | https://www.geeksforgeeks.org/cpp/static-keyword-cpp/               |
| general              | Documentation guidelines | https://developer.lsst.io/cpp/api-docs.html                         |
| general              | Usage of auto           | https://www.geeksforgeeks.org/cpp/type-inference-in-c-auto-and-decltype/ |
| general              | Exceptions              | https://www.geeksforgeeks.org/cpp/how-to-throw-custom-exception-in-cpp/ |
| dna_sequence         | Helper functions        | https://www.w3tutorials.net/blog/what-are-helper-functions-in-c/     |
| dna_sequence         | size_t                  | https://www.geeksforgeeks.org/cpp/difference-between-int-and-size_t-in-cpp/ |
| dna_sequence         | Switch and cases        | https://www.w3schools.com/cpp/cpp_switch.asp                        |
| kmer_encoding        | Bitshift and masking    | https://www.geeksforgeeks.org/cpp/left-shift-right-shift-operators-c-cpp/ |
| open_addressing_table | Iterator implementation | https://stackoverflow.com/questions/46431762/how-to-implement-standard-iterators-in-class |
| de_bruijn_graph      | Usage of auto           | Reference source from Module 1                                       |
| de_bruijn_graph      | Use of [[nodiscard]]    | https://stackoverflow.com/questions/76489630/explanation-of-nodiscard-in-c17 |
| de_bruijn_graph      | Use of explicit         | https://www.geeksforgeeks.org/cpp/use-of-explicit-keyword-in-cpp/    |
| general              | Static cast             | https://www.geeksforgeeks.org/cpp/static_cast-in-cpp/   |
| graphics             | ImGui                   | http://imgui.net/articles/RoadMap.html  |
| gui                  | ImGui commands          | http://imgui.net/api/index.html  |
| general              | uint128 conversion      | https://www.geeksforgeeks.org/solidity/explicit-conversions-in-solidity/   |
| general              | Max value handling      | https://www.geeksforgeeks.org/cpp/integer-literal-in-c-cpp-prefixes-suffixes/   |

