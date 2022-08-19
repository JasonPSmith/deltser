# deltser
A modified version of flagser that computes persistent homology of delta complexes, and associated python wrapper.

To install just deltser (without python wrapper), then run:
```sh
git clone https://github.com/JasonPSmith/tournser.git
make
```

If the python wrapper is required, first ensure python packages numpy and pybind11 are installed, then run:
```sh
git clone --recursive https://github.com/JasonPSmith/tournser.git
make
pip install .
``` 

To run deltser simply do
```sh
./deltser in_address out_address approx_val
```
where in_address is the address of a file containing the delta complex (see format below), out_address is where to save the output, and approx_val is an optional entry that approximates the homology. See [flagser-count]([)https://github.com/JasonPSmith/flagser-count) for more details on approximate.

For example:
```sh
./deltser ./Examples/simple_example.dlt ./Examples/simple_example.dlt
```

The format of the input file is as follows, first line is "dim 0", second line is a list of the filtration values of the vertices, then for each subsequent dimension there is a line "dim i" followed by all faces of dimension i, where each line is an i-face given by the boundary of that face by it's location in the dim i-1 list. 

For example, consider the complex constructed from the complete bidirectional digraph on 3 vertices, where the 3-cycles are the 2-simplices (so two 2-simplices with the same vertices, but different edges). The input to deltser would be:

dim 0
0 0 0
dim 1
0 1
1 2
2 0
1 0
2 1
0 2
dim 2
0 1 2
3 4 5

where the zeroes under dim 0 indicate the number of vertices, under dim 1 is a list of directed edges, and under dim 2 are the 2 simplices where "0 1 2" means the simplex whose boundary is the 1st, 2nd and 3rd edge, which is (0,1),(1,2),(2,0). More examples can be seen in the Examples folder.


To use pydeltser, in python run, for example:
```sh
from pydeltser import *
faces=[[[0],[0],[0]],[[0, 1, 0],[1, 2, 0],[2, 0, 0],[1, 0, 0],[2, 1, 0],[0, 2, 0]],[[0, 1, 2, 0],[3, 4, 5, 0]]]
deltser(faces)
```
This will compute the same complex as in Examples/simple_example.dlt. So the entry is a list of lists, where faces[0] is the filtration values of the vertices, each within a list of length 1, and faces[i] contains all the entries under dim i of the format described above for C++ deltser.

There is also deltser_file, which takes as input the address of dlt file. For example
```sh
from pydeltser import *
deltser('./Examples/simple_example.dlt')
```

The output of both pydeltser functions is a dictionary with entries: `cell\_counts', `finite\_pairs', `infinte\_pairs' and `bettis'. Where the i'th entry of each corresponds to the i'th dimension. And an entry of infinite pairs is a single number denoting the birth time of the pair.
