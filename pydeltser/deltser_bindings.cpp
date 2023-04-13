#include <stdio.h>
#include <iostream>

#include "../deltser.cpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(pydelt, m) {

  m.doc() = "Python interface for deltser";

  m.def("run_deltser", [](std::vector<std::vector<std::vector<value_t>>>& faces,
                                bool approx, uint32_t approx_val,
                                bool in_file, char* filename) {

    // Save std::cout status and disable
    auto cout_buff = std::cout.rdbuf();
    std::cout.rdbuf(nullptr);

    std::vector<size_t> cell_counts;
    py::dict output;
    delta_complex_t complex = delta_complex_t();

    if(in_file){
        complex = delta_complex_t(filename);
    } else{
        complex = delta_complex_t(faces);
    }
    complex.compute_oldest_cofaces();

    //initialise approximate functionality
    size_t max_entries = std::numeric_limits<size_t>::max();
    if (approx) max_entries = approx_val;

    deltser thedeltser = deltser(&complex,filename,max_entries,true);
    thedeltser.compute_barcodes();

    for(int i = 0; i < complex.cells.size(); i++){
        cell_counts.push_back(complex.cells[i].size());
    }

    output["cell_counts"] = cell_counts;
    output["finite_pairs"] = thedeltser.finite_pairs;
    output["infinite_pairs"] = thedeltser.infinite_pairs;
    output["bettis"] = thedeltser.num_infinite_pairs();

    // Re-enable again cout
    std::cout.rdbuf(cout_buff);

    return output;
  });
}
