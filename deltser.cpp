//-------------------------------------------------------------------------//
//BEGIN header

//#define INDICATE_PROGRESS
#define SORT_COLUMNS_BY_PIVOT

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <thread>
#include <unordered_set>
#include <functional>

typedef float value_t;
typedef int64_t index_t;
typedef int16_t coefficient_t;
typedef index_t entry_t;
typedef std::pair<value_t, index_t> filtration_index_t;

typedef std::deque<index_t> pivot_column_index_t;
const index_t INVALID_INDEX = std::numeric_limits<index_t>::max();

float string_to_float(std::string s) { return atof(s.c_str()); }

const index_t get_index(const entry_t& i) { return i; }
index_t get_coefficient(const entry_t& i) { return 1; }
entry_t make_entry(index_t _index, coefficient_t _value) { return entry_t(_index); }

const entry_t& get_entry(const entry_t& e) { return e; }

value_t get_filtration(const filtration_index_t& i) { return i.first; }
index_t get_index(const filtration_index_t& i) { return i.second; }

class filtration_entry_t : public std::pair<value_t, entry_t> {
public:
	filtration_entry_t() {}
	filtration_entry_t(const entry_t& e) : std::pair<value_t, entry_t>(0, e) {}
    filtration_entry_t(index_t _filtration, index_t _index)
        : std::pair<value_t, entry_t>(_filtration, entry_t(_index)) {}
	filtration_entry_t(value_t _filtration, index_t _index, coefficient_t _coefficient)
	    : std::pair<value_t, entry_t>(_filtration, make_entry(_index, _coefficient)) {}
	filtration_entry_t(const filtration_index_t& _filtration_index, coefficient_t _coefficient)
	    : std::pair<value_t, entry_t>(get_filtration(_filtration_index),
	                                 make_entry(get_index(_filtration_index), _coefficient)) {}
	filtration_entry_t(const filtration_index_t& _filtration_index)
	    : filtration_entry_t(_filtration_index, 1) {}
};

template <typename Heap> filtration_entry_t pop_pivot(Heap& column) {
    if (column.empty())
        return filtration_entry_t(-1);
    else {
        auto pivot = column.top();
        column.pop();
        while (!column.empty() &&
               get_index(column.top()) == get_index(pivot)) {
            column.pop();

            if (column.empty())
                return filtration_entry_t(-1);
            else {
                pivot = column.top();
                column.pop();
            }
        }
        return pivot;
    }
}

const entry_t& get_entry(const filtration_entry_t& p) { return p.second; }
const index_t get_index(const filtration_entry_t& p) { return get_index(get_entry(p)); }
const coefficient_t get_coefficient(const filtration_entry_t& p) { return get_coefficient(get_entry(p)); }
const value_t& get_filtration(const filtration_entry_t& p) { return p.first; }

struct greater_filtration_or_smaller_index {
	bool operator()(const filtration_index_t a, const filtration_index_t b) {
		return (get_filtration(a) > get_filtration(b)) ||
		       ((get_filtration(a) == get_filtration(b)) && (get_index(a) < get_index(b)));
	}
};

template <typename Entry> struct smaller_index {
	bool operator()(const Entry& a, const Entry& b) { return get_index(a) < get_index(b); }
};

class filtered_union_find {
	std::vector<index_t> parent;
	std::vector<std::vector<index_t>> rank;
	const std::vector<value_t> filtration;

public:
	filtered_union_find(const std::vector<value_t>& _filtration)
	    : rank(_filtration.size()), filtration(_filtration), parent(_filtration.size()) {
		for (index_t i = 0; i < _filtration.size(); ++i){
			parent[i] = i;
			rank[i] = std::vector<index_t>{i};
		}
	}
	index_t find(index_t x) {
		return parent[x];
	}
	value_t link(index_t x, index_t y) {
		x = find(x);
		y = find(y);
		if (x == y) return -1;
		if (filtration[x] < filtration[y] || (filtration[x] == filtration[y] && rank[x].size() > rank[y].size())){
			for(auto i : rank[y]) parent[i] = x;
			rank[x].insert(rank[x].end(),rank[y].begin(),rank[y].end());
			return filtration[y];
		} else {
			for(auto i : rank[x]) parent[i] = y;
			rank[y].insert(rank[y].end(),rank[x].begin(),rank[x].end());
			return filtration[x];
		}
	}
};

template <typename ValueType> class compressed_sparse_matrix {
	std::deque<size_t> bounds;
	std::deque<ValueType> entries;

public:
	size_t size() const { return bounds.size(); }

	void clear() {
		bounds.clear();
		bounds.shrink_to_fit();
		entries.clear();
		entries.shrink_to_fit();
	}

	typename std::deque<ValueType>::const_iterator cbegin(size_t index) const {
		assert(index < size());
		return index == 0 ? entries.cbegin() : entries.cbegin() + bounds[index - 1];
	}

	typename std::deque<ValueType>::const_iterator cend(size_t index) const {
		assert(index < size());
		return entries.cbegin() + bounds[index];
	}

	template <typename Iterator> void append_column(Iterator begin, Iterator end) {
		for (Iterator it = begin; it != end; ++it) { entries.push_back(*it); }
		bounds.push_back(entries.size());
	}

	void append_column() { bounds.push_back(entries.size()); }

	void push_back(ValueType e) {
		assert(0 < size());
		entries.push_back(e);
		++bounds.back();
	}

	void pop_back() {
		assert(0 < size());
		entries.pop_back();
		--bounds.back();
	}

	template <typename Collection> void append_column(const Collection collection) {
		append_column(collection.cbegin(), collection.cend());
	}
};

//END header
//-------------------------------------------------------------------------//
//BEGIN delta_complex




class delta_complex_t;

template <typename t>
std::vector<t> split(const std::string& s, char delim, const std::function<t(std::string)>& transform) {
	std::vector<t> elems;
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) elems.push_back(transform(item));
	return elems;
}

class delta_complex_cell_t {
public:
    std::vector<delta_complex_cell_t*> boundary;
    std::vector<index_t> coboundary;
	std::unordered_set<delta_complex_cell_t*> children;
    index_t location;
    value_t filtration;
    std::vector<index_t> vertices;
    int dimension;
	index_t oldest_coface;

    //initialise class
    delta_complex_cell_t(int _dim, index_t _v, value_t _filt) : dimension(_dim), filtration(_filt), location(_v), oldest_coface(-1) {vertices.push_back(_v);}
    delta_complex_cell_t(int _dim, std::vector<delta_complex_cell_t*> _cells,
        value_t _filt, value_t _loc) : dimension(_dim), boundary(_cells), filtration(_filt), location(_loc), oldest_coface(-1) {}

    void set_children(){
        for ( auto b : boundary){
            children.insert(b);
            children.insert(b->children.begin(),b->children.end());
            b->coboundary.push_back(location);
        }
    }

    void set_vertices(){
        for ( auto c : children ){
            if ( c->dimension == 0 ){
                vertices.push_back(c->vertices.front());
            }
        }
    }
    size_t coboundary_size(){
        return coboundary.size();
    }
    void set_filtration(value_t filter_val){
        filtration = filter_val;
    }
	void compute_oldest_coface(delta_complex_t* complex);
}; //END delta_complex_cell_t

class delta_complex_t {
public:
    //stores all cells as a vector of vectors, each vector being all cells of that dimension
    std::vector<std::vector<delta_complex_cell_t>> cells;

    //initialised with a string s with the address of a file containing all simplices
    //the format of this list is:
    //dim 0, followed by a line with an int for each vertex, the value of which is the filtration
    //after any subsequent dim i each line contained a list of facets of the simplex, followed
    //by a filtration value, the facets are indexed by their position in the list for dim i-1
    delta_complex_t(){};
	delta_complex_t(std::string s){
        std::ifstream infile;
        infile.open(s);

        index_t current_dimension = 0;
        int val = 0;
        std::string line;

        //read through all lines of input file, if a dim i is encountered then set
        //current dimension to i
        while (std::getline(infile,line)){
            if (line.length() == 0) continue;
            if (line[0] == 'd' && line[1] == 'i' && line[2] == 'm') {
                current_dimension = (int)line[4] - '0';
                cells.push_back( std::vector<delta_complex_cell_t>() );
            }
            //for dim 0 create a vertex cell for each entry on that line with filtration
            //value given by the entry
            else if (current_dimension == 0) {
                std::vector<value_t> vertex_filtration = split<value_t>(line, ' ', string_to_float);
                for( auto v : vertex_filtration ){
                    cells.back().push_back(delta_complex_cell_t(0,cells.back().size(),v));
                }
            //for all larger dimensions create a cell by giving the boundary faces and filtration value
        	} else {
        		std::vector<int> faces = split<int>(line, ' ', string_to_float);
                val = faces.back();
                faces.pop_back();
                std::vector<delta_complex_cell_t*> new_face;
                for (auto p : faces){
                    new_face.push_back(&cells[current_dimension-1][p]);
                }
                value_t loc = cells.back().size();
                cells.back().push_back(delta_complex_cell_t(current_dimension,new_face,val,loc));
                cells.back().back().set_children();
                cells.back().back().set_vertices();
        	}
        }
	}
	delta_complex_t(std::vector<std::vector<std::vector<value_t>>>& faces){
        int val = 0;
        cells.resize(faces.size());

        //for dim 0 create a vertex cell for each entry on that line with filtration
	    //value given by the entry
        for( auto v : faces[0] ){
            cells[0].push_back(delta_complex_cell_t(0,cells[0].size(),v[0]));
        }
        //read through all lines of input file, if a dim i is encountered then set
        //current dimension to i
        for (index_t current_dimension = 1; current_dimension < faces.size(); current_dimension++){
            for (auto f : faces[current_dimension]){
                val = f.back();
                f.pop_back();
                std::vector<delta_complex_cell_t*> new_face;
                for (auto p : f){
                    new_face.push_back(&cells[current_dimension-1][p]);
                }
                value_t loc = cells[current_dimension].size();
                cells[current_dimension].push_back(delta_complex_cell_t(current_dimension,new_face,val,loc));
                cells[current_dimension].back().set_children();
                cells[current_dimension].back().set_vertices();
        	}
        }
	}
    index_t number_of_cells(index_t dimension) const {
        if ( dimension >= cells.size() ) { return 0; }
        return cells[dimension].size();
    }
    bool is_top_dimension(index_t dimension){
        return cells.size() == dimension-1;
    }
    index_t top_dimension(){
        return cells.size()-1;
    }
    value_t filtration(index_t dimension,index_t index){
        return cells[dimension][index].filtration;
    }
	std::vector<value_t> vertex_filtration(){
		std::vector<value_t> out;
		for( auto c : cells[0] ){
			out.push_back(c.filtration);
		}
        return out;
    }
    delta_complex_cell_t* get(index_t dimension,index_t index){
        return &cells[dimension][index];
    }
	void compute_oldest_cofaces(){
		for (auto p : cells){
			for(auto q : p){
				q.compute_oldest_coface(this);
			}
		}
	}
//END delta_complex_t
};

void delta_complex_cell_t::compute_oldest_coface(delta_complex_t* complex){
	 value_t oldest = -1;
	for( auto c : coboundary ){
		value_t f = complex->get(dimension+1,c)->filtration;
		if( f > oldest ){
			oldest = f;
			oldest_coface = c;
		}
	}
}

class simplex_coboundary_enumerator {
private:
	index_t  idx_above, dim;
	delta_complex_cell_t* simplex;
	delta_complex_t* complex;

public:
	simplex_coboundary_enumerator(const filtration_entry_t _simplex, index_t _dim,
	                              delta_complex_t* _complex)
	    : idx_above(0), simplex(_complex->get(_dim,_simplex.second)), dim(_dim), complex(_complex) {}

	bool has_next() {
		return idx_above < simplex->coboundary_size();
	}

	filtration_entry_t next() {
        idx_above++;
        return filtration_entry_t(complex->get(dim+1,simplex->coboundary[idx_above-1])->filtration, simplex->coboundary[idx_above-1]);
	}
};

#ifdef SORT_COLUMNS_BY_PIVOT
struct greater_filtration_or_better_pivot_or_smaller_index {
	greater_filtration_or_better_pivot_or_smaller_index(delta_complex_t* _complex, index_t _dimension) : complex(_complex), dimension(_dimension) {}
	bool operator()(filtration_index_t a, filtration_index_t b) const {
		// First order by the filtration value
		if (get_filtration(a) > get_filtration(b)) return true;
		if (get_filtration(a) < get_filtration(b)) return false;

		auto ta = get_coboundary_size_and_gap_and_pivot(get_index(a));
		auto tb = get_coboundary_size_and_gap_and_pivot(get_index(b));

		// Then the number of non-trivial coboundary entries
		if (std::get<0>(ta) < std::get<0>(tb)) return true;
		if (std::get<0>(ta) > std::get<0>(tb)) return false;

		// Then order by the better pivoting
		if (std::get<2>(ta) < std::get<2>(tb)) return true;
		if (std::get<2>(ta) > std::get<2>(tb)) return false;

		if (std::get<1>(ta) > std::get<1>(tb)) return true;
		if (std::get<1>(ta) < std::get<1>(tb)) return false;

		// Finally, order by their indices
		return get_index(a) < get_index(b);
	}

private:
	delta_complex_t* complex;
	index_t dimension;

	// A column is considered to be a better pivot if the jump from pivot to the next
	// non-trivial element is as big as possible. This prevents accidentally inserting
	// non-trivial elements just below the pivot, which sometimes creates very long
	// reduction chains.
	// The second sort criterium is for it to be small because the small pivots will be
	// used the most.
	std::tuple<size_t, size_t, index_t> get_coboundary_size_and_gap_and_pivot(index_t a) const {
		// Look at the first two gaps of the pivot and the next element
		index_t pivot = 0;
		size_t gap_after_pivot = 0;
		simplex_coboundary_enumerator iterator(a,dimension,complex);
		size_t coboundary_size = 0;
		while (iterator.has_next()) {
			coboundary_size++;
			index_t next_index = get_index(iterator.next().second);
			if (next_index > pivot) {
				gap_after_pivot = next_index - pivot;
				pivot = next_index;
			}
		}

		return std::make_tuple(coboundary_size, gap_after_pivot, pivot);
	}
};
#endif

// This class is just an ordinary priority queue, but once the
// queue gets too long (because a lot of faces are inserted multiple
// times) it starts collecting the coefficients and only inserting each
// new face once
template <class Container, class Comparator>
class priority_queue_t : public std::priority_queue<filtration_entry_t, Container, Comparator> {
	std::unordered_map<index_t, bool> coefficients;
	static const filtration_entry_t dummy;
	bool use_dense_version = false;
	size_t dense_threshold;

public:
	priority_queue_t(size_t _dense_threshold)
	    : dense_threshold(_dense_threshold) {}

	void push(const filtration_entry_t& value) {
		if (use_dense_version) {
			// If we already have this value: update the count and don't push it again
			auto p = coefficients.find(get_index(value));
			if (p != coefficients.end()) {
				p->second = !p->second;
				return;
			}
		}

		std::priority_queue<filtration_entry_t, Container, Comparator>::push(value);

		if (use_dense_version) coefficients.insert(std::make_pair(get_index(value), get_coefficient(value)));

		if (!use_dense_version &&
		    std::priority_queue<filtration_entry_t, Container, Comparator>::size() >= dense_threshold)
			use_dense_version = true;
	}

	void pop() {
		// Don't use this, only allow get_pivot
		throw std::exception();
	}

	filtration_entry_t pop_pivot() {
		remove_trivial_coefficient_entries();
		if (std::priority_queue<filtration_entry_t, Container, Comparator>::empty())
			return dummy;
		else {
			auto pivot = get_top();
			safe_pop();
			while (!std::priority_queue<filtration_entry_t, Container, Comparator>::empty() &&
			       get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()) ==
			           get_index(pivot)) {
				safe_pop();
				remove_trivial_coefficient_entries();

				if (std::priority_queue<filtration_entry_t, Container, Comparator>::empty())
					return dummy;
				else {
					pivot = get_top();
					safe_pop();
				}
			}
			return pivot;
		}
	}

	filtration_entry_t get_pivot() {
		filtration_entry_t result = pop_pivot();
		if (get_index(result) != -1) { push(result); }
		return result;
	}

private:
	inline filtration_entry_t get_top() {
		auto pivot = std::priority_queue<filtration_entry_t, Container, Comparator>::top();
		return pivot;
	}

	inline void safe_pop() {
		if (use_dense_version) {
			auto e =
			    coefficients.find(get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()));
			if (e != coefficients.end()) coefficients.erase(e);
		}
		std::priority_queue<filtration_entry_t, Container, Comparator>::pop();
	}

	inline void remove_trivial_coefficient_entries() {
		if (use_dense_version) {
			auto p = coefficients.find(get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()));
			while (p != coefficients.end() && p->second == false) {
				coefficients.erase(p);
				std::priority_queue<filtration_entry_t, Container, Comparator>::pop();
				p = coefficients.find(get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()));
			}
		}
	}
};
template <class Container, class Comparator>
const filtration_entry_t priority_queue_t<Container, Comparator>::dummy(filtration_entry_t(0.0, -1));



//END delta_complex
//-------------------------------------------------------------------------//
//BEGIN deltser


class deltser {
	delta_complex_t* complex;
	std::ofstream outfile;
	index_t n, dim_max;
	bool python;
	mutable std::vector<filtration_entry_t> coface_entries;
	size_t max_entries;
	std::vector<size_t> skipped_entries;

public:
	std::vector<std::vector<std::pair<value_t,value_t>>> finite_pairs;
	std::vector<std::vector<value_t>> infinite_pairs;
	deltser(	delta_complex_t* _complex, char* _outname, size_t _max_entries, bool _python)
	    : complex(_complex), n(complex->number_of_cells(0)),
	      dim_max(complex->top_dimension()),
		  max_entries(_max_entries),
		  python(_python) {
				if(!python) outfile.open(_outname);
                skipped_entries.assign(dim_max+1,0);
				infinite_pairs.resize(dim_max+1);
				finite_pairs.resize(dim_max+1);
			  }

	value_t compute_filtration(const index_t index, index_t dim) const {
        return complex->get(dim,index)->filtration;
	}

	void assemble_columns_to_reduce(std::vector<filtration_index_t>& simplices,
	                                std::vector<filtration_index_t>& columns_to_reduce,
	                                pivot_column_index_t& pivot_column_index, index_t dim, index_t num_simplices);

	void compute_dim_0_pairs(std::vector<filtration_index_t>& edges,
	                         std::vector<filtration_index_t>& columns_to_reduce) {

        //Get all edges and sort them
		filtered_union_find dset(complex->vertex_filtration());
		edges = get_edges();
		std::sort(edges.rbegin(), edges.rend(),
		          greater_filtration_or_smaller_index());

		if(!python) { outfile << "persistence intervals in dim 0:" << std::endl; }

		for (auto e : edges) {
			value_t birth = dset.link(complex->get(1,get_index(e))->vertices[0], complex->get(1,get_index(e))->vertices[1]);
			if (birth != -1) {
				if (get_filtration(e) > birth) {
					if(!python){
						outfile << " [" << birth << ", " << get_filtration(e) << ")" << std::endl;
					} else {
						finite_pairs[0].push_back(std::make_pair(birth,get_filtration(e)));
					}
				}
			} else {
				columns_to_reduce.push_back(e);
			}
		}
		std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

		for (index_t i = 0; i < n; ++i)
			if (dset.find(i) == i){
				if(!python){
				    outfile << " [0, )" << std::endl << std::flush;
				}
				infinite_pairs[0].push_back(0);
			}
	}

	template <typename Column, typename Iterator>
	filtration_entry_t add_coboundary_and_get_pivot(Iterator column_begin, Iterator column_end,
	                                              Column& working_coboundary, const index_t& dim,
                                                  std::priority_queue<filtration_entry_t, std::deque<filtration_entry_t>, smaller_index<filtration_entry_t>>&);

	void sort_columns(std::vector<filtration_index_t>& columns_to_reduce, index_t dimension) {
#ifdef SORT_COLUMNS_BY_PIVOT
		std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
			greater_filtration_or_better_pivot_or_smaller_index(complex,dimension));
#else
		std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
	        greater_filtration_or_smaller_index());
#endif
}

	void compute_pairs(std::vector<filtration_index_t>& columns_to_reduce,
	                   pivot_column_index_t& pivot_column_index, index_t dim) {

		if(!python) { outfile << "# persistence intervals in dim " << dim << ":" << std::endl; }

		std::cout << "Computing Dimension " << dim << std::endl;
        compressed_sparse_matrix<filtration_entry_t> reduction_coefficients;

		std::vector<filtration_entry_t> coface_entries;

		for (index_t index_column_to_reduce = 0; index_column_to_reduce < columns_to_reduce.size();
		     ++index_column_to_reduce) {
			auto column_to_reduce = columns_to_reduce[index_column_to_reduce];
            std::priority_queue<filtration_entry_t, std::deque<filtration_entry_t>, smaller_index<filtration_entry_t>> reduction_column;

			priority_queue_t<std::vector<filtration_entry_t>,
			                    greater_filtration_or_smaller_index>
			    working_coboundary(columns_to_reduce.size());

			value_t filtration = get_filtration(column_to_reduce);

#ifdef INDICATE_PROGRESS
			if ((index_column_to_reduce + 1) % 1000000 == 0)
				std::cout << "\033[K"
				          << "reducing column " << index_column_to_reduce + 1 << "/"
				          << columns_to_reduce.size() << " (filtration " << filtration << ")"
				          << std::flush << "\r";
#endif

			index_t index_column_to_add = index_column_to_reduce;

			filtration_entry_t pivot;

            reduction_coefficients.append_column();
            reduction_coefficients.push_back(filtration_entry_t(column_to_reduce, 1));

			while (true) {
                auto reduction_column_begin = reduction_coefficients.cbegin(index_column_to_add);
                auto reduction_column_end = reduction_coefficients.cend(index_column_to_add);

				pivot = add_coboundary_and_get_pivot(
				    reduction_column_begin, reduction_column_end,
				    working_coboundary, dim, reduction_column);

				if (get_index(pivot) > -1) {
                    auto pivot_column_idx = pivot_column_index[get_index(pivot)];

                    if (pivot_column_idx != INVALID_INDEX) {
                        index_column_to_add = pivot_column_idx;
                        continue;
                    } else {
						value_t death = get_filtration(pivot);
						if (death > filtration) {
							if(!python){
							    outfile << " [" << filtration << ", " << death << ")" << std::endl << std::flush;
						    }
						    else { finite_pairs[dim].push_back(std::make_pair(filtration,death)); }
						}

                        pivot_column_index[get_index(pivot)] =  index_column_to_reduce;
                            reduction_coefficients.pop_back();
            				while (true) {
            					filtration_entry_t e = pop_pivot(reduction_column);
            					if (get_index(e) == -1) break;
            					reduction_coefficients.push_back(e);
            				}
						break;
					}
				} else if(get_index(pivot) == -1) {
					if(!python){
						outfile << " [" << filtration << ", )" << std::endl << std::flush;
					}
					infinite_pairs[dim].push_back(filtration);
					break;
				}else {
					//outfile << "[?" << filtration << ", " << ", ?)" << std::endl;
					skipped_entries[dim]++;
					break;
				}
			}
		}
#ifdef INDICATE_PROGRESS
		std::cout << "\033[K";
#endif
	}

	std::vector<filtration_index_t> get_edges();

	std::vector<value_t> num_infinite_pairs(){
		std::vector<value_t> out;
		for(auto i : infinite_pairs) out.push_back(i.size());
		return out;
	}

	void print_summary(){
		if(!python){
			std::vector<value_t> inf_pairs = num_infinite_pairs();
			outfile << std::endl;
			outfile << "# Betti Numbers:" << std::endl;
			for(index_t i = 0; i <= dim_max; i++){
				outfile << "#        dim H_" << i << " = " << inf_pairs[i];
				if( skipped_entries[i] > 0){
					outfile << " : (" << skipped_entries[i] << " entries skipped)";
				}
				outfile << std::endl;
			}
			outfile << std::endl;
			outfile << "# Cell Counts:" << std::endl;
			for(index_t i = 0; i <= dim_max; i++){
				outfile << "#        dim C_" << i << " = " << complex->number_of_cells(i) << std::endl;
			}
		}
	}

	void compute_barcodes() {

		std::vector<filtration_index_t> simplices, columns_to_reduce;

		compute_dim_0_pairs(simplices, columns_to_reduce);

		for (index_t dim = 1; dim <= dim_max; ++dim) {
			pivot_column_index_t pivot_column_index(complex->number_of_cells(dim + 1), INVALID_INDEX);

			sort_columns(columns_to_reduce,dim);

			compute_pairs(columns_to_reduce, pivot_column_index, dim);
			if (dim < dim_max) {
				assemble_columns_to_reduce(simplices, columns_to_reduce, pivot_column_index,
				                           dim + 1, complex->number_of_cells(dim+1));
			}
		}
		print_summary();
	}
};

template <typename Column, typename Iterator>
filtration_entry_t deltser::add_coboundary_and_get_pivot(
    Iterator column_begin, Iterator column_end,
    Column& working_coboundary, const index_t& dim,
    std::priority_queue<filtration_entry_t, std::deque<filtration_entry_t>, smaller_index<filtration_entry_t>>& reduction_column) {
	index_t iterations = 0;
	for (auto it = column_begin; it != column_end; ++it) {
		filtration_entry_t simplex = *it;

        reduction_column.push(simplex);

		coface_entries.clear();
		simplex_coboundary_enumerator cofaces(simplex, dim, complex);
		while (cofaces.has_next()) {
			filtration_entry_t coface = cofaces.next();

            iterations++;
            working_coboundary.push(coface);
		}
		if (iterations > max_entries) {
			return filtration_entry_t(0,-2);
		}
	}

	return working_coboundary.get_pivot();
}

//returns a vector of all the edges where each is representated as a pair:
//(filtration,index), where the edge is the simplex at complex.get(1,index)
std::vector<filtration_index_t> deltser::get_edges() {
	std::vector<filtration_index_t> edges;
    int n = complex->number_of_cells(1);
	for ( index_t index = 0; index < n; index++) {
		edges.push_back(std::make_pair(complex->get(1,index)->filtration, index));
	}
	return edges;
}

void deltser::assemble_columns_to_reduce(
    std::vector<filtration_index_t>& simplices, std::vector<filtration_index_t>& columns_to_reduce,
    pivot_column_index_t& pivot_column_index, index_t dim, index_t num_simplices ) {

	columns_to_reduce.clear();

	for (index_t index = 0; index < num_simplices; ++index) {
		if (pivot_column_index[index] == INVALID_INDEX) {
			value_t filtration = compute_filtration(index, dim);
			columns_to_reduce.push_back(std::make_pair(filtration, index));
		}
	}

	//std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
	//          greater_filtration_or_smaller_index());
}


//END deltser
//-------------------------------------------------------------------------//
//BEGIN main


int main(int argc, char** argv) {
	//input takes the form: ./deltser in_address out_address approx_val

    //read in complex
	const char* filename = argv[1];
	char* outname = argv[2];
    delta_complex_t complex(filename);
	complex.compute_oldest_cofaces();

	//initialise approximate functionality
	size_t max_entries = std::numeric_limits<size_t>::max();
	if(argc > 3) max_entries = atoi(argv[3]);

    //create deltser object and compute persistent homology
	deltser(&complex,outname,max_entries,false).compute_barcodes();
}

//END main
//-------------------------------------------------------------------------//
