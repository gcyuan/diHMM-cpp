//
// Created by Stephanos Tsoucas on 6/15/17.
//

#ifndef BOOSTPYTHONHELLOWORLD_MODEL_PY_H
#define BOOSTPYTHONHELLOWORLD_MODEL_PY_H

#include "Annotation_Py.h"
#include "Arma_Numpy_Conversions.h"
#include "Forward_Backward.h"
#include "Model.h"
#include "Viterbi.h"

namespace np = boost::python::numpy;

void test_gpu() {
  mat m1 = randu(100, 100);
  mat m2 = randu(100, 100);
  mat m3 = kron(m1, m2);

  std::cout << m3 << std::endl;
}

class Model_Py {

  const Model m;

public:
  Model_Py(const Model m) : m(m) {}

  np::ndarray domain_transition_probabilities() {
    return from_mat(m.domain_transition_probabilities);
  }

  np::ndarray bin_transition_probabilities() {
    return from_cube(m.bin_transition_probabilities);
  }

  np::ndarray emission_probabilities() {
    return from_mat(m.emission_probs.emissions_probabilities_m());
  }

  const uword domain_size() const { return m.domain_size; }

  const uword n_bin_states() const { return m.n_bin_states; }

  const uword n_domain_states() const { return m.n_domain_states; }

  const Model &model() const { return m; }
};

mat bin_data(std::string filepath, const uword domain_size) {
  std::ifstream infile(filepath);
  std::cout << "FILE PATH " << filepath << std::endl;
  // TODO: see if file exists and exit gracefully
  std::string genome_chromosome;
  std::getline(infile, genome_chromosome);

  std::string chromosome_names;
  std::vector<std::string> histone_marks;
  std::istringstream iss(chromosome_names);
  std::getline(infile, chromosome_names);

  for (std::string name; iss >> chromosome_names;) {
    histone_marks.push_back(name);
  }

  std::string bin_data_string;
  std::vector<std::vector<unsigned int>> bin_d;

  while (std::getline(infile, bin_data_string)) {
    std::istringstream iss(bin_data_string);
    std::vector<unsigned int> bin_data_row;
    for (int histone; iss >> histone;) {
      bin_data_row.push_back(histone);
    }
    bin_d.push_back(bin_data_row);
  }

  const uword n_histone_marks = bin_d.at(0).size();

  const unsigned long n_bins = bin_d.size();

  const unsigned long n_bins_cropped = n_bins - (n_bins % domain_size);

  mat bin_data_m(n_histone_marks, n_bins_cropped);

  // Initialize (histone mark, bin position) matrix
  std::vector<std::vector<unsigned int>>::const_iterator it;
  int column_index = 0;
  for (it = bin_d.begin(); it != bin_d.end() && column_index < n_bins_cropped;
       it++, column_index++) {
    std::vector<double> col;
    std::vector<unsigned int>::const_iterator c;
    for (c = (*it).begin(); c < (*it).end(); c++) {
      col.push_back(double(*c));
    }
    bin_data_m.col(column_index) = vec(col);
  }

  return bin_data_m;
}

bool mat_comparator(const mat &a, const mat &b) { return a.n_cols > b.n_cols; }

Model_Py run_dihmm_existing_py(Model_Py &m, const int max_iter,
                               const double tolerance,
                               boost::python::list filepaths) {

  std::vector<mat> input_chromosomes;
  const int domain_size = m.model().domain_size;
  for (int i = 0; i < boost::python::len(filepaths); ++i) {
    const char *p = boost::python::extract<const char *>(filepaths[i]);
    const std::string path(p);
    input_chromosomes.push_back(bin_data(path, domain_size));
  }

  std::sort(input_chromosomes.begin(), input_chromosomes.end(), mat_comparator);
  for (auto &a : input_chromosomes) {
    std::cout << a.n_cols << std::endl;
  }

  return Model_Py(run_dihmm<Forward_Backward>(m.model(), input_chromosomes,
                                              max_iter, tolerance));
}

Model_Py run_dihmm_py(const int n_bin_states, const int n_domain_states,
                      const int domain_size, const int max_iter,
                      const double tolerance, boost::python::list filepaths) {

  std::vector<mat> input_chromosomes;

  for (int i = 0; i < boost::python::len(filepaths); ++i) {
    const char *p = boost::python::extract<const char *>(filepaths[i]);
    const std::string path(p);
    input_chromosomes.push_back(bin_data(path, domain_size));
  }

  std::sort(input_chromosomes.begin(), input_chromosomes.end(), mat_comparator);
  for (auto &a : input_chromosomes) {
    std::cout << a.n_cols << std::endl;
  }

  return Model_Py(run_dihmm<Forward_Backward>(n_bin_states, n_domain_states,
                                              domain_size, max_iter, tolerance,
                                              input_chromosomes));
}

Model_Py load_model_py(const std::string &directory, const int domain_size) {
  return load_model(directory, domain_size);
}

bool save_model_py(const Model_Py &m, const std::string &directory) {
  return save_model(m.model(), directory);
}

boost::python::list annotate_py(const Model_Py &m,
                                boost::python::list filepaths) {

  std::vector<mat> bin_d;

  for (int i = 0; i < boost::python::len(filepaths); i++) {
    const char *p = boost::python::extract<const char *>(filepaths[i]);
    std::string path(p);
    bin_d.push_back(bin_data(path, m.domain_size()));
  }

  std::vector<Annotation> as = annotate(m.model(), bin_d);

  boost::python::list ret;

  for (const auto &a : as) {
    ret.append(Annotation_Py(a));
  }

  return ret;
}

#endif // BOOSTPYTHONHELLOWORLD_MODEL_PY_H
