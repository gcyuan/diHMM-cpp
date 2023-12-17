//
// Created by Stephanos Tsoucas on 6/14/17.
//

#include "Model.h"

using namespace arma;

double squared_euclid(const vec &p1, const vec &p2) {
  return accu(square(p1 - p2));
}

/*
   This implementation clusters the domains by number of each bin state that
   occurs in the domain.
   TODO: cluster by the actual vector of bin states for each domain
 */
mat initial_domain_transition_probabilities(
    const std::vector<unsigned int> &domains_to_states, int n_domain_states) {

  mat domain_transition_probabilities = zeros(n_domain_states, n_domain_states);

  // For each domain transition, we count the number of times it appears
  const unsigned long n_domains = domains_to_states.size();

  for (int i = 1; i != n_domains; i++) {
    const unsigned int from = domains_to_states.at(i - 1);
    const unsigned int to = domains_to_states.at(i);

    double counts_so_far = domain_transition_probabilities(from, to);

    domain_transition_probabilities(from, to) = counts_so_far + 1.0;
  }

  // Get rid of zero counts- replace with a value smaller than the min for that
  // col In reference implementation, they take the global min- this takes the
  // min element of the column
  domain_transition_probabilities.replace(0.0, 0.5);

  mat total_counts = sum(domain_transition_probabilities, 1);
  domain_transition_probabilities.each_col() /= total_counts;

  return domain_transition_probabilities;
}

cube initial_bin_transition_probabilities(
    const std::vector<unsigned int> &bins_to_states,
    const std::vector<unsigned int> &domains_to_states,
    unsigned int n_bin_states, unsigned int n_domain_states,
    unsigned int domain_size) {

  // bin_transition_probabilities(dest_bin_state, curr_bin_state,
  // dest_domain_state) First, we fill up bin_transition_probabilities with the
  // counts of transitions Then, we divide by the total number of outcomes.
  cube bin_transition_probabilities =
      zeros(n_bin_states, n_bin_states, n_domain_states);

  const unsigned long n_bins = bins_to_states.size();

  for (int i = 1; i < n_bins; i++) {
    const int dest_domain_index = i / domain_size;
    const int dest_domain_state = domains_to_states[dest_domain_index];

    const int curr_bin_state = bins_to_states.at(i - 1);
    const int dest_bin_state = bins_to_states.at(i);

    bin_transition_probabilities(curr_bin_state, dest_bin_state,
                                 dest_domain_state) += 1;
  }

  bin_transition_probabilities.replace(0.0, 0.5);

  bin_transition_probabilities.each_slice([](mat &slice) {
    mat total_transition_counts = sum(slice, 1);
    slice.each_col() /= total_transition_counts;
  });

  return bin_transition_probabilities;
}

// TODO make sure the right matrices are initialized to zero
// TODO pass in matrix instead of vector vector
Model initial_model(const mat &bin_data_m, unsigned int n_bin_states,
                    unsigned int n_domain_states, unsigned int n_histone_marks,
                    const unsigned int domain_size) {

  // Using k-means clustering to initialize the bin-level states
  // Find the indexes of the columns where there is at least 1 non-zero element
  std::vector<int> non_zero_indexes;

  const uword n_bins_cropped = bin_data_m.n_cols;

  for (int i = 0; i < n_bins_cropped; i++) {
    double col_sum = accu(bin_data_m.col(i));
    if (col_sum > 0) {
      non_zero_indexes.push_back(i);
    }
  }

  // Get the non-zero columns in a new matrix- this matrix is used for k-means.
  mat non_zero_m = mat(n_histone_marks, non_zero_indexes.size());

  std::vector<int>::iterator i;
  int counter = 0;
  for (i = non_zero_indexes.begin(); i < non_zero_indexes.end(); i++) {
    non_zero_m.col(counter) = bin_data_m.col((*i));
    counter++;
  }

  // We only cluster with n - 1 states, since all columns with zero are assigned
  // state 0.
  const int non_zero_states = n_bin_states - 1;

  // mat centroids = non_zero_m.cols(0, non_zero_states);
  mat centroids;

  // means, data, k, seed_mode, n_iter, print_mode
  kmeans(centroids, non_zero_m, non_zero_states, static_subset, 11, false);

  // For each bin, find the index i of the closest centroid
  // (i + 1) becomes the state of the bin, since we already reserve state 0.
  std::vector<unsigned int> bins_to_bin_states_nz;

  const unsigned long n_nz_indexes = non_zero_indexes.size();

  non_zero_m.each_col(
      [centroids, non_zero_states, &bins_to_bin_states_nz](vec &c) mutable {
        double min_dist = squared_euclid(c, centroids.col(0));
        unsigned int min_index = 0;

        for (unsigned int i = 1; i < non_zero_states; i++) {
          double dist = squared_euclid(c, centroids.col(i));

          if (dist < min_dist) {
            min_dist = dist;
            min_index = i;
          }
        }
        bins_to_bin_states_nz.push_back(min_index);
      });

  std::vector<unsigned int>::iterator nz_indexes;

  // We have a vector of indexes that are not zero, and a vector storing the
  // corresponding values of the indexes.
  std::vector<unsigned int> bins_to_bin_states(n_bins_cropped, 0);
  std::vector<unsigned int> bin_state_counts(n_bin_states, 0);

  std::vector<unsigned int>::iterator values = bins_to_bin_states_nz.begin();
  std::vector<int>::iterator indexes;
  for (indexes = non_zero_indexes.begin(); indexes != non_zero_indexes.end();
       indexes++) {
    bins_to_bin_states[*indexes] = (*values) + 1;
    bin_state_counts[*values] += 1;
    values++;
  }

  const unsigned long n_domains = n_bins_cropped / domain_size;

  // Replace with the actual
  mat domains_to_bin_state_counts = zeros(domain_size, n_domains);

  // i is the index of the domain
  for (int i = 0; i < n_domains; i++) {

    const int domain_start_index = i * domain_size;

    vec domain = zeros(domain_size);

    for (int j = 0; j < domain_size; j++) {
        const unsigned long bin_state =
          bins_to_bin_states[domain_start_index + j];

        domain(j) = bin_state;
    }

    domains_to_bin_state_counts.col(i) = domain;
  }

  mat domain_centroids;

  kmeans(domain_centroids, domains_to_bin_state_counts, n_domain_states,
         static_subset, 11, false);

  // Assign the domains back to their states
  // TODO this is duplicated code with assigning bins to states
  std::vector<unsigned int> domains_to_states;

  for (int i = 0; i < n_domains; i++) {
    vec bin = domains_to_bin_state_counts.col(i);
    double min_dist = squared_euclid(bin, domain_centroids.col(0));
    int min_index = 0;

    for (int j = 1; j < n_domain_states; j++) {
      double dist = squared_euclid(bin, domain_centroids.col(j));
      if (dist < min_dist) {
        min_dist = dist;
        min_index = j;
      }
    }

    domains_to_states.push_back(min_index);
  }

  mat domain_transition_probabilities = initial_domain_transition_probabilities(
      domains_to_states, n_domain_states);

  cube bin_transition_probabilities = initial_bin_transition_probabilities(
      bins_to_bin_states, domains_to_states, n_bin_states, n_domain_states,
      domain_size);

  Emissions_Dec emission_probabilities = Emissions_Dec::emissions_from_states(
      n_bin_states, n_histone_marks, bins_to_bin_states, bin_data_m);

  mat initial_probabilities = mat(n_bin_states, n_domain_states)
                                  .fill(1.0 / (n_domain_states * n_bin_states));

  return Model(n_bin_states, n_domain_states, domain_size,
               domain_transition_probabilities, bin_transition_probabilities,
               emission_probabilities, initial_probabilities);
}

bool save_model(const Model &m, const std::string &directory) {

  bool success =

  // Save domain transitions
  m.domain_transition_probabilities.save(directory + DOMAIN_TRANSITIONS) &&
  
  // Save bin transitions
  m.bin_transition_probabilities.save(directory + BIN_TRANSITIONS) &&

  // Save emission probabilities
  m.emission_probs.emissions_probabilities_m().save(directory +
                                                    EMISSION_PROBABILITIES) &&

  // Save initial probabilities
  m.initial_probabilities.save(directory + INITIAL_PROBABILITIES);

  return success;
}

Model load_model(const std::string &directory, const int domain_size) {

  mat domain_transitions;
  domain_transitions.load(directory + DOMAIN_TRANSITIONS);

  const uword n_domain_states = domain_transitions.n_cols;

  if (domain_transitions.n_rows != n_domain_states) {
    // Fail
  }

  cube bin_transitions;
  bin_transitions.load(directory + BIN_TRANSITIONS);

  const uword n_bin_states = bin_transitions.n_cols;

  if (n_bin_states != bin_transitions.n_rows ||
      n_domain_states != bin_transitions.n_slices) {
    // Fail
  }

  mat emission_probabilities;
  emission_probabilities.load(directory + EMISSION_PROBABILITIES);

  mat initial_probabilities;
  initial_probabilities.load(directory + INITIAL_PROBABILITIES);

  // domain transitions: n_domains * n_domains
  // bin transitions: n_bins * n_bins * n_domains
  // emissions_probabilities: histone_marks * n_bins
  // initial_probabilities: n_bins, n_domains

  return Model(n_bin_states, n_domain_states, domain_size, domain_transitions,
               bin_transitions, emission_probabilities, initial_probabilities);
}
