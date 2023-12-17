//
// Created by Stephanos Tsoucas on 6/14/17.
//

#ifndef BOOSTPYTHONHELLOWORLD_MODEL_H
#define BOOSTPYTHONHELLOWORLD_MODEL_H

#define ARMA_NO_DEBUG true

#include "Emissions_Dec.h"
#include <armadillo>
#include <stdio.h>

using namespace arma;

const std::string BIN_TRANSITIONS = "bin_transitions.mat";
const std::string DOMAIN_TRANSITIONS = "domain_transitions.mat";
const std::string EMISSION_PROBABILITIES = "emission_probabilities.mat";
const std::string INITIAL_PROBABILITIES = "initial_probabilities.mat";

struct Model {
  const unsigned int n_bin_states;
  const unsigned int n_domain_states;
  const unsigned int domain_size;
  // domain_transition_probabilities(from, to)
  const mat domain_transition_probabilities;
  // bin_transition_probabilities(from, to, dest_domain)
  const cube bin_transition_probabilities;

  // TODO: change to emissions_dec
  const Emissions_Dec emission_probs;
  // initial_probabilities(bin_state, domain_state)
  const mat initial_probabilities;

  const uword n_histone_marks;

  double emission_probability(const int bin_state,
                              const uword observation) const {
    return emission_probs.emission_probability(bin_state, observation);
  }

  rowvec emission_probabilities_row(const uword observation) const {
    return emission_probs.emission_probabilities_row(observation);
  }

  colvec emission_probabilities_col(const uword observation) const {
    return emission_probs.emission_probabilities_col(observation);
  }

  cube bin_transition_probabilities_t() const {
    cube output(bin_transition_probabilities.n_cols,
                bin_transition_probabilities.n_rows,
                bin_transition_probabilities.n_slices);

    for (int i = 0; i < bin_transition_probabilities.n_slices; i++) {
      output.slice(i) = bin_transition_probabilities.slice(i).t();
    }

    return output;
  }

  Model(const unsigned int n_bin_states, const unsigned int n_domain_states,
        const unsigned int domain_size, mat domain_transition_probabilities,
        cube bin_transition_probabilities, Emissions_Dec emission_probabilities,
        mat initial_probabilities)
      : n_bin_states(n_bin_states), n_domain_states(n_domain_states),
        domain_size(domain_size),
        domain_transition_probabilities(domain_transition_probabilities),
        bin_transition_probabilities(bin_transition_probabilities),
        emission_probs(emission_probabilities),
        initial_probabilities(initial_probabilities),
        n_histone_marks(
            emission_probabilities.emissions_probabilities_m().n_rows) {}
};

Model initial_model(const mat &bin_data, unsigned int n_bin_states,
                    unsigned int n_domain_states, unsigned int n_histone_marks,
                    unsigned int domain_size);

Model load_model(const std::string &directory, const int domain_size);

bool save_model(const Model &m, const std::string &directory);

#endif // BOOSTPYTHONHELLOWORLD_MODEL_H
