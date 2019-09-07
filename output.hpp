#ifndef EV_OUTPUT_HPP
#define EV_OUTPUT_HPP

#include <fstream>
#include <string>
#include <vector>

#include "motherbase.hpp"
#include "matrix.hpp"

class Output : protected Motherbase
{
private:

public:

  Output(DSMC*);
  ~Output() = default;

  // Testing outputs
  void output_kernel(void);
  void output_weights(void);
  void output_forces(void);
  void output_density(void);
  void output_majorants(void);
  void output_collisions(void);

  // Output sample
  template <class data_type>
  void output_sample(ev_matrix::MaskMatrix<data_type>& sample, const DefaultString& file_name)
  {
    std::ofstream file1(file_name);
    file1 << sample;
    file1.close();
  }

  // Output sample, with label
  template <class data_type, class tag_type>
  void output_sample(ev_matrix::MaskMatrix<data_type>& sample, const DefaultString& file_name, tag_type label)
  {
    DefaultString file_name_tag = file_name + "_t=" + std::to_string(label) + ".txt";
    std::ofstream file1(file_name_tag);
    file1 << sample;
    file1.close();
  }

  // Output vector
  template <class data_type>
  void output_vector(std::vector<data_type>& sample, const DefaultString& file_name)
  {
    std::ofstream file1(file_name);
    for (auto it = sample.cbegin(); it!=sample.cend(); ++it)
      file1 << *it << "\n";
    file1.close();
  }

  // Output collisions statistics
  void output_collisions_stat(void);

  // Output function form vectors
  void output_fun_vec(const std::vector<real_number>&,
    const std::vector<real_number>&);

};

#endif /* EV_OUTPUT_HPP */
