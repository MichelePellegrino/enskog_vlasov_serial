#ifndef EV_OUTPUT_HPP
#define EV_OUTPUT_HPP

#include <fstream>
#include <string>

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
  void output_forces(void);
  void output_density(void);
  void output_majorants(void);
  void output_collisions(void);

  // Output sample
  template <class data_type>
  void output_sample(ev_matrix::MaskMatrix<data_type>& sample, const std::string& file_name)
  {
    std::ofstream file1(file_name);
    file1 << sample;
    file1.close();
  }

  // Output fucntion form vectors
  void output_fun_vec(const std::vector<real_number>&,
    const std::vector<real_number>&);

};

#endif /* EV_OUTPUT_HPP */
