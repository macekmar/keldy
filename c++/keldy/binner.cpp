#include "binner.hpp"

namespace keldy::binner {

void mpi_broadcast(continuous_axis_t &in, mpi::communicator c, int root) {
  mpi::broadcast(in.xmin, c, root);
  mpi::broadcast(in.xmax, c, root);
  mpi::broadcast(in.nr_bins, c, root);
  mpi::broadcast(in.bin_size, c, root);
};

void mpi_broadcast(discreet_axis_t &in, mpi::communicator c, int root) { mpi::broadcast(in.nr_bins, c, root); };

bool operator==(continuous_axis_t const &a, continuous_axis_t const &b) {
  return (a.xmin == b.xmin) && (a.xmax == b.xmax) && (a.nr_bins == b.nr_bins) && (a.bin_size == b.bin_size);
};

}; // namespace keldy::binner
