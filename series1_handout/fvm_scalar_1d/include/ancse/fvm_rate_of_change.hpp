#ifndef FVMSCALAR1D_FVM_RATE_OF_CHANGE_HPP
#define FVMSCALAR1D_FVM_RATE_OF_CHANGE_HPP

#include <iostream>

#include <Eigen/Dense>
#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/rate_of_change.hpp>
#include <ancse/simulation_time.hpp>
#include <memory>

//----------------FVMRateOfChangeBegin----------------
/// Compute the rate of change due to FVM.
/** The semidiscrete approximation of a PDE using FVM is
 *      du_i/dt = - (F_{i+0.5} - F_{i-0.5}) / dx.
 *  This computes the right hand side of the ODE.
 *
 * @tparam NumericalFlux see e.g. `CentralFlux`.
 * @tparam Reconstruction see e.g. `PWConstantReconstruction`.
 */
template <class NumericalFlux, class Reconstruction>
class FVMRateOfChange : public RateOfChange {
  public:
    FVMRateOfChange(const Grid &grid,
                    const NumericalFlux &numerical_flux,
                    const Reconstruction &reconstruction)
        : grid(grid),
          numerical_flux(numerical_flux),
          reconstruction(reconstruction) {}

    virtual void operator()(Eigen::VectorXd &dudt,
                            const Eigen::VectorXd &u0) const override {
        // implement the flux loop here.
	//int n_cells = dudt.size();
	int n_cells = grid.n_cells;
	double dx = grid.dx;
	double FL,FR;
	double uLL,uLR,uRL,uRR;

	for (int j=2;j<n_cells-2;j++) {
	    //std::cout<<"j: "<<j<<std::endl;
            //std::cout<<"u0 size: "<<u0.size()<<std::endl;
	    //std::cout<<u0[0]<<std::endl;
	    //std::cout<<u0[1]<<std::endl;
	    auto recL = reconstruction(u0,j-1);
	    //std::cout<<"FVMRateOfChange () operator just after first rec call"<<std::endl;
	    auto recR = reconstruction(u0,j);

	    uLL = std::get<0>(recL);
	    uLR = std::get<1>(recL);
	    uRL = std::get<0>(recR);
	    uRR = std::get<1>(recR);

	    FL = numerical_flux(uLL,uLR); //-1/2
	    //std::cout<<"FVMRateOfChange () operator just after first nf call"<<std::endl;
	    FR = numerical_flux(uRL,uRR); //+1/2
	    dudt[j] = -1.0*(FR-FL)/dx;
	}
    }

  private:
    Grid grid;
    NumericalFlux numerical_flux;
    Reconstruction reconstruction;
};
//----------------FVMRateOfChangeEnd----------------

std::shared_ptr<RateOfChange>
make_fvm_rate_of_change(const Grid &grid,
                        const Model &model,
                        const std::shared_ptr<SimulationTime> &simulation_time);

#endif // FVMSCALAR1D_FVM_RATE_OF_CHANGE_HPP
