#ifndef FVMSCALAR1D_NUMERICAL_FLUX_HPP
#define FVMSCALAR1D_NUMERICAL_FLUX_HPP

#include <memory>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/simulation_time.hpp>

#include <cmath>
#include <algorithm>

/// Central flux.
/** This flux works does not depend on the model. It is also unconditionally a
 * bad choice.
 */
class CentralFlux {
  public:
    // Note: the interface for creating fluxes will give you access to the
    //       following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit CentralFlux(const Model &model) : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    double operator()(double uL, double uR) const {
        auto fL = model.flux(uL);
        auto fR = model.flux(uR);

        return 0.5 * (fL + fR);
    }

  private:
    Model model;
};


//----------------FluxLFBegin----------------
/// Lax-Friedrichs numerical flux.
/** This flux works for any model. */
class LaxFriedrichs {
  public:
    // Note: This version is a bit tricky. A numerical flux should be
    //       a function of the two trace values at the interface, i.e. what we
    //       call `uL`, `uR`. However, it requires 'dt' and 'dx'. Therefore,
    //       these need to be made available to the flux. This is one of the
    //       reasons why `SimulationTime`.
    LaxFriedrichs(const Grid &grid,
                  const Model &model,
                  std::shared_ptr<SimulationTime> simulation_time)
        : simulation_time(std::move(simulation_time)),
          grid(grid),
          model(model) {}

    double operator()(double uL, double uR) const {
        double dx = grid.dx;
        double dt = simulation_time->dt;

        // Implement the LaxFriedrichs flux
	auto fL = model.flux(uL);
	auto fR = model.flux(uR);
	double Flxf = 0.5*(fL + fR) - 0.5*(dx/dt)*(uR-uL);
        return Flxf;

    }

  private:
    std::shared_ptr<SimulationTime> simulation_time;
    Grid grid;
    Model model;
};
//----------------FluxLFEnd----------------

class Rusanov {
  public:
    Rusanov(const Grid &grid,
                  const Model &model,
                  std::shared_ptr<SimulationTime> simulation_time)
        : simulation_time(std::move(simulation_time)),
          grid(grid),
          model(model) {}

    double operator()(double uL, double uR) const {
        double dx = grid.dx;
        double dt = simulation_time->dt;

        // Implement the LaxFriedrichs flux
	auto fL = model.flux(uL);
	auto fR = model.flux(uR);

	auto sL = model.dflux_du(uL);
	auto sR = model.dflux_du(uR);

        auto sM = std::max(sL,sR); //speed middle
	double Frus = 0.5*(fL + fR) - 0.5*sM*(uR-uL);
	return Frus;

    }

  private:
    std::shared_ptr<SimulationTime> simulation_time;
    Grid grid;
    Model model;
};


class Roe {
  public:
    Roe(const Grid &grid,
                  const Model &model,
                  std::shared_ptr<SimulationTime> simulation_time)
        : simulation_time(std::move(simulation_time)),
          grid(grid),
          model(model) {}

    double operator()(double uL, double uR) const {
        // Implement the LaxFriedrichs flux
	auto fL = model.flux(uL);
	auto fR = model.flux(uR);

	auto AM;
	if (uL!=uR) {
	    AM = (model.flux(uR)-model.flux(uL))/(uR-uL);
	}
	else {
	    AM = model.dflux_du(uL);
	}

	auto Froe;
	if (AM>=0) {
	    Froe = model.flux(uL);
	}
	else {
            Froe = model.flux(uR);
	}
	return Froe;
    }

  private:
    std::shared_ptr<SimulationTime> simulation_time;
    Grid grid;
    Model model;
};

class Godunov {
  public:
    Godunov(const Grid &grid,
                  const Model &model,
                  std::shared_ptr<SimulationTime> simulation_time)
        : simulation_time(std::move(simulation_time)),
          grid(grid),
          model(model) {}

    double operator()(double uL, double uR) const {

	double omega = model.flux_omega();
	auto f1 = model.flux(std::max(uL,omega))
	auto f2 = model.flux(std::min(uR,omega))
	double Fgod = std::max(f1,f2)
	return Fgod;
    }

  private:
    std::shared_ptr<SimulationTime> simulation_time;
    Grid grid;
    Model model;
};


class EngquistOsher {
  public:
    EngquistOsher(const Grid &grid,
                  const Model &model,
                  std::shared_ptr<SimulationTime> simulation_time)
        : simulation_time(std::move(simulation_time)),
          grid(grid),
          model(model) {}

    double operator()(double uL, double uR) const {

	double omega = model.flux_omega();
	auto f1 = model.flux(std::max(uL,omega))
	auto f2 = model.flux(std::min(uR,omega))
	double Feo = f1 + f2;
	return Feo;
    }

  private:
    std::shared_ptr<SimulationTime> simulation_time;
    Grid grid;
    Model model;
};

#endif // FVMSCALAR1D_NUMERICAL_FLUX_HPP
