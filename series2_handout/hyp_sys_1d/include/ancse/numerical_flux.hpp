#ifndef HYPSYS1D_NUMERICAL_FLUX_HPP
#define HYPSYS1D_NUMERICAL_FLUX_HPP

#include <memory>
#include <iostream>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/simulation_time.hpp>

/// Central flux.
/** This flux works does not depend on the model.
 * It is also unconditionally a bad choice.
 */
class CentralFlux {
  public:
    // Note: the interface for creating fluxes will give you access
    //       to the following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit CentralFlux(const std::shared_ptr<Model> &model)
        : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
                               const Eigen::VectorXd &uR) const
    {
        auto fL = model->flux(uL);
        auto fR = model->flux(uR);

        return 0.5 * (fL + fR);
    }

  private:
    std::shared_ptr<Model> model;
};


/// Lax-Friedrichs numerical flux.
/** This flux works for any model. */
//----------------FluxLFBegin----------------
class LaxFriedrichs {
  public:
    // Note: This version is a bit tricky. A numerical flux should be
    //       a function of the two trace values at the interface,
    //       i.e. what we call `uL`, `uR`.
    //       However, it requires 'dt' and 'dx'. Therefore,
    //       these need to be made available to the flux.
    //       This is one of the reasons why `SimulationTime`.
    LaxFriedrichs(const Grid &grid,
                  const std::shared_ptr<Model> &model,
                  std::shared_ptr<SimulationTime> simulation_time)
        : simulation_time(std::move(simulation_time)),
          grid(grid),
          model(model) {}

    Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
                               const Eigen::VectorXd &uR) const {
        double dx = grid.dx;
        double dt = simulation_time->dt;

        auto fL = model->flux(uL);
        auto fR = model->flux(uR);

        return 0.5*(fL+fR) - 0.5*(dx/dt)*(uR-uL);
        //return Eigen::VectorXd();
    }

  private:
    std::shared_ptr<SimulationTime> simulation_time;
    Grid grid;
    std::shared_ptr<Model> model;
};
//----------------FluxLFEnd----------------


/// Rusanov's flux (or local Lax-Friedrichs).
/** This flux works for any model. */
//----------------FluxRusanovBegin----------------
class Rusanov {
  public:
    explicit Rusanov(const std::shared_ptr<Model> &model)
        : model(model) {}

    Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
                               const Eigen::VectorXd &uR) const
    {
        auto fL = model->flux(uL);
        auto fR = model->flux(uR);
       
        double max_eigvalL = model->max_eigenvalue(uL);
        double max_eigvalR = model->max_eigenvalue(uR);
        double s = std::max(max_eigvalL,max_eigvalR);

        return 0.5 * (fL + fR) - 0.5*s*(uR-uL);
    }

  private:
    std::shared_ptr<Model> model;
};
//----------------FluxRusanovEnd----------------

/// Roe flux.
/** This requires knowledge about the model.
 *  It is also well-known for generating unphysical weak solutions.
 */
//----------------FluxRoeBegin----------------
class Roe{
  public:
    explicit Roe(const std::shared_ptr<Model> &model)
        : model(model) {}

    Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
                               const Eigen::VectorXd &uR) const
    {
        auto fL = model->flux(uL);
        auto fR = model->flux(uR);

        Eigen::VectorXd uRoe = model->roe_avg(uL,uR);

        int n_vars = model->get_nvars();
        Eigen::VectorXd uRoe_eigvals = model->eigenvalues(uRoe);

        Eigen::MatrixXd Lambda = Eigen::MatrixXd::Zero(n_vars,n_vars);
        Lambda(0,0) = std::abs(uRoe_eigvals(0));
        Lambda(1,1) = std::abs(uRoe_eigvals(1));
        Lambda(2,2) = std::abs(uRoe_eigvals(2));

        Eigen::MatrixXd R = model->eigenvectors(uRoe);

        Eigen::MatrixXd s = R*Lambda*(R.inverse());

        return 0.5 * (fL + fR) - 0.5*s*(uR-uL);
        //return Eigen::VectorXd();
    }

  private:
    std::shared_ptr<Model> model;
};
//----------------FluxRoeEnd----------------



/// HLL flux.
/** This requires knowledge about the model. */
//----------------FluxHLLBegin---------------- 
class HLL {
  public:
    explicit HLL(const std::shared_ptr<Model> &model) : model(model) {}

    Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
                               const Eigen::VectorXd &uR) const
    {
        auto fL = model->flux(uL);
        auto fR = model->flux(uR);

        Eigen::VectorXd uRoe = model->roe_avg(uL,uR);
        
        int n_vars = model->get_nvars();
        Eigen::VectorXd uL_eigvals = model->eigenvalues(uL);
        Eigen::VectorXd uR_eigvals = model->eigenvalues(uR);
        Eigen::VectorXd uRoe_eigvals = model->eigenvalues(uRoe);

        double sL = std::numeric_limits<double>::max();
        for (int i=0;i<n_vars;++i){
            sL = std::min(sL,uL_eigvals(i));
            sL = std::min(sL,uRoe_eigvals(i));
        }

        double sR = std::numeric_limits<double>::lowest();
        for (int i=0;i<n_vars;++i){
            sR = std::max(sR,uR_eigvals(i));
            sR = std::max(sR,uRoe_eigvals(i));
        }

        Eigen::VectorXd result;
        if (sL>=0) {
            result = fL;
        }
        else if (sL<0 && sR>0) {
            result = (sR*fL-sL*fR+sR*sL*(uR-uL))/(sR-sL);
        }
        else if (sR<=0) {
            result = fR;
        }

        return result;
        //return Eigen::VectorXd();
    }

  private:
    std::shared_ptr<Model> model;
};
//----------------FluxHLLEnd---------------- 

/// HLLC flux
/** This requires knowledge about the model.
 *  This version is for the Euler equation.
 */
//----------------FluxHLLCEulerBegin----------------  
class HLLCEuler {
  public:
    explicit HLLCEuler(const std::shared_ptr<Euler> &model)
        : model(model) {
        n_vars = model->get_nvars();
    }

    Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
                               const Eigen::VectorXd &uR) const
    {
        auto fL = model->flux(uL);
        auto fR = model->flux(uR);

        Eigen::VectorXd uRoe = model->roe_avg(uL,uR);
        
        int n_vars = model->get_nvars();
        Eigen::VectorXd uL_eigvals = model->eigenvalues(uL);
        Eigen::VectorXd uR_eigvals = model->eigenvalues(uR);
        Eigen::VectorXd uRoe_eigvals = model->eigenvalues(uRoe);

        double sL = std::numeric_limits<double>::max();
        for (int i=0;i<n_vars;++i){
            sL = std::min(sL,uL_eigvals(i));
            sL = std::min(sL,uRoe_eigvals(i));
        }

        double sR = std::numeric_limits<double>::lowest();
        for (int i=0;i<n_vars;++i){
            sR = std::max(sR,uR_eigvals(i));
            sR = std::max(sR,uRoe_eigvals(i));
        }

        double vL = uL(1)/uL(0);
        double vR = uR(1)/uR(0);

        double pL = model->pressure(uL(0),vL,uL(2));
        double pR = model->pressure(uR(0),vR,uR(2));
        //9.33
        double sM = uR(1)*(sR-vR)-uL(1)*(sL-vL)-(pR-pL);
        sM = sM/(uR(0)*(sR-vR)-uL(0)*(sL-vL));

        //9.31
        double rho_intL = uL(0)*(vL-sL)/(sM-sL);
        double rho_intR = uR(0)*(vR-sR)/(sM-sR);

        //9.34
        double p_int = pL+uL(0)*(vL-sM)*(vL-sL);

        double gamma = model->get_gamma();

        Eigen::VectorXd u_intL(n_vars);
        u_intL(0) = rho_intL;
        u_intL(1) = rho_intL*sM;
        u_intL(2) = p_int/(gamma-1) + 0.5*rho_intL*sM*sM;

        Eigen::VectorXd u_intR(n_vars);
        u_intR(0) = rho_intR;
        u_intR(1) = rho_intR*sM;
        u_intR(2) = p_int/(gamma-1) + 0.5*rho_intR*sM*sM;

        Eigen::VectorXd f_intL = fL + sL*(u_intL-uL);
        Eigen::VectorXd f_intR = fR + sR*(u_intR-uR);

        Eigen::VectorXd result;
        if (sL>0) {
            result = fL;
        }
        else if (sL<=0 && sM>0) {
            result = f_intL;
        }
        else if (sM<=0 && sR>0) {
            result = f_intR;
        }
        else {
            result = fR;
        }

        return result;       
        //return Eigen::VectorXd();
    }

  private:
    std::shared_ptr<Euler> model;
    int n_vars;
};
//----------------FluxHLLCEulerEnd----------------  



#endif // HYPSYS1D_NUMERICAL_FLUX_HPP
