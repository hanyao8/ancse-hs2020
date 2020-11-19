#include <ancse/cfl_condition.hpp>


//----------------StandardCFLConditionDefnBegin----------------
StandardCFLCondition::StandardCFLCondition(const Grid &grid,
                                           const Model &model,
                                           double cfl_number)
    : grid(grid), model(model), cfl_number(cfl_number) {}

double StandardCFLCondition::operator()(const Eigen::VectorXd &u) const {

    auto n_cells = grid.n_cells;
    auto n_ghost = grid.n_ghost;

    auto dx = grid.dx;
    double dt_bound;

    double max_speed = 0.0;
    double tmp_speed = 0.0;
    for (int j=0;j<n_cells;j++) {
	tmp_speed = std::abs(model.dflux_du(u[j]));
        if (tmp_speed > max_speed) {
	    max_speed = tmp_speed;
	}
    }

    dt_bound = cfl_number * dx * max_speed;

    return dt_bound;
}
//----------------StandardCFLConditionDefnEnd----------------


std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid, const Model &model, double cfl_number) {
    // implement this 'factory' for your CFL condition.
    return std::make_shared<StandardCFLCondition>(grid, model, cfl_number);
    return nullptr;
}
