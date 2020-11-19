#ifndef FVMSCALAR1D_RECONSTRUCTION_HPP
#define FVMSCALAR1D_RECONSTRUCTION_HPP

#include <Eigen/Dense>
#include <cmath>
#include <memory>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/rate_of_change.hpp>
#include <ancse/simulation_time.hpp>

inline double sign(double a) { return copysign(1.0, a); }


//----------------SlopeLimiterABegin----------------
inline double minmod(double a, double b) {
    return 0.5 * (sign(a) + sign(b)) * std::min(std::abs(a), std::abs(b));
}

inline double maxmod(double a, double b) {
    return 0.5 * (sign(a) + sign(b)) * std::max(std::abs(a), std::abs(b));
}

inline double minabs(double a, double b) {
    double result;
    if (std::abs(a) < std::abs(b)) {
        result = a;
    }
    else if (std::abs(a) > std::abs(b)) {
        result = b;
    }
    else {
        result = 0.5*(a+b);
    }
    return result;
}

inline double minmod3(double a, double b, double c) {
    double result;
    double epsilon = 0.001;
    if ( sign(a)+sign(b)+sign(c)+epsilon > 3.0 ) {
        result = std::min(std::abs(a),std::abs(b));
	result = 1.0 * std::min(result,std::abs(c));
    }
    else if ( sign(a)+sign(b)+sign(c)-epsilon < 3.0 ) {
        result = std::min(std::abs(a),std::abs(b));
	result = -1.0 * std::min(result,std::abs(c));
    }
    else {
        result = 0.0;
    }
    return result;
}

inline double vanleer(double a, double b) {
    double w1 = std::abs(b)/(std::abs(b)+std::abs(a));
    double w2 = std::abs(a)/(std::abs(b)+std::abs(a));
    return w1*a+w2*b;
}


//----------------SlopeLimiterAEnd----------------

//----------------SlopeLimiterBBegin----------------
struct MinMod {
    inline double operator()(double sL, double sR) const {
        return minmod(sL, sR);
    }
};

struct MinAbs {
    inline double operator()(double sL, double sR) const {
        return minabs(sL, sR);
    }
};

struct SuperBee {
    inline double operator()(double sL, double sR) const {
	double minmodL = minmod(2.0*sL, 1.0*sR);
	double minmodR = minmod(1.0*sL, 2.0*sR);
        return maxmod(minmodL, minmodR);
    }
};

struct MonotonizedCentral {
    inline double operator()(double sL, double sR) const {
        return minmod3(2.0*sR,0.5*(sL+sR),2.0*sL);
    }
};

struct VanLeer {
    inline double operator()(double sL, double sR) const {
        return vanleer(sL,sR);
    }
};

//----------------SlopeLimiterBEnd----------------


class PWConstantReconstruction {
  public:
    /// Compute the left and right trace at the interface i + 1/2.
    /** Note: This API is agnostic to the number of cell-averages required
     *        by the method. Therefore, reconstructions with different stencil
     *        sizes can implement this API; and this call can be used in parts
     *        of the code that do not need to know about the details of the
     *        reconstruction.
     */
    inline std::pair<double, double> operator()(const Eigen::VectorXd &u,
                                                int i) const {
        return (*this)(u[i], u[i + 1]);
    }

    /// Compute the left and right trace at the interface.
    /** Piecewise constant reconstruction of the left and right trace only
     *  requires the cell-average to the left and right of the interface.
     *
     *  Note: Compared to the other overload this reduces the assumption on
     *        how the cell-averages are stored. This is useful when testing and
     *        generally makes the function useful in more situations.
     */
    inline std::pair<double, double> operator()(double ua, double ub) const {
        return {ua, ub};
    }
};


//----------------LinearRCBegin----------------
template <class SlopeLimiter>
class PWLinearReconstruction {
  public:
    explicit PWLinearReconstruction(
		    const SlopeLimiter &slope_limiter)
        : slope_limiter(slope_limiter) {}

    std::pair<double, double> operator()(const Eigen::VectorXd &u,
                                         int i) const {
        return (*this)(u[i - 1], u[i], u[i + 1], u[i + 2]);
    }

    std::pair<double, double>
    operator()(double ua, double ub, double uc, double ud) const {

        //double uL = 0.0;
        //double uR = 0.0;

	auto sb = ub-ua; //j
	auto sc = uc-ub; //j+1
	auto sd = ud-uc; //j+2

	auto sigmaL = slope_limiter(sc,sb); //j
	auto sigmaR = slope_limiter(sd,sc); //j+1

	double uL = ub + 0.5*sigmaL; //j
	double uR = uc - 0.5*sigmaR; //j+1


        return {uL, uR};
    }

  private:
    SlopeLimiter slope_limiter;
};
//----------------LinearRCEnd----------------


#endif // FVMSCALAR1D_RATE_OF_CHANGE_HPP
