#pragma once
#include <Eigen/Core>
#include <string>
#include <igl/readMESH.h>
#include <igl/readSTL.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include "stiffness_matrix_assembly.hpp"
#include "load_vector_assembly.hpp"
#include "dirichlet_boundary.hpp"

typedef Eigen::VectorXd Vector;

//----------------solveBegin----------------
//! Solve the FEM system.
//!
//! @param[out] u will at the end contain the FEM solution.
//! @param[in] vertices list of triangle vertices for the mesh
//! @param[in] triangles list of triangles (described by indices)
//! @param[in] f the RHS f (as in the exercise)
//! return number of degrees of freedom (without the boundary dofs)
int solveFiniteElement(Vector& u,
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& triangles,
    const std::function<double(double, double)>& f)
{
    SparseMatrix A;
    //// ANCSE_START_TEMPLATE
    assembleStiffnessMatrix(A, vertices, triangles);
    //// ANCSE_END_TEMPLATE

    Vector F;
    //// ANCSE_START_TEMPLATE
    assembleLoadVector(F, vertices, triangles, f);
    //// ANCSE_END_TEMPLATE

    u.resize(vertices.rows());
    u.setZero();
    Eigen::VectorXi interiorVertexIndices;

    auto zerobc = [](double x, double y){ return 0;};
    // set homogeneous Dirichlet Boundary conditions
    //// ANCSE_START_TEMPLATE
    setDirichletBoundary(u, interiorVertexIndices, vertices, triangles, zerobc);
    F -= A * u;
    //// ANCSE_END_TEMPLATE

    SparseMatrix AInterior;

    igl::slice(A, interiorVertexIndices, interiorVertexIndices, AInterior);
    Eigen::SimplicialLDLT<SparseMatrix> solver;

    Vector FInterior;

    igl::slice(F, interiorVertexIndices, FInterior);

    //initialize solver for AInterior
    //// ANCSE_START_TEMPLATE
    solver.compute(AInterior);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Could not decompose the matrix");
    }
    //// ANCSE_END_TEMPLATE

    //solve interior system
    //// ANCSE_START_TEMPLATE
    Vector uInterior = solver.solve(FInterior);
    igl::slice_into(uInterior, interiorVertexIndices, u);
    //// ANCSE_END_TEMPLATE

    return interiorVertexIndices.size();

}
//----------------solveEnd----------------
