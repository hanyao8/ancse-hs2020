#pragma once
#include <Eigen/Core>
#include "load_vector.hpp"

//----------------AssembleVectorBegin----------------
//! Assemble the load vector into the full right hand side
//! for the linear system
//!
//! @param[out] F will at the end contain the RHS values for each vertex.
//! @param[in] vertices a list of triangle vertices
//! @param[in] triangles a list of triangles
//! @param[in] f the RHS function f.
void assembleLoadVector(Eigen::VectorXd& F,
                           const Eigen::MatrixXd& vertices,
                           const Eigen::MatrixXi& triangles,
                           const std::function<double(double, double)>& f)
{
     const int numberOfElements = triangles.rows();

     F.resize(vertices.rows());
     F.setZero();
     //// ANCSE_START_TEMPLATE
     for (int i = 0; i < numberOfElements; ++i) {
         const auto& indexSet = triangles.row(i);

         const auto& a = vertices.row(indexSet(0));
         const auto& b = vertices.row(indexSet(1));
         const auto& c = vertices.row(indexSet(2));

         Eigen::Vector3d elementVector;
         computeLoadVector(elementVector, a, b, c, f);

         for (int i = 0; i < 3; ++i) {
            F(indexSet(i)) += elementVector(i);
         }
     }
     //// ANCSE_END_TEMPLATE
}
//----------------AssembleVectorEnd----------------
