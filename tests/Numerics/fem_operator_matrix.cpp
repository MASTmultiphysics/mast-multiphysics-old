/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Fixtures
#include <boost/test/unit_test.hpp>

struct Massive
{
    int m;
    
    Massive() : m(2)
    {
        BOOST_TEST_MESSAGE("setup mass");
    }
    
    ~Massive()
    {
        BOOST_TEST_MESSAGE("teardown mass");
    }
};

BOOST_FIXTURE_TEST_SUITE(Physics, Massive)

BOOST_AUTO_TEST_CASE(specialTheory)
{
    int e = 32;
    int c = 4;
    
    BOOST_CHECK(e == m * c * c);
}

BOOST_AUTO_TEST_CASE(newton2)
{
    int f = 10;
    int a = 5;
    
    BOOST_CHECK(f == m * a);
}

BOOST_AUTO_TEST_SUITE_END()


/*#include "Numerics/fem_operator_matrix.h"

int main_fem_operator (int argc, char* const argv[])
{
    DenseRealVector vec1, vec2, shp;
    DenseRealMatrix mat1, mat2, bmat1, bmat2, tmp;
    
    vec1.resize(12); shp.resize(4);
    
    for (unsigned int i=0; i<shp.size(); i++)
        shp(i) = i+1;
    
    for (unsigned int i=0; i<vec1.size(); i++)
        vec1(i) = i+100;
    
    shp.print(libMesh::out);
    vec1.print(libMesh::out);
    mat1.print(libMesh::out);
    
    FEMOperatorMatrix b1, b2;
    b1.reinit(3, shp); b2.reinit(3, shp);
    bmat1.resize(3, shp.size()*3); bmat2.resize(3, shp.size()*3);
    for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<shp.size(); j++) {
            bmat1(i, i*shp.size()+j) = shp(j);
            bmat2(i, i*shp.size()+j) = shp(j);
        }
    
    
    vec2.resize(3);
    b1.vector_mult(vec2, vec1);
    libMesh::out << "FEMOperator: " << std::endl; vec2.print(libMesh::out);
    bmat1.vector_mult(vec2, vec1);
    libMesh::out << "Matrix:      " << std::endl; vec2.print(libMesh::out);
    
    vec1.resize(3);
    for (unsigned int i=0; i<vec1.size(); i++)
        vec1(i) = 100+i;
    vec2.resize(12);
    b1.vector_mult_transpose(vec2, vec1);
    libMesh::out << "FEMOperator: " << std::endl; vec2.print(libMesh::out);
    bmat1.vector_mult_transpose(vec2, vec1);
    libMesh::out << "Matrix:      " << std::endl; vec2.print(libMesh::out);
    
    
    mat1.resize(5,3);
    mat2.resize(5, 12);
    for (unsigned int i=0; i<5; i++)
        for (unsigned int j=0; j<3; j++)
            mat1(i,j) = (i+1)*(j+1);
    libMesh::out << "multiply matrix" << std::endl;  mat1.print();
    b1.left_multiply(mat2, mat1);
    libMesh::out << "FEMOperator: " << std::endl;  mat2.print();
    tmp = bmat1;
    tmp.left_multiply(mat1);
    libMesh::out << "Matrix:      " << std::endl;  tmp.print();
    
    
    mat1.resize(5,12);
    mat2.resize(5, 3);
    for (unsigned int i=0; i<5; i++)
        for (unsigned int j=0; j<12; j++)
            mat1(i,j) = (i+1)*(j+1);
    libMesh::out << "multiply matrix" << std::endl;  mat1.print();
    b1.left_multiply_transpose(mat2, mat1);
    libMesh::out << "FEMOperator: " << std::endl;  mat2.print();
    bmat1.get_transpose(tmp);
    tmp.left_multiply(mat1);
    libMesh::out << "Matrix:      " << std::endl;  tmp.print();
    
    
    mat1.resize(12,5);
    mat2.resize(3, 5);
    for (unsigned int i=0; i<12; i++)
        for (unsigned int j=0; j<5; j++)
            mat1(i,j) = (i+1)*(j+1);
    libMesh::out << "multiply matrix" << std::endl;  mat1.print();
    b1.right_multiply(mat2, mat1);
    libMesh::out << "FEMOperator: " << std::endl;  mat2.print();
    tmp = bmat1;
    tmp.right_multiply(mat1);
    libMesh::out << "Matrix:      " << std::endl;  tmp.print();
    
    
    mat1.resize(3, 5);
    mat2.resize(12, 5);
    for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<5; j++)
            mat1(i,j) = (i+1)*(j+1);
    libMesh::out << "multiply matrix" << std::endl;  mat1.print();
    b1.right_multiply_transpose(mat2, mat1);
    libMesh::out << "FEMOperator: " << std::endl;  mat2.print();
    bmat1.get_transpose(tmp);
    tmp.right_multiply(mat1);
    libMesh::out << "Matrix:      " << std::endl;  tmp.print();
    
    
    mat2.resize(12, 12);
    b1.right_multiply_transpose(mat2, b1);
    libMesh::out << "FEMOperator: " << std::endl;  mat2.print();
    bmat1.get_transpose(tmp);
    tmp.right_multiply(bmat1);
    libMesh::out << "Matrix:      " << std::endl;  tmp.print();
    
    return 0;
}
*/

