// Copyright (c) 2011  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s) : Kan-Le Shi
//
// This class is used to avoid repetitious cells while outputing.
// For example, cell<1,3,2,4> and cell<1,2,3,4> are the same cell, 
// we sorted the identities of the vertices of a cell and stored them
// in a sorted list.


#ifndef CGAL_ANISOTROPIC_MESH_3_OUTPUT_CELLS_H
#define CGAL_ANISOTROPIC_MESH_3_OUTPUT_CELLS_H

#include <iostream>
#include <fstream>
#include <utility>


#include "utility"
#include "vector"

namespace CGAL{
namespace Anisotropic_mesh_3{

class Output_cells {
public:
    struct Cell {
        int vertices[4];
        Cell(int a, int b, int c, int d) 
        {
            vertices[0] = a;
            vertices[1] = b;
            vertices[2] = c;
            vertices[3] = d;	
            for (int i = 3; i >= 0; i--)
                for (int j = 1; j <= i; j++)
                    if (vertices[j] < vertices[j - 1])
                        std::swap(vertices[j], vertices[j - 1]);
        }
        bool operator<(const Cell &other) const 
        {
            for (int i = 0; i < 4; i++) {
                if (vertices[i] < other.vertices[i])
                    return true;
                else if (vertices[i] > other.vertices[i])
                    return false;
            }
            return false;
        }
        bool operator==(const Cell &other) const 
        {
            for (int i = 0; i < 4; i++)
                if (vertices[i] != other.vertices[i])
                    return false;
            return true;
        }
    };

public:
	typedef std::vector<Cell> CellList;
    typedef CellList::iterator Cell_handle;

protected:
    CellList cells;

protected:
    bool insert(const Cell &cell) 
    {
        if (cells.size() == 0) {
            cells.push_back(cell);
            return true;
        }
        int lp = 0, rp = (int)cells.size() - 1;
        int mp = (lp + rp) / 2;
        while (rp - lp > 1) {
            if (cells[mp] == cell)
                return false;
            else if (cells[mp] < cell)
                lp = mp;
            else
                rp = mp;
            mp = (lp + rp) / 2;
        }
        if (cells[lp] == cell)
            return false;
        if (cells[rp] == cell)
            return false;
        if (rp == lp)
            if (cells[lp] < cell)
                cells.insert(cells.begin() + lp + 1, cell);
            else
                cells.insert(cells.begin() + lp, cell);
        else
            if (cell < cells[lp])
                cells.insert(cells.begin() + lp, cell);
            else if (cells[rp] < cell)
                cells.insert(cells.begin() + rp + 1, cell);
            else
                cells.insert(cells.begin() + rp, cell);
        return true;
    }

public:
    Cell_handle begin() { return cells.begin(); }
    Cell_handle end() { return cells.end(); }
    int size() { return (int)cells.size(); }
    bool insert(int a, int b, int c, int d) {
        return insert(Cell(a, b, c, d));
    }
    
    Output_cells() : cells() { }
    ~Output_cells() { }
};
}
}

#endif // CGAL_ANISOTROPIC_MESH_3_OUTPUT_CELLS_H