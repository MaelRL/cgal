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
// This class is used to avoid repetitious facets while outputing.
// For example, facet<1,3,2> and facet<1,2,3> are the same facet, 
// we sorted the identities of the vertices of a facet and stored them
// in a sorted list.

#ifndef CGAL_ANISOTROPIC_MESH_3_OUTPUT_FACETS_H
#define CGAL_ANISOTROPIC_MESH_3_OUTPUT_FACETS_H

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>

namespace CGAL{
namespace Anisotropic_mesh_3{

class Output_facets {
public:
    struct Facet {
        int vertices[3];
        Facet(int a, int b, int c) {
            vertices[0] = a;
            vertices[1] = b;
            vertices[2] = c;
            for (int i = 2; i >= 0; i--)
                for (int j = 1; j <= i; j++)
                    if (vertices[j] < vertices[j - 1])
                        std::swap(vertices[j], vertices[j - 1]);
        }
        bool operator<(const Facet &other) const {
            for (int i = 0; i < 3; i++) {
                if (vertices[i] < other.vertices[i])
                    return true;
                else if (vertices[i] > other.vertices[i])
                    return false;
            }
            return false;
        }
        bool operator==(const Facet &other) const {
            for (int i = 0; i < 3; i++)
                if (vertices[i] != other.vertices[i])
                    return false;
            return true;
        }
    };

public:
	typedef std::vector<Facet> FacetList;
    typedef FacetList::iterator Facet_handle;

protected:
    FacetList facets;

protected:
    bool insert(const Facet &facet) {
        if (facets.size() == 0) {
            facets.push_back(facet);
            return true;
        }
        int lp = 0, rp = (int)facets.size() - 1;
        int mp = (lp + rp) / 2;
        while (rp - lp > 1) {
            if (facets[mp] == facet)
                return false;
            else if (facets[mp] < facet)
                lp = mp;
            else
                rp = mp;
            mp = (lp + rp) / 2;
        }
        if (facets[lp] == facet)
            return false;
        if (facets[rp] == facet)
            return false;
        if (rp == lp)
            if (facets[lp] < facet)
                facets.insert(facets.begin() + lp + 1, facet);
            else
                facets.insert(facets.begin() + lp, facet);
        else
            if (facet < facets[lp])
                facets.insert(facets.begin() + lp, facet);
            else if (facets[rp] < facet)
                facets.insert(facets.begin() + rp + 1, facet);
            else
                facets.insert(facets.begin() + rp, facet);
        return true;
    }

public:
    Facet_handle begin() { return facets.begin(); }
    Facet_handle end() { return facets.end(); }
    int size() { return (int)facets.size(); }
    bool insert(int a, int b, int c) {
        return insert(Facet(a, b, c));
    }
    
    Output_facets() : facets() { }
    ~Output_facets() { }
};
}
}

#endif // CGAL_ANISOTROPIC_MESH_3_OUTPUT_FACETS_H
