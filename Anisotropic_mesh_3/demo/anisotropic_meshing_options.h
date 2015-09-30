#ifndef ANISOTROPIC_MESHING_OPTIONS_H
#define ANISOTROPIC_MESHING_OPTIONS_H

enum Domain_type   { POLYHEDRAL_SURFACE = 0, 
                     IMPLICIT_SURFACE 
                   };

enum Metric_options{ EUCLIDEAN = 0
                   , POLYHEDRON_CURVATURE 
                   , IMPLICIT_CURVATURE
                   , HYPERBOLIC_SHOCK
                   , EXTERNAL_MF
                   };

#endif
