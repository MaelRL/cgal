#ifndef CGAL_POLYGON_SOUP_H
#define CGAL_POLYGON_SOUP_H

#include <vector>

template<typename Point_>
class Polygon_soup {
public:
    typedef std::size_t Point_handle;
    typedef std::size_t Polygon_handle;

    typedef Point_ Point;
    typedef std::vector<Point_handle> Polygon;

    typedef std::vector<Point> Point_range;
    typedef std::vector<Polygon> Polygon_range;

    typedef typename Point_range::iterator Point_iterator;
    typedef typename Polygon_range::iterator Polygon_iterator;

public:
    Polygon_soup() = default;

    Polygon_handle add_polygon(const Point& p0, const Point& p1, const Point& p2) {
        const Point_handle p0_idx = points.size();

        points.push_back(p0);
        points.push_back(p1);
        points.push_back(p2);

        Polygon t;
        t.push_back(p0_idx + 0);
        t.push_back(p0_idx + 1);
        t.push_back(p0_idx + 2);

        const Point_handle poly_idx = polygons.size();
        polygons.push_back(t);
        return poly_idx;
    }

    const Point& point(Point_handle p) const { return points[p]; }

    Point_range& point_range() { return points; }
    const Point_range& point_range() const { return points; }

    Point_iterator points_begin() { return points.begin(); }
    Point_iterator points_end() { return points.end(); }

    Polygon_range& polygon_range() { return polygons; }
    const Polygon_range& polygon_range() const { return polygons; }

    Polygon_iterator polygons_begin() { return polygons.begin(); }
    Polygon_iterator polygons_end() { return polygons.end(); }

private:
    Point_range points;
    Polygon_range polygons;
};

#endif // CGAL_POLYGON_SOUP_H
