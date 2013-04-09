template<typename PlaneType, typename PointType>
struct Is_between
{
public:
  bool operator()(const PointType& p) const
  {
    if(xcondition && p.x() < 0) 
      return false;
    PointType pp1 = plane1.projection(p);
    PointType pp2 = plane2.projection(p);
    double scal = (p.x() - pp1.x()) * (p.x() - pp2.x())
      + (p.y() - pp1.y()) * (p.y() - pp2.y())
      + (p.z() - pp1.z()) * (p.z() - pp2.z());
    return (scal <= 0.);
  }
public:
  Is_between()
    : plane1(), plane2(), xcondition(false){}
  Is_between(const PlaneType& p1, const PlaneType& p2, const bool xc)
    : plane1(p1), plane2(p2), xcondition(xc){}
  Is_between(const Is_between& ib)
    : plane1(ib.plane1), plane2(ib.plane2), xcondition(ib.xcondition){}
private:
  PlaneType plane1;
  PlaneType plane2;
  bool xcondition;
};
