#ifndef _PFACET_H_
#define _PFACET_H_

// functor for priority queue
template <class Pfacet>
struct more_pfacet
{
    bool operator()(const Pfacet& f1,
                    const Pfacet& f2) const
    {
        return f1.priority() < f2.priority();
    }
};


template <class Facet_handle>
class CPfacet
{
private:
    Facet_handle m_facet;
    double m_priority;

public:

    // life cycle
    CPfacet(const Facet_handle& facet,
            const double& priority)
    {
        m_facet = facet;
        m_priority = priority;

    }

    // copy constructor
    CPfacet(const CPfacet& op)
    {
        m_facet = op.facet();
        m_priority = op.priority();
    }
    ~CPfacet() {}

public:

    // accessors
    Facet_handle& facet() {
        return m_facet;
    }
    const Facet_handle& facet() const {
        return m_facet;
    }

    double& priority() {
        return m_priority;
    }
    const double& priority() const {
        return m_priority;
    }


};

#endif // _PFACET_H_