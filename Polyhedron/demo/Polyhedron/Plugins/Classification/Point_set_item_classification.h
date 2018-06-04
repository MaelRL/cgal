#ifndef POINT_SET_ITEM_CLASSIFICATION_H
#define POINT_SET_ITEM_CLASSIFICATION_H

//#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#define CGAL_CLASSIFICATION_VERBOSE
#define CGAL_CLASSTRAINING_VERBOSE

#include <CGAL/Three/Scene_item.h>

#include "Scene_points_with_normal_item.h"
#include "Item_classification_base.h"
#include "Kernel_type.h"
#include "Point_set_3.h"

#include <CGAL/Classification.h>

#include <iostream>

#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif



// This class represents a point set in the OpenGL scene
class Point_set_item_classification : public Item_classification_base
{
 public:
  typedef Kernel::Point_3 Point_3;
  typedef Kernel::Vector_3 Vector_3;
  typedef CGAL::Classification::RGB_Color Color;
  
  typedef Point_set::Point_map Point_map;
  typedef Point_set::Vector_map Vector_map;

  typedef CGAL::Classification::Point_set_feature_generator<Kernel, Point_set, Point_map>               Generator;
  
  struct Cluster
  {
    std::vector<Point_set::Index> inliers;
    int training;
    int label;
    Cluster() : training (-1), label (-1) { }

    std::size_t size() const { return inliers.size(); }
    const Point_set::Index& operator[] (std::size_t i) const { return inliers[i]; }
  };
  
  struct Cluster_neighborhood
  {
    Point_set* point_set;
    Point_set::Property_map<int> cluster_id;
    std::vector<Cluster>* clusters;
    
    Cluster_neighborhood (Point_set* point_set,
                          std::vector<Cluster>& clusters)
      : point_set (point_set)
      , clusters (&clusters)
    {
      cluster_id = point_set->property_map<int>("shape").first;
    }
    
    template <typename OutputIterator>
    OutputIterator operator() (const Point_set::Index& idx,
                               OutputIterator output) const
    {
      int c = cluster_id[idx];
      if (c == -1)
        *(output ++) = idx;
      else
      {
        std::copy ((*clusters)[c].inliers.begin(),
                   (*clusters)[c].inliers.end(),
                   output);
      }
      return output;
    }
  };
  
 public:
  
  Point_set_item_classification(Scene_points_with_normal_item* points);
  ~Point_set_item_classification();

  CGAL::Three::Scene_item* item() { return m_points; }
  void erase_item() { m_points = NULL; }

  CGAL::Bbox_3 bbox()
  {
    if (m_points->point_set()->nb_selected_points() == 0)
      return m_points->bbox();

    CGAL::Bbox_3 bb = CGAL::bbox_3 (boost::make_transform_iterator
                                    (m_points->point_set()->first_selected(),
                                     CGAL::Property_map_to_unary_function<Point_set::Point_map>
                                     (m_points->point_set()->point_map())),
                                    boost::make_transform_iterator
                                    (m_points->point_set()->end(),
                                     CGAL::Property_map_to_unary_function<Point_set::Point_map>
                                     (m_points->point_set()->point_map())));

    double xcenter = (bb.xmax() + bb.xmin()) / 2.;
    double ycenter = (bb.ymax() + bb.ymin()) / 2.;
    double zcenter = (bb.zmax() + bb.zmin()) / 2.;

    double dx = bb.xmax() - bb.xmin();
    double dy = bb.ymax() - bb.ymin();
    double dz = bb.zmax() - bb.zmin();
    
    dx *= 10.;
    dy *= 10.;
    dz *= 10.;

    return CGAL::Bbox_3 (xcenter - dx, ycenter - dy, zcenter - dz,
                         xcenter + dx, ycenter + dy, zcenter + dz);
  }

  void compute_features (std::size_t nb_scales);
  void add_remaining_point_set_properties_as_features();
  
  void select_random_region();

  template <typename Type>
  bool try_adding_simple_feature (const std::string& name)
  {
    typedef typename Point_set::template Property_map<Type> Pmap;
    bool okay = false;
    Pmap pmap;
    boost::tie (pmap, okay) = m_points->point_set()->template property_map<Type>(name.c_str());
    if (okay)
      m_features.template add<CGAL::Classification::Feature::Simple_feature <Point_set, Pmap> >
        (*(m_points->point_set()), pmap, name.c_str());

    return okay;
  }
  
  void add_selection_to_training_set (std::size_t label)
  {
    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
    {
      m_training[*it] = int(label);
      m_classif[*it] = int(label);
    }

    m_points->resetSelection();
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void reset_training_set(std::size_t label)
  {
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
      if (m_training[*it] == int(label))
        m_training[*it] = -1;
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void reset_training_set_of_selection()
  {
    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
    {
      m_training[*it] = -1;
      m_classif[*it] = -1;
    }
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void reset_training_sets()
  {
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
      m_training[*it] = -1;
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void validate_selection ()
  {
    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
      m_training[*it] = m_classif[*it];

    m_points->resetSelection();
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void train(int classifier, unsigned int nb_trials,
             std::size_t num_trees, std::size_t max_depth);
  bool run (int method, int classifier, std::size_t subdivisions, double smoothing);

  void update_color () { change_color (m_index_color); }
  void change_color (int index);
  CGAL::Three::Scene_item* generate_one_item (const char* name,
                                              int label) const
  {
    Scene_points_with_normal_item* points_item
      = new Scene_points_with_normal_item;
    
    points_item->setName (QString("%1 (%2)").arg(name).arg(m_labels[label]->name().c_str()));
    points_item->setColor (m_label_colors[label]);
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
    {
      int c = m_classif[*it];
      if (c == label)
        points_item->point_set()->insert (m_points->point_set()->point(*it));
    }
    return points_item;
  }
  void generate_one_item_per_label(std::vector<CGAL::Three::Scene_item*>& items,
                                   const char* name) const
  {
    std::vector<Scene_points_with_normal_item*> points_item
      (m_labels.size(), NULL);
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      {
        points_item[i] = new Scene_points_with_normal_item;
        points_item[i]->setName (QString("%1 (%2)").arg(name).arg(m_labels[i]->name().c_str()));
        points_item[i]->setColor (m_label_colors[i]);
        items.push_back (points_item[i]);
      }

    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
      {
        int c = m_classif[*it];
        if (c != -1)
          points_item[c]->point_set()->insert (m_points->point_set()->point(*it));
      }
  }
  
  bool write_output(std::ostream& out);

  QColor add_new_label (const char* name)
  {
    QColor out = Item_classification_base::add_new_label (name);
    update_comments_of_point_set_item();
    return out;
  }

  void remove_label (std::size_t position)
  {
    Item_classification_base::remove_label (position);

    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
    {
      if (m_training[*it] == int(position))
        m_training[*it] = -1;
      else if (m_training[*it] > int(position))
        m_training[*it] --;
      
      if (m_classif[*it] == int(position))
        m_classif[*it] = -1;
      else if (m_classif[*it] > int(position))
        m_classif[*it] --;
    }
    update_comments_of_point_set_item();
  }
  
  int real_index_color() const;
  void reset_indices();
  void backup_existing_colors_and_add_new();
  void reset_colors();

 private:

  void update_comments_of_point_set_item()
  {
    std::string& comments = m_points->comments();
    
    // Remove previously registered labels from comments
    std::string new_comment;
      
    std::istringstream stream (comments);
    std::string line;
    while (getline(stream, line))
    {
      std::stringstream iss (line);
      std::string tag;
      if (iss >> tag && tag == "label")
        continue;
      new_comment += line + "\n";
    }
    comments = new_comment;

    comments += "label -1 unclassified\n";
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
    {
      std::ostringstream oss;
      oss << "label " << i << " " << m_labels[i]->name() << std::endl;
      comments += oss.str();
    }
  }
  
  template <typename Classifier>
  bool run (int method, const Classifier& classifier,
            std::size_t subdivisions, double smoothing)
  {
    std::vector<int> indices (m_points->point_set()->size(), -1);

    if (method == 0)
      CGAL::Classification::classify<Concurrency_tag> (*(m_points->point_set()),
                                                       m_labels, classifier,
                                                       indices);
    else if (method == 1)
    {
      if (m_clusters.empty()) // Use real local smoothing
        CGAL::Classification::classify_with_local_smoothing<Concurrency_tag>
          (*(m_points->point_set()), m_points->point_set()->point_map(), m_labels, classifier,
           m_generator->neighborhood().sphere_neighbor_query(m_generator->radius_neighbors()),
           indices);
      else // Smooth on clusters
      {
        std::cerr << "Smoothing on clusters" << std::endl;
        CGAL::Classification::classify_with_local_smoothing<Concurrency_tag>
          (*(m_points->point_set()),
           CGAL::Identity_property_map<Point_set::Index>(),
           m_labels, classifier,
           Cluster_neighborhood(m_points->point_set(),
                                m_clusters),
           indices);
      }
    }
    else if (method == 2)
      CGAL::Classification::classify_with_graphcut<Concurrency_tag>
        (*(m_points->point_set()), m_points->point_set()->point_map(),
         m_labels, classifier,
         m_generator->neighborhood().k_neighbor_query(12),
         smoothing, subdivisions, indices);

    std::vector<int> ground_truth(m_points->point_set()->size(), -1);
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->first_selected(); ++ it)
      {
        m_classif[*it] = indices[*it];
        ground_truth[*it] = m_training[*it];
      }
  
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);

    std::cerr << "Precision, recall, F1 scores and IoU:" << std::endl;
    
    CGAL::Classification::Evaluation eval (m_labels, ground_truth, indices);
  
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      {
        std::cerr << " * " << m_labels[i]->name() << ": "
                  << eval.precision(m_labels[i]) << " ; "
                  << eval.recall(m_labels[i]) << " ; "
                  << eval.f1_score(m_labels[i]) << " ; "
                  << eval.intersection_over_union(m_labels[i]) << std::endl;
      }

    std::cerr << "Accuracy = " << eval.accuracy() << std::endl
              << "Mean F1 score = " << eval.mean_f1_score() << std::endl
              << "Mean IoU = " << eval.mean_intersection_over_union() << std::endl;


    return true;
  }

  Scene_points_with_normal_item* m_points;

  std::vector<Cluster> m_clusters;

  std::vector<std::vector<float> > m_label_probabilities;

  Point_set::Property_map<unsigned char> m_red;
  Point_set::Property_map<unsigned char> m_green;
  Point_set::Property_map<unsigned char> m_blue;
  Point_set::Property_map<Color> m_color;
  Point_set::Property_map<int> m_training;
  Point_set::Property_map<int> m_classif;

  Generator* m_generator;
  
  int m_index_color;
  
}; // end class Point_set_item_classification




#endif // POINT_SET_ITEM_CLASSIFICATION_H
