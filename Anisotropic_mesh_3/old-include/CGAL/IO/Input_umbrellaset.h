#ifndef INPUT_UMBRELLASET_H
#define INPUT_UMBRELLASET_H

/** Reads a Umbrella_set type structure and
 * outputs a mesh file to the assigned out stream
 */

#include <iostream>
#include <map>
#include <iterator>
#include <vector>
#include <ctime>

#include <CGAL/Conflict_surface_triangles.h>///TODO to nest it in the umbrella_set

template <typename Umbrella_set>
bool init_umbrella_set_from_file(Umbrella_set &umbs, std::istream &input_file){
	typedef typename Umbrella_set::Point		Point;

	std::vector<Point>	all_points;
	std::vector<bool>	all_on_constraint;
	unsigned points_size;
	input_file>>points_size;
	all_points.reserve(points_size);
	all_on_constraint.reserve(points_size);
	
#ifdef DEBUG_READING_FROM_FILE
	std::cout<<points_size<<" points expected from file."<<std::endl;
#endif //DEBUG_READING_FROM_FILE

	unsigned point_index = 0;
	for(;point_index<points_size && ( (!input_file.eof()) &&(!input_file.bad()));){
		std::string char_cache;
		getline(input_file,char_cache);
		if(char_cache.size()>4 && char_cache[0]!='#'){
			double x_axis,y_axis,z_axis;
			int on_constraint;
			std::istringstream char_stream(char_cache);
			char_stream>>x_axis>>y_axis>>z_axis;
			if(char_stream.fail()){
				std::cerr<<"Warning: error in reading data for point "<<point_index<<std::endl;
			}else{
				all_points.push_back(Point(x_axis,y_axis,z_axis));
				point_index ++;
#ifdef DEBUG_READING_FROM_FILE
				std::cout<<"("<<x_axis<<", "<<y_axis<<", "<<z_axis<<" )"<<std::endl;
#endif //DEBUG_READING_FROM_FILE
				if(!char_stream.eof()){
					char_stream >> on_constraint;
					if(!char_stream.fail()){
						all_on_constraint.push_back((on_constraint == 0));
					}
				}
			}
		}
	}
	if(point_index < points_size){
		std::cerr<<"Error in reading file."<<std::endl;
		return false;
	}

	unsigned umb_index =0;
	unsigned percentage = 0;
	while(umb_index<all_points.size()){
		std::string char_cache;
		getline(input_file,char_cache);
		if(input_file.eof() || input_file.fail()){
			break;
		}
		
		unsigned neighbor_index;
		std::vector<unsigned> neighbor_index_list;
		neighbor_index_list.reserve(32);

		if(char_cache.size()>4 && char_cache[0]!='#'){
			unsigned size ;
			std::istringstream char_stream(char_cache);
			char_stream>>size;
#ifdef DEBUG_READING_FROM_FILE
			std::cout<<"size of neighbor of umb "<<umb_index<<":"<<size<<std::endl;
#endif //DEBUG_READING_FROM_FILE
			neighbor_index_list.clear();
			for(unsigned j=0;j<size;j++){
				if(char_stream.eof()){
					std::cerr<<"Warning: Missing neighbor data for umbrella "<<umb_index<<std::endl;
					return false;
				}
				char_stream>>neighbor_index;
				if(char_stream.fail()){
					std::cerr<<"Warning: Error in reading neighbor data for umbrella "<<umb_index<<std::endl;
					return false;
				}
				neighbor_index_list.push_back(neighbor_index);
#ifdef DEBUG_READING_FROM_FILE
				std::cout<<neighbor_index<<"  "<<std::endl;
#endif //DEBUG_READING_FROM_FILE
			}
			if(!char_stream.fail()){
				if(umb_index==0){
#ifdef DEBUG_READING_FROM_FILE
					std::cout<<"Init all points for umbrella with "<<all_points.size()<<" points"<<std::endl;
#endif
					if(all_on_constraint.size()==all_points.size()){
						umbs.set_all_points(all_points,true,&all_on_constraint);
					}else{
						umbs.set_all_points(all_points,true);	//true for flag of empty umbs
					}
				}
				umbs.init_umbrella(umb_index,neighbor_index_list);
				umb_index++;
				if(  (umb_index * 100 / all_points.size()) == percentage + 5){
					percentage += 5;
					std::cout<<"\nInitialization finished "<<percentage<<"\%";
				}
				if(umb_index%50==0){
					std::cout<<".";
					std::cout.flush();
				}
			}
		}//if line is valid
	}

	if(umb_index == 0){
		std::cout<<"Default initialization"<<std::endl;
		if(all_on_constraint.size()==all_points.size()){
			umbs.set_all_points(all_points,true,&all_on_constraint);
		}else{
			umbs.set_all_points(all_points);
		}
	}else if(umb_index < all_points.size()){
		std::cerr <<"Warning: End of data meet before completing umbrella set initialization."<<std::endl;
		return false;
	}
	return true;
}


#endif //INPUT_UMBRELLASET_H
