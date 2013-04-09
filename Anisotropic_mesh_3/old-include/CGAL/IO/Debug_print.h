#ifndef DEBUG_PRINT_H
#define DEBUG_PRINT_H

#ifndef NDEBUG

#include <iostream>

template<class Point_3>
void print_point(const Point_3 &p){
	std::cout<<" "<<p.x()<<" "<<p.y()<<" "<<p.z()<<"  0 "<<std::endl;
}

#endif //NDEBUG

#endif //DEBUG_PRINT_H
