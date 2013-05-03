//
//  kdTreeStruct.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 17/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

// Most of the class is nanoflann kdtree example, with some modifications

#ifndef FFT_Chain_kdTreeStruct_h
#define FFT_Chain_kdTreeStruct_h
#include "SetupMacros.h"
#include <math.h>
#define DIMEN FNMELBANDS * MFCCGROUPING
struct PointCloud
{
	struct Point{
		double  array[DIMEN];
        int original_index;
	};
    
	std::vector<Point>  pts;
	inline size_t kdtree_get_point_count() const { return pts.size(); }
	inline double kdtree_distance(const double *p1, const size_t idx_p2,size_t size) const{
        double sum = 0;
       // printf("changed");
        double runningSum = 0;
        for(int i= 0; i <DIMEN; i++){
            sum = ((Point*)p1)->array[i] - pts[idx_p2].array[i];
            sum *= sum;
            runningSum += sum;
        }
        //runningSum /= (double)DIMEN;
        
		return sqrt(runningSum);
	}
    inline void add_database(double **mfcc_vec,const int num_coef,const long num_frames){
        for(int i = 0 ; i < num_frames; i++){
            Point new_point;
            for(int j = 0 ; j < num_coef; j++){
                new_point.array[j] = mfcc_vec[i][j];
            }
            new_point.original_index = i;
            pts.push_back(new_point);
        }
    }
	inline double kdtree_get_pt(const size_t idx, int dim) const{
        return pts[idx].array[dim];
        //		if (dim==0) return pts[idx].x;
        //		else if (dim==1) return pts[idx].y;
        //		else return pts[idx].z;
	}
    inline double kdtree_get_orig_index(const size_t idx) const{
        return pts[idx].original_index;
    }
	template <class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const { return false; }
    
};

#endif
