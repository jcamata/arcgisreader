/* 
 * File:   InterpolateUniformGridData.h
 * Author: camata
 *
 * Created on October 19, 2015, 12:07 PM
 */

#ifndef INTERPOLATEUNIFORMGRIDDATA_H
#define	INTERPOLATEUNIFORMGRIDDATA_H

#include <vector>
#include <map>


class InterpolatedUniformGridData {
public:
    InterpolatedUniformGridData(const  std::vector< std::pair<double, double> > &interval_points, 
                                const  std::vector<int>   &n_subinterval,
                                const  std::vector< std::vector <double> > &data_values
            );
    double value(double x, double y);
    
private:
    std::vector< std::pair<double, double> > interval_points; 
    std::vector<int>                         n_subintervals;
    std::vector< std::vector <double> >      data_values; 
    double interpolate(std::vector<int> ix, double x, double y);
};


#endif	/* INTERPOLATEUNIFORMGRIDDATA_H */

