
#include "InterpolateUniformGridData.h"


InterpolatedUniformGridData::InterpolatedUniformGridData(const  std::vector< std::pair<double, double> > &interval_points, 
                                const  std::vector<int>   &n_subintervals,
                                const  std::vector< std::vector <double> > &data_values)
        : interval_points(interval_points),
          n_subintervals(n_subintervals),
          data_values(data_values)
{
    
}

double InterpolatedUniformGridData::interpolate(std::vector<int> ix, double x, double y)
{
     return (((1-x)*data_values[ix[0]][ix[1]]
               +
               x*data_values[ix[0]+1][ix[1]])*(1-y)
              +
              ((1-x)*data_values[ix[0]][ix[1]+1]
               +
               x*data_values[ix[0]+1][ix[1]+1])*y);
}


double InterpolatedUniformGridData::value(double x, double y)
{
    // find out where this data point lies, relative to the given
    // subdivision points
    
    std::vector<int> ix(2);
    std::vector<double> p(2);
    p[0] = x;
    p[1] = y;
    
    for (unsigned int d=0; d<2; ++d)
      {
        const double delta_x = ((interval_points[d].second - interval_points[d].first) /
                                n_subintervals[d]);
        if (p[d] <= interval_points[d].first)
          ix[d] = 0;
        else if (p[d] >= interval_points[d].second-delta_x)
          ix[d] = n_subintervals[d]-1;
        else
          ix[d] = (unsigned int)((p[d]-interval_points[d].first) / delta_x);
      }

    // now compute the relative point within the interval/rectangle/box
    // defined by the point coordinates found above. truncate below and
    // above to accommodate points that may lie outside the range
    std::vector<double> p_unit(2);
    for (unsigned int d=0; d<2; ++d)
      {
        const double delta_x = ((interval_points[d].second - interval_points[d].first) /
                                n_subintervals[d]);
        p_unit[d] = std::max(std::min((p[d]-interval_points[d].first-ix[d]*delta_x)/delta_x,1.), 0.);
      }

    return interpolate (ix, p_unit[0], p_unit[1]);
    
}