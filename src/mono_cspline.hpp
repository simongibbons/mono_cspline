#ifndef _MONO_CSPLINE_HPP
#define _MONO_CSPLINE_HPP

#include <vector>
#include <math>

class MonoCspline {

    double xmin, xmax;

    unsigned int npts;

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> m;

public:
    MonoCspline(std::vector<double> x, std::vector<double> y);
    MonoCspline( double *x, double *y, int npts );

    double operator() ( const double xpoint ) const;

private:
    inline double basis_00( const double t ) const {
        return (1 + 2*t) * pow( 1 - t, 2.0 );
    }

    inline double basis_10( const double t ) const {
        return t * pow( 1-t, 2.0 );
    }

    inline double basis_01( const double t ) const {
        return pow(t, 2.0) * (3 - 2*t);
    }

    inline double basis_11( const double t ) const {
        return pow(t, 2.0) * (t - 1.0);
    }

    void make_tangents();

    int find_index( double x );

};


#endif
