/* 
 * File:   main.cpp
 * Author: camata
 *
 * Created on May 6, 2015, 3:48 PM
 */

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cmath>
using namespace std;

typedef struct {
    int nrows;
    int ncols;
    double xllcorner;
    double yllcorner;
    double cellsize;
    vector<double>  bathymetry;
} arcgis_data_t; 


void ArcGisReader(arcgis_data_t* grid, const char* fname)
{
    FILE *fin = fopen(fname,"r");
    if(!fin)
    {
        printf("Error: it was not possible to open %s\n", fname);
        exit(1);
    }
    double scale = 1.0E-05;
    
    fscanf(fin,"ncols       %d\n", &grid->ncols);
    fscanf(fin,"nrows       %d\n", &grid->nrows);
    fscanf(fin,"xllcorner   %lf\n", &grid->xllcorner);
    fscanf(fin,"yllcorner   %lf\n", &grid->yllcorner);
    fscanf(fin,"cellsize    %lf\n", &grid->cellsize);
    grid->bathymetry.resize(grid->ncols*grid->nrows);
    for(int i=0; i < grid->nrows; ++i)
    {
        for(int j = 0; j < grid->ncols; ++j)
        {
            int tmp;
            fscanf(fin,"%d", &tmp);
            grid->bathymetry[i*grid->ncols+j] = (double) scale*tmp;
        }
    }
    
    fclose(fin);
}


void ArcGisWriteSTL(arcgis_data_t* grid, const char* fname)
{
    FILE *fout = fopen(fname, "w");
     if(!fout)
    {
        printf("Error: it was not possible to open %s\n", fname);
        exit(1);
    }
    
    fprintf(fout,"solid bathymetry\n");
    double x[2];
    double y[2];
    double z[2][2];
    double V[3][3];
    double scale = 1;
    
    //cell loop
    for(int i=0; i < (grid->nrows-1); ++i)
    {
        y[0] = scale*(grid->yllcorner+(i  )*grid->cellsize);
        y[1] = scale*(grid->yllcorner+(i+1)*grid->cellsize);
        
        for(int j = 0; j < (grid->ncols-1); ++j)
        {
            x[0] = scale*(grid->xllcorner+(j  )*grid->cellsize);
            x[1] = scale*(grid->xllcorner+(j+1)*grid->cellsize);
                   
            z[0][0] = grid->bathymetry[(i  )*grid->ncols+(j  )];
            z[1][0] = grid->bathymetry[(i  )*grid->ncols+(j+1)];
            z[0][1] = grid->bathymetry[(i+1)*grid->ncols+(j  )];
            z[1][1] = grid->bathymetry[(i+1)*grid->ncols+(j+1)];       
                   
            //facet 1:
            //  z[0][1]
            //     +----------+ z[1][1]
            //     |          |
            //     |          |
            //     |          |
            //     +----------+
            //     z[0][0]    z[1][0]
            {
                V[0][0] = x[0];
                V[0][1] = y[0];
                V[0][2] = z[0][0];
                
                V[1][0] = x[1];
                V[1][1] = y[0];
                V[1][2] = z[1][0];
                
                V[2][0] = x[1];
                V[2][1] = y[1];
                V[2][2] = z[1][1];
                
                double vx = V[1][0] - V[0][0];
                double vy = V[1][1] - V[0][1];
                double vz = V[1][2] - V[0][2];

                double wx = V[2][0] - V[0][0];
                double wy = V[2][1] - V[0][1];
                double wz = V[2][2] - V[0][2];

                double nx = vy*wz-vz*wy;
                double ny = vz*wx-vx*wz;
                double nz = vx*wy-vy*wx;
                double norm = sqrt(nx*nx+ny*ny+nz*nz);
                nx /= norm;
                ny /= norm;
                nz /= norm;
                
                fprintf(fout,"facet normal %8.8f %8.8f %8.8f\n", nx, ny, nz);
                fprintf(fout,"  outer loop\n");
                fprintf(fout,"     vertex %8.8f %8.8f %8.8f\n", V[0][0], V[0][1], V[0][2]);
                fprintf(fout,"     vertex %8.8f %8.8f %8.8f\n", V[1][0], V[1][1], V[1][2]);
                fprintf(fout,"     vertex %8.8f %8.8f %8.8f\n", V[2][0], V[2][1], V[2][2]);
                fprintf(fout,"  endloop\n");
                fprintf(fout,"endfacet\n");
            }
            
            {
                V[0][0] = x[0];
                V[0][1] = y[0];
                V[0][2] = z[0][0];
                
                V[1][0] = x[1];
                V[1][1] = y[1];
                V[1][2] = z[1][1];
                
                V[2][0] = x[0];
                V[2][1] = y[1];
                V[2][2] = z[0][1];
                
                double vx = V[1][0] - V[0][0];
                double vy = V[1][1] - V[0][1];
                double vz = V[1][2] - V[0][2];

                double wx = V[2][0] - V[0][0];
                double wy = V[2][1] - V[0][1];
                double wz = V[2][2] - V[0][2];

                double nx = vy*wz-vz*wy;
                double ny = vz*wx-vx*wz;
                double nz = vx*wy-vy*wx;
                
                double norm = sqrt(nx*nx+ny*ny+nz*nz);
                nx /= norm;
                ny /= norm;
                nz /= norm;
                
                fprintf(fout,"facet normal %8.8f %8.8f %8.8f\n", nx, ny, nz);
                fprintf(fout,"  outer loop\n");
                fprintf(fout,"     vertex %8.8f %8.8f %8.8f\n", V[0][0], V[0][1], V[0][2]);
                fprintf(fout,"     vertex %8.8f %8.8f %8.8f\n", V[1][0], V[1][1], V[1][2]);
                fprintf(fout,"     vertex %8.8f %8.8f %8.8f\n", V[2][0], V[2][1], V[2][2]);
                fprintf(fout,"  endloop\n");
                fprintf(fout,"endfacet\n");
            }
        }
    }
    fprintf(fout,"endsolid bathymetry\n");
    fclose(fout);
    
}


void ArcGisWriteGmshGEO(arcgis_data_t* grid, const char* fname)
{
    int np = 0;
    FILE *fout = fopen(fname, "w");
     if(!fout)
    {
        printf("Error: it was not possible to open %s\n", fname);
        exit(1);
    }
    
     fprintf(fout,"lc = %f;\n", grid->cellsize);
    
    //cell loop
    for(int i=0; i < grid->nrows; ++i)
    {
        double y = grid->yllcorner+(i)*grid->cellsize;
       
        for(int j = 0; j < grid->ncols; ++j)
        {
            np++;
            
            double x = grid->xllcorner+(j)*grid->cellsize;            
            double z = grid->bathymetry[i*grid->ncols+j];
                          
            fprintf(fout,"Point(%d) = {%8.8f, %8.8f, %8.8f, lc};\n", np, x, y, z);

        }
    }

    fclose(fout);
    
}



/*
 * 
 */
int main(int argc, char** argv) {
  
    arcgis_data_t* grid = new arcgis_data_t();
    
    ArcGisReader(grid, "etopo1_bedrock.asc");
    ArcGisWriteSTL(grid, "etopo1.stl");
    //ArcGisWriteGmshGEO(grid,"etopo.geo");
    
            
    delete grid;
    return 0;
}

