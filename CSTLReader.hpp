/*******************************************************************************
 * Copyright (c) 2017  IBM Corporation, Carnegie Mellon University and others
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *******************************************************************************/

#ifndef CSTLReader_hpp
#define CSTLReader_hpp

#include <stdio.h>
#include <fstream>
#include "opencv2/opencv.hpp"

#endif /* CSTLReader_hpp */

// Author: Marco Visentini-Scarzanella
// www.commsp.ee.ic.ac.uk/~marcovs


using namespace std;
using namespace cv;

//TODO: use indexes for vertices and have an array of faces
class STLReader
{
public:
    STLReader(void);
    ~STLReader(void);
    

    //input
    void addPointsToCloud(const vector<Point3d> &pcl);
    void addPointsToCloud(const vector<Point3d> &pcl, const vector<Vec3i> &colour);
    void addPointsToCloud(const vector<Matx31d> &pcl);
    void addPointsToCloud(const vector<Matx31d> &pcl, const vector<Vec3i> &colour);
    
    //geometric transforms
    void applyRegistration(const Matx33d &R, const Vec3d &t);
    void normaliseVolume();
    void scaleVolume(double scaleF);
    void centerVolume();
    
    //file IO
    void writePLYPointCloud(const string fileName);
    
private:
    void updateCentroid();
    
    vector<Point3d> _volume;
    vector<Vec3d> _normal;
    vector<Vec3d> _colour;
    
    int _nPoints;
    double _maxNorm;
    Point3d _centroid;
};
