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

#include "CSTLReader.hpp"



STLReader::STLReader(void)
{
    _centroid = {0,0,0};
    _nPoints = 0;
    _maxNorm = 0;
}

STLReader::~STLReader(void)
{
}

void STLReader::addPointsToCloud(const vector<Point3d> &pcl) {
    _volume.reserve(_volume.size() + pcl.size());
    _colour.reserve(_colour.size() + pcl.size());
    
    for (int i = 0; i < pcl.size(); i++) {
        _volume.push_back(pcl[i]);
        _colour.push_back(Vec3d(128,128,128));
        _maxNorm = (_maxNorm > norm(pcl[i])) ? _maxNorm : norm(pcl[i]);
    }
    
    _nPoints += pcl.size();
    updateCentroid();
}

void STLReader::addPointsToCloud(const vector<Matx31d> &pcl) {
    _volume.reserve(_volume.size() + pcl.size());
    _colour.reserve(_colour.size() + pcl.size());
    
    for (int i = 0; i < pcl.size(); i++) {
        _volume.push_back(Point3d(pcl[i].val[0],pcl[i].val[1],pcl[i].val[2]));
        _colour.push_back(Vec3d(128,128,128));
        _maxNorm = (_maxNorm > norm(pcl[i])) ? _maxNorm : norm(pcl[i]);
    }
    
    _nPoints += pcl.size();
    updateCentroid();
}

void STLReader::addPointsToCloud(const vector<Point3d> &pcl, const vector<Vec3i> &colour) {
    _volume.reserve(_volume.size() + pcl.size());
    _colour.reserve(_colour.size() + pcl.size());
    
    for (int i = 0; i < pcl.size(); i++) {
        _volume.push_back(pcl[i]);
        _colour.push_back(colour[i]);
        _maxNorm = (_maxNorm > norm(pcl[i])) ? _maxNorm : norm(pcl[i]);
    }
    
    _nPoints += pcl.size();
    updateCentroid();
}

void STLReader::addPointsToCloud(const vector<Matx31d> &pcl, const vector<Vec3i> &colour) {
    _volume.reserve(_volume.size() + pcl.size());
    _colour.reserve(_colour.size() + pcl.size());
    
    for (int i = 0; i < pcl.size(); i++) {
        _volume.push_back(Point3d(pcl[i].val[0],pcl[i].val[1],pcl[i].val[2]));
        _colour.push_back(colour[i]);
        _maxNorm = (_maxNorm > norm(pcl[i])) ? _maxNorm : norm(pcl[i]);
    }
    
    _nPoints += pcl.size();
    updateCentroid();
}

//scales the volume to be beween -0.5 and 0.5
void STLReader::normaliseVolume() {
    for (int i = 0; i < _volume.size(); i++)
        _volume[i] = (_volume[i] - _centroid)*(1.0/_maxNorm);
    
    _maxNorm = 1;
}

void STLReader::updateCentroid() {
    //calculate centroid based on median to avoid outliers
    _centroid = {0,0,0};
    vector<double> x, y, z;
    x.reserve(_volume.size());
    y.reserve(_volume.size());
    z.reserve(_volume.size());
    
    for (int i = 0; i < _nPoints; i++) {
        x.push_back(_volume[i].x);
        y.push_back(_volume[i].y);
        z.push_back(_volume[i].z);
    }
    
    sort(x.begin(),x.end());
    sort(y.begin(),y.end());
    sort(z.begin(),z.end());
    int idx = round(_nPoints/2.0);
    _centroid = {x[idx],y[idx],z[idx]};
}

void STLReader::centerVolume() {
    updateCentroid();
    
    for (int i = 0; i < _nPoints; i++) {
        _volume[i] -= _centroid;
    }
    updateCentroid();
}

void STLReader::writePLYPointCloud(const string fileName) {
    //output point cloud, no connectivity or normal information
    ofstream fout;
    fout.open(fileName);
    //header
    fout << "ply" << endl;
    fout << "format ascii 1.0" << endl;
    fout << "comment VCGLIB generated" << endl;
    fout << "element vertex " << _nPoints << endl;
    fout << "property float x" << endl;
    fout << "property float y" << endl;
    fout << "property float z" << endl;
    fout << "property uchar red" << endl;
    fout << "property uchar green" << endl;
    fout << "property uchar blue" << endl;
    fout << "property uchar alpha" << endl;
    fout << "element face 0" << endl;
    fout << "property list uchar int vertex_indices" << endl;
    fout << "end_header" << endl;
    
    for (int i = 0; i < _nPoints; i++) {
        fout << _volume[i].x << " " << _volume[i].y << " " << _volume[i].z << " " << _colour[i][0] << " " << _colour[i][1] << " " << _colour[i][2] << " 255" << endl;
    }

    fout.close();
}

void STLReader::applyRegistration(const Matx33d &R, const Vec3d &t) {
    //compose transformation matrix
    Matx44d P(R(0,0), R(0,1), R(0,2), t(0), R(1,0), R(1,1), R(1,2), t(1), R(2,0), R(2,1), R(2,2), t(2), 0, 0, 0, 1);
    //convert points to homogeneous coordinates
    Mat pts;
    convertPointsToHomogeneous(_volume, pts);
    convertPointsFromHomogeneous(P*pts, _volume);
    
    //apply to centroid
    _centroid = R*_centroid + Point3d(t);
    
    //if normals are present, apply rotation to normals
    if (_normal.size() != 0) {
        P = Matx44d(R(0,0), R(0,1), R(0,2), 0, R(1,0), R(1,1), R(1,2), 0, R(2,0), R(2,1), R(2,2), 0, 0, 0, 0, 1);
        Mat n;
        convertPointsFromHomogeneous(_normal,n);
        convertPointsToHomogeneous(P*n, _normal);
    }
}

void STLReader::scaleVolume(double scaleF) {
    for (int i = 0; i < _nPoints; i++) {
        _volume[i] *= scaleF;
    }
    _centroid *= scaleF;
    _maxNorm *= scaleF;
}

