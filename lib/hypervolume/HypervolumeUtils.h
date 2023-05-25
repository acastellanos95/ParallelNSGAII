//
// Created by andre on 8/9/22.
//

#ifndef NSGAII_LIB_HYPERVOLUME_HYPERVOLUMEUTILS_H_
#define NSGAII_LIB_HYPERVOLUME_HYPERVOLUMEUTILS_H_

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <bitset> // for calculation of trellis
#include <vector>

using namespace std;



/* function invoked by main */
void stream(double regLow[], double regUp[], const vector<double*>& cubs, int lev, double cov);
bool cmp(double* a, double* b);
/* function invoked by stream */
/*
inline bool cmp(double* a, double* b);
inline bool covers(const double* cub, const double regLow[]);
inline bool partCovers(const double* cub, const double regUp[]);
inline int containsBoundary(const double* cub, const double regLow[], const int split);
inline double getMeasure(const double regLow[], const double regUp[]);
inline int isPile(const double* cub, const double regLow[], const double regUp[]);
inline double computeTrellis(const double regLow[], const double regUp[], const double trellis[]);
inline double getMedian(vector<double>& bounds);
*/


/* global variables */
static int dataNumber;
static int dimension;
static double dSqrtDataNumber;
static double volume;

inline bool cmp(double* a, double* b) {
  return (a[dimension-1] < b[dimension-1]);
}


inline bool covers(const double* cub, const double regLow[]) {
  static int i;
  for (i=0; i<dimension-1; i++) {
    if (cub[i] > regLow[i]) {
      return false;
    }
  }
  return true;
}


inline bool partCovers(const double* cub, const double regUp[]) {
  static int i;
  for (i=0; i<dimension-1; i++) {
    if (cub[i] >= regUp[i]) {
      return false;
    }
  }
  return true;
}


inline int containsBoundary(const double* cub, const double regLow[], const int split) {
  // condition only checked for split>0
  if (regLow[split] >= cub[split]){
    // boundary in dimension split not contained in region, thus
    // boundary is no candidate for the splitting line
    return -1;
  }
  else {
    static int j;
    for (j=0; j<split; j++) { // check boundaries
      if (regLow[j] < cub[j]) {
        // boundary contained in region
        return 1;
      }
    }
  }
  // no boundary contained in region
  return 0;
}


inline double getMeasure(const double regLow[], const double regUp[]) {
  static double vol;
  static int i;
  vol = 1.0;
  for (i=0; i<dimension-1; i++) {
    vol *= (regUp[i] - regLow[i]);
  }
  return vol;
}


inline int isPile(const double* cub, const double regLow[], const double regUp[]) {
  static int pile;
  static int k;

  pile = dimension;
  // check all dimensions of the node
  for (k=0; k<dimension-1; k++) {
    // k-boundary of the node's region contained in the cuboid?
    if (cub[k] > regLow[k]) {
      if (pile != dimension) {
        // second dimension occured that is not completely covered
        // ==> cuboid is no pile
        return -1;
      }
      pile = k;
    }
  }
  // if pile == this.dimension then
  // cuboid completely covers region
  // case is not possible since covering cuboids have been removed before

  // region in only one dimenison not completly covered
  // ==> cuboid is a pile
  return pile;
}



inline double computeTrellis(const double regLow[], const double regUp[], const double trellis[]) {

  static int i,j;
  static double vol;
  static int numberSummands;
  static double summand;
  static bitset<16> bitvector;

  vol= 0.0;
  summand = 0.0;
  numberSummands = 0;

  // calculate number of summands
  static bitset<16> nSummands;
  for (i=0; i<dimension-1; i++) {
    nSummands[i] = 1;
  }
  numberSummands = nSummands.to_ulong();

  static double* valueTrellis = new double[dimension-1];
  static double* valueRegion = new double[dimension-1];
  for (i=0; i<dimension-1; i++) {
    valueTrellis[i] = trellis[i] - regUp[i];
  }
  for (i=0; i<dimension-1; i++) {
    valueRegion[i] = regUp[i] - regLow[i];
  }


  static double* dTemp = new double[numberSummands/2 + 1];

  // sum
  for (i=1; i<=numberSummands/2; i++) {

    // set bitvector length to fixed value 16
    // TODO Warning: dimension-1 <= 16 is assumed
    bitvector = (long)i;

    // construct summand
    // 0: take factor from region
    // 1: take factor from cuboid
    summand = 1.0;
    for (j=0; j<dimension-2; j++) {
      if (bitvector[j]) {
        summand *= valueTrellis[j];
      }
      else {
        summand *= valueRegion[j];
      }
    }
    summand *= valueRegion[dimension-2];

    // determine sign of summand
    vol -= summand;
    dTemp[i] =- summand;

    // add summand to sum
    // sign = (int) pow((double)-1, (double)counterOnes+1);
    // vol += (sign * summand);
  }


  bitvector = (long)i;
  summand = 1.0;
  for (j=0; j<dimension-1; j++) {
    if (bitvector[j]) {
      summand *= valueTrellis[j];
    }
    else {
      summand *= valueRegion[j];
    }
  }
  vol -= summand;

  for (i=1; i<=numberSummands/2; i++) {
    summand = dTemp[i];
    summand *= regUp[dimension-2] - trellis[dimension-2];
    summand /= valueRegion[dimension-2];
    vol -= summand;
  }

  //delete[] valueTrellis;
  //delete[] valueRegion;
  return vol;
}



// return median of the list of boundaries considered as a set
// TODO linear implementation
inline double getMedian(vector<double>& bounds) {
  // do not filter duplicates
  static unsigned int i;
  if (bounds.size()==1) {
    return bounds[0];
  }
  else if (bounds.size()==2) {
    return bounds[1];
  }
  vector<double>::iterator median;
  median = bounds.begin();
  for(i=0;i<=bounds.size()/2;i++){
    median++;
  }
  partial_sort(bounds.begin(),median,bounds.end());
  return bounds[bounds.size()/2];
}



// recursive calculation of hypervolume
inline void stream(double regionLow[], double regionUp[], const vector<double*>& points, int split, double cover) {

  //--- init --------------------------------------------------------------//

  static double coverOld;
  coverOld = cover;
  unsigned int coverIndex = 0;
  static int c;

  //--- cover -------------------------------------------------------------//

  // identify first covering cuboid
  double dMeasure = getMeasure(regionLow, regionUp);
  while (cover == coverOld && coverIndex < points.size()) {
    if ( covers(points[coverIndex], regionLow) ) {
      // new cover value
      cover = points[coverIndex][dimension-1];
      volume += dMeasure * (coverOld - cover);
    }
    else coverIndex++;
  }

  /* coverIndex shall be the index of the first point in points which
   * is ignored in the remaining process
   *
   * It may occur that that some points in front of coverIndex have the same
   * d-th coordinate as the point at coverIndex. This points must be discarded
   * and therefore the following for-loop checks for this points and reduces
   * coverIndex if necessary.
   */
  for (c=coverIndex; c>0; c--) {
    if (points[c-1][dimension-1] == cover) {
      coverIndex--;
    }
  }

  // abort if points is empty
  if (coverIndex == 0) {
    return;
  }
  // Note: in the remainder points is only considered to index coverIndex



  //--- allPiles  ---------------------------------------------------------//

  bool allPiles = true;
  unsigned int i;

  static int* piles = new int[coverIndex];
  for (i = 0; i<coverIndex; i++) {
    piles[i] = isPile(points[i], regionLow, regionUp);
    if (piles[i] == -1) {
      allPiles = false;
      //delete[] piles;
      break;
    }
  }

  /*
   * trellis[i] contains the values of the minimal i-coordinate of
   * the i-piles.
   * If there is no i-pile the default value is the upper bpund of the region.
   * The 1-dimensional KMP of the i-piles is: reg[1][i] - trellis[i]
   *
   */

  if (allPiles) { // sweep

    // initialize trellis with region's upper bound
    static double* trellis = new double[dimension-1];
    for (c=0; c<dimension-1; c++) {
      trellis[c] = regionUp[c];
    }

    double current = 0.0;
    double next = 0.0;
    i = 0;
    do { // while(next != coverNew)
      current = points[i][dimension-1];
      do { // while(next == current)
        if (points[i][piles[i]] < trellis[piles[i]]) {
          trellis[piles[i]] = points[i][piles[i]];
        }
        i++; // index of next point
        if (i < coverIndex) {
          next = points[i][dimension-1];
        }
        else {
          next = cover;
        }

      } while(next == current);
      volume += computeTrellis(regionLow, regionUp, trellis) * (next - current);
    } while(next != cover);
  }


    //--- split -------------------------------------------------------------//
    // inner node of partition tree
  else{
    double bound = -1.0;
    vector<double> boundaries;
    vector<double> noBoundaries;

    do {
      for (i=0; i<coverIndex; i++) {
        int contained = containsBoundary(points[i], regionLow, split);
        if (contained == 1) {
          boundaries.push_back(points[i][split]);
        } else if (contained == 0) {
          noBoundaries.push_back(points[i][split]);
        }
      }

      if (boundaries.size() >  0) {
        bound = getMedian(boundaries);
        //bound = getRandom(boundaries);
      }
      else if (noBoundaries.size() >  dSqrtDataNumber) {
        bound = getMedian(noBoundaries);
        //bound = getRandom(noBoundaries);
      }
      else {
        split++;
      }
    } while (bound == -1.0);

    double dLast;
    vector<double*> pointsChild;
    pointsChild.reserve(coverIndex);

    // left child
    // reduce maxPoint
    dLast = regionUp[split];
    regionUp[split] = bound;
    for (i=0; i<coverIndex; i++) {
      if (partCovers(points[i], regionUp)) {
        pointsChild.push_back(points[i]);
      }
    }
    if (!pointsChild.empty()) {
      stream(regionLow, regionUp, pointsChild, split, cover);
    }

    // right child
    // increase minPoint
    pointsChild.clear();
    regionUp[split] = dLast;
    dLast = regionLow[split];
    regionLow[split] = bound;
    for (i=0; i<coverIndex; i++) {
      if (partCovers(points[i], regionUp)) {
        pointsChild.push_back(points[i]);
      }
    }
    if (!pointsChild.empty()) {
      stream(regionLow, regionUp, pointsChild, split, cover);
    }
    regionLow[split] = dLast;

  }// end inner node

} // end stream

#endif //NSGAII_LIB_HYPERVOLUME_HYPERVOLUMEUTILS_H_
