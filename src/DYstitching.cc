#include "Tools/Tools/interface/DYstitching.h"

// Constructors

DYstitching::DYstitching (int year) {
  if (year == 2016) {
    stitchWeights = stitchWeights_2016;
    
  } else if (year == 2017) {
    stitchWeights = stitchWeights_2017;

  } else {
    stitchWeights = stitchWeights_2018;
  }

};

// Destructor
DYstitching::~DYstitching() {}

float DYstitching::get_stitching_weight(int LHE_Nb, int LHE_Njets, float LHE_HT) {
  int njets = LHE_Njets;
  int nb = LHE_Nb;
  // these protections should be useless
  if (njets < 0) njets = 0;
  if (njets > 4) njets = 4;
  if (nb < 0)    nb = 0;
  if (nb > 4)    nb = 4;

  float ht = LHE_HT;
  int nht = 0;
  if      (ht  < 0                ) nht = 0;
  else if (ht >= 0    && ht < 70  ) nht = 0;
  else if (ht >= 70   && ht < 100 ) nht = 1;
  else if (ht >= 100  && ht < 200 ) nht = 2;
  else if (ht >= 200  && ht < 400 ) nht = 3;
  else if (ht >= 400  && ht < 600 ) nht = 4;
  else if (ht >= 600  && ht < 800 ) nht = 5;
  else if (ht >= 800  && ht < 1200) nht = 6;
  else if (ht >= 1200 && ht < 2500) nht = 7;
  else  /* ht >= 2500 */            nht = 8;

  return stitchWeights[njets][nb][nht];
}
