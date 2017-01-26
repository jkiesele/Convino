/*
 * graphCosmetics.h
 *
 *  Created on: 28 Oct 2016
 *      Author: jkiesele
 */

#ifndef INCLUDE_GRAPHCOSMETICS_H_
#define INCLUDE_GRAPHCOSMETICS_H_

#include "TString.h"

enum graphcosmetics {gc_nominal, gc_scancombined, gc_minchi2};

class TGraphAsymmErrors;

void applyGraphCosmetics(TGraphAsymmErrors* g, graphcosmetics, double xlow, double xhigh, const TString & name, const TString & scanned, double scaler=1);

#endif /* INCLUDE_GRAPHCOSMETICS_H_ */
