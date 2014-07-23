/*
 *    This file is part of HMMoC 0.5, a hidden Markov model compiler.
 *    Copyright (C) 2006 by Gerton Lunter, Oxford University.
 *
 *    HMMoC is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    HMMOC is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with HMMoC; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
\*/
#include <cstdlib>
#include <cstring>
#include "homology.h"


void run(std::string& sequence, std::string& prediction, const Params& params ) 
{

  // The parameters of the model
  Params iPar = params;

  //
  // Next, build an input emission sequence by sampling the emitted symbols according to true path
  //

  int iPathLength = sequence.length() ;
  char* aSequence = new char[ iPathLength ];
  memcpy(aSequence, sequence.data(), iPathLength );

  // Decode the emission sequence using Viterbi, and compute posteriors and Baum Welch counts using Forward and Backward
  HomologyDPTable *pViterbiDP, *pFWDP, *pBWDP;
  HomologyBaumWelch bw;

  bfloat iFWProb = Forward(&pFWDP, iPar, aSequence, iPathLength );
  bfloat iBWProb = Backward(bw, pFWDP, &pBWDP, iPar, aSequence, iPathLength );

  prediction.resize(iPathLength);
  for (int i=0; i<iPathLength; i++) {

    double iPosterior = pFWDP->getProb("homologous",i+1)*pBWDP->getProb("homologous",i+1)/iFWProb;
//    if (iViterbiPath.toState(i) == iVHomologous) {
    if (iPosterior >= 0.9) {
      prediction[i] = 'H';
    } else {
      prediction[i] = 'N';
    }
//    cout << " " << iPosterior << endl;

  }
  //clean up aSequence, does this do any good? 
  delete[] aSequence;
  delete pFWDP;
  delete pBWDP;

}


