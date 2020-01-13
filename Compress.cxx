// $Id$ 
// Author: Pierre Lasorak <plasorak@Pierres-MBP>
// Update: 2020-01-09 09:58:39+0000
// Copyright: 2020 (C) Pierre Lasorak

#ifndef __CINT__
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TInterpreter.h"
#include <vector>
#include <numeric>
#include <iostream>
#endif

std::map<int, size_t> FibNumbers;


// template<typename T>
// std::vector<bool> EncodeToBool(T number) {
//   std::vector<bool> encoded;
//   // std::cout << "number : " << number << "\n";

//   T n = number;
//   for (int i=sizeof(T)*8-1; i>=0; --i) {
//     encoded.insert(encoded.begin(),floor(n / TMath::Power(2,i)));
//     n -= TMath::Power(2,i) * encoded.front();
//     // std::cout << "i : " << i << "\t" << encoded.back() << "\n";
//   }
  
//   return encoded;
// }

// template<typename T>
// T DecodeFromBool(std::vector<bool> encoded) {
//   std::bitset<sizeof(T)*8> bits;
//   for (int i=0; i<sizeof(T)*8; ++i) {
//     bits[i] = encoded[i];
//   }
//   // std::cout << "bits : " << bits << "\n";
//   short r = (short)bits.to_ulong();
//   // std::cout << "r : " << r << "\n";

//   return r;
// }

std::vector<short> CompressThis(const std::vector<short> wf, bool verb=false) {

  // The format it has to be in the end
  std::vector<short> comp_short;
  // First number is not encoded (baseline)
  if (wf.size() > std::numeric_limits<short>::max()) {
    std::cout << "WOW! You are trying to compress a " << wf.size() << " long waveform";
    std::cout << " in a vector of shorts.\nUnfortunately, the encoded waveform needs to store the size";
    std::cout << " of the wavefom (in its first number) and " << wf.size() << " is bigger than the maximum\n";
    std::cout << " of the shorts (" << std::numeric_limits<short>::max() << ").\nBailing out disgracefully to avoid massive trouble.\n";
    throw;
  }
  comp_short.push_back(wf.size());
  comp_short.push_back(*wf.begin());

  // The format we are working with
  std::vector<bool> comp;

  // The input is changed to be the difference between ticks
  std::vector<short> diff;
  diff.reserve(wf.size());
  std::adjacent_difference(wf.begin(), wf.end(), std::back_inserter(diff));

  // Start from 1 to avoid the first number
  for (size_t iSample=1; iSample<diff.size(); ++iSample) {
    short d = diff[iSample];

    // Store that for posterity and checks
    const short d_orig_nozigzag = d;

    // Fibonacci number are positive, so need a way to get negative number
    // We use the ZigZag approach (even numbers are positive, odd are negative)
    if (d > 0) d =  2 * d;
    else       d = -2 * d + 1;

    // Store this zigzag number for posterity
    const short d_orig = d;

    // This vector stores the indices of Fibonacci number used
    std::vector<int> fib_sequence;

    while (true) {
      short max_fib = 0;
      size_t iFib=0;

      if (d == 1) {
        fib_sequence.push_back(1);
        break;
      } else if (d == 0) {
        break;
      }

      // Find the first Fibonacci number that exceed the number
      for (; iFib<FibNumbers.size(); ++iFib) {
        if ((size_t)d < FibNumbers[iFib]) {
          break;
        }
      }
      // And remove it from the initial number
      d = d - FibNumbers[iFib-1];
      // Store the previous Fibonacci number
      fib_sequence.push_back(iFib-1);

      // And continue until you get 0 or 1
    }

    short sum=0;
    std::vector<bool> including;

    // Little cross check
    for (size_t it=0; it<fib_sequence.size(); ++it) {
      // std::cout << "Fibonacci number index: " << fib_sequence[it] << " (which is equal to " << FibNumbers[fib_sequence[it]] << ")\n";
      sum += FibNumbers[fib_sequence[it]];
    }

    // fib_sequence are in decreasing order
    // Go from 0 to the last number and check for each of them if the number is in the sequence
    for (int it=1; it<=(*fib_sequence.begin()); ++it) {
      auto found = std::find(fib_sequence.begin(),fib_sequence.end(),it);
      // If it is, store a 1
      if (found != fib_sequence.end()) {
        including.push_back(1);
        comp.push_back(1);
      // If it isn't, store a 0
      } else {
        including.push_back(0);
        comp.push_back(0);
      }
    }
    // Always add a 1 at the end.
    including.push_back(1);
    comp.push_back(1);
    if (verb) {
      std::cout << d_orig_nozigzag << "\t-> zigzag version: " << d_orig << "\t-> encoded version: " << including.size() << " bits [ ";
      for (size_t it=0; it<including.size(); ++it)
        std::cout << including[it] << " ";
      std::cout << "]\n";
    }
    // Some cross-checks
    if (sum != d_orig) {
      std::cerr << "#$%*#J@@$(S!!!\n";
      std::cerr << "diff[" << iSample << "]: " << diff[iSample] << "\t";
      std::cerr << "d: " << d_orig << "\n";
      std::cerr << "fib sequence size : " << fib_sequence.size() << "\n";
      short sum=0;
      for (size_t it=0; it<fib_sequence.size(); ++it) {
        std::cerr << fib_sequence[it] << " -> " << FibNumbers[fib_sequence[it]] << "\n";
      }
      throw;
    }
  }

  // Now convert all this to a vector of short and it's just another day in paradise for larsoft
  size_t n_vector = comp.size();
  while (n_vector>sizeof(short)*8) {
    // Create a bitset of the size of the short to simplify things
    std::bitset<8*sizeof(short)> this_encoded;
    // Set the bitset to match the stream
    if (verb) std::cout <<"[ ";
    for (size_t it=0; it<8*sizeof(short); ++it) {
      if (comp[it]) this_encoded.set(it);
      if (verb) std::cout << comp[it] << " ";
    }
    if (verb) std::cout <<"]\n";
    // Get rid of stuff we've dealt with 
    comp.erase(comp.begin(), comp.begin()+8*sizeof(short));
    // Cast the bitset to short
    short comp_s = (short)this_encoded.to_ulong();
    
    // Store the short in the output waveform
    comp_short.push_back(comp_s);
    n_vector = comp.size();
  }
  
  // Deal with the last part
  std::bitset<8*sizeof(short)> this_encoded;
  for (size_t it=0; it<comp.size(); ++it) {
    if(comp[it]) this_encoded.set(it);
  }
  short comp_s = (short)this_encoded.to_ulong();
  comp_short.push_back(comp_s);
  // std::cout << "yay\n";
  // std::cout << "Printing comp_short size: " << comp_short.size() << "\n";
  // for (size_t it=0; it<comp_short.size(); ++it) 
  //   std::cout << "comp_short[" << it << "] = " << comp_short[it] << "\n";
  // std::cout << "Done\n";

  return comp_short;
}

std::vector<short> UncompressThis(const std::vector<short> comp_short, bool verb=false) {
  // First compressed sample is the size
  size_t n_samples = comp_short[0];

  // The output vector
  std::vector<short> wf;
  // The seecond compressed sample is the first uncompressed sample
  wf.push_back(comp_short[1]);
  //std::cout << "wf[0] : " << wf[0] << "\n";
  if (verb) std::cout << "Decoding\n";

  // The thing that we want to decode (rather than jumbled short vector)
  std::vector<bool> comp;
  if(verb) std::cout << "[ ";
  for (size_t i=2; i<comp_short.size(); ++i) {

    std::bitset<8*sizeof(short)> this_encoded(comp_short[i]);
    for (size_t i2=0; i2<this_encoded.size(); ++i2) {
      if (verb) std::cout << this_encoded[i2] << " ";
      comp.push_back(this_encoded[i2]);
    }
  }
  if (verb) std::cout <<"]\n";

  // The bit which has to be decoded ("chunk")
  std::vector<bool> current_number;
  for (size_t it=0; it<comp.size(); ++it) {
    current_number.push_back(comp[it]);
    // If we have at least 2 numbers in the current chunk
    // and if the last 2 number are one
    // and if we are not at the end
    if ((current_number.size()>=2 && it >= 1 && comp[it-1] == 1 && comp[it] == 1) || it == comp.size()-1) {
      short zigzag_number = 0;
      for (size_t it2=0; it2<current_number.size()-1; ++it2) {
        if (current_number[it2])
          //if (verb) std::cout << "Fibonacci number index: " << it2 << " (which is equal to " << FibNumbers[it2] << ")\n";
        zigzag_number += current_number[it2] * FibNumbers[it2+1];
      }
      if (verb) {
        std::cout << current_number.size() << " bits [ ";
        for (size_t it2=0; it2<current_number.size(); ++it2) {
          std::cout << current_number[it2] << " ";
        }
        std::cout << " ]\n";
      }
      short decoded = 0;
      if (zigzag_number%2 == 0) decoded = zigzag_number / 2;
      else                      decoded = -(zigzag_number - 1) / 2;
      short baseline = wf.back();
      // std::cout << "baseline : " << baseline << "\n";
      wf.push_back(baseline + decoded);
      // std::cout << "]\t-> decoded zigzag: " << zigzag_number << "\t-> decoded: " << decoded<< "\n";
      current_number.clear();
      if (wf.size() == n_samples) break;
    }
  }
  
  return wf;
}


int Compress()
{
  // DEFINE YOUR MAIN FUNCTION HERE
  TFile f("~/Desktop/full_waveform.root");
  TTree* t=(TTree*)f.Get("waveformtree/WaveformTree");
  std::vector<std::vector<int>>* waveforms = nullptr;
  std::vector<int>* channels = nullptr;
  gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
  t->SetBranchAddress("waveforms", &waveforms);
  t->SetBranchAddress("chans",  &channels);
  t->GetEntry(0);
  
  FibNumbers[0] = 1;
  FibNumbers[1] = 1;
  for (int i=2; i<25; ++i) {
    FibNumbers[i] = FibNumbers[i-1] + FibNumbers[i-2];
  }
  FibNumbers.erase(0);
  // FibNumbers[0] = 0;
  
  // std::cout << "Printing FibNumbers size: " << FibNumbers.size() << "\n";
  // for (auto const & it: FibNumbers)
  //   std::cout << "FibNumbers[" << it.first << "] = " << it.second << "\n";
  // std::cout << "Done\n";
  
  // for (size_t i=0; i<waveforms->size(); ++i) {
  TH1D* original_decode  = new TH1D("diff", "Original-Decoded;Original-Decode(Encoded);N ticks", 20, -20, 20);
  TH1D* compression_fact = new TH1D("Compression_Factor", "Original size / encoded size;Original size/Encoded size;N waveforms (2.2 ms)", 200, 0, 5);
  
  for (int i=0; i<waveforms->size(); ++i) {
    // for (int i=0; i<20; ++i) {
     if (i%100==0) {
       double progress = 100. * i / waveforms->size();
       std::cout << "progress: " << progress << "%\n";
       // if (progress > 1) break;
     }
    std::vector<short> wf;
    wf.reserve((*waveforms)[i].size());
    size_t size=0;
    for (auto const& adc: (*waveforms)[i]) {
      wf.push_back(adc);
      if (size++>5000000) break;
    }
    
    // std::vector<short> compressed = CompressThis(wf, i==9);
    // std::vector<short>  uncompressed = UncompressThis(compressed, i==9);
    std::vector<short> compressed = CompressThis(wf);
    std::vector<short>  uncompressed = UncompressThis(compressed);

    for (size_t it=0; it<wf.size(); ++it) {
      double diff = wf[it] - uncompressed[it];
      if (diff != 0) {
        std::cout << "i " << i << "\n";
        std::cout << "it : " << it << "\n";
        std::cout << "wf[it] : " << wf[it] << "\n";
        std::cout << "uncompressed[it] : " << uncompressed[it] << "\n";
        std::cout << "Printing compressed size: " << compressed.size() << "\n";
        for (size_t it=0; it<4; ++it)
          std::cout << "compressed[" << it << "] = " << compressed[it] << "\n";
        std::cout << "Done\n";

        throw;
      }
      
      if (diff > 20) diff = 19;
      else if (diff < -20) diff=-19;
      original_decode->Fill(diff);
    }

    compression_fact->Fill((double)wf.size()/compressed.size());
  }

  TCanvas c;
  c.Print("compression.pdf[");
  gPad->SetLogy(1);
  original_decode->SetStats(0);
  original_decode->SetMinimum(0.1);
  original_decode->Draw();
  c.Print("compression.pdf");
  gPad->SetLogy(0);
  compression_fact->Draw();
  c.Print("compression.pdf");
  c.Print("compression.pdf]");
  
  
  return 0;
}

#ifndef __CINT__
int main()
{
  return Compress();
}
#endif

//____________________________________________________________________ 
//  
// EOF
//
