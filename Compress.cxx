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
#include "Timing.hpp"
#include "TH1D.h"
#include "TInterpreter.h"
#include <vector>
#include <numeric>
#include <iostream>
#include <gperftools/profiler.h>
#endif

std::map<int, size_t> FibNumbers;
TH1D* timing = nullptr;
TH1D* diff_histo = nullptr;
bool standard_fibonacci;
   
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
inline void add_to_sequence_terminate(std::vector<bool> array, std::vector<bool>& cmp) {
  cmp.insert(cmp.end(), array.begin(), array.end());
  cmp.push_back(1);
}

std::vector<short> CompressThis(const std::vector<short> wf, std::function<void(int, std::vector<std::vector<bool>>&)> add_to_table) {
  // The format it has to be in the end
  static std::vector<std::vector<bool>> table;
  MeasureTimeFrom();
  std::vector<short> comp_short;

  // First numbers are not encoded (size and baseline)
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
  std::vector<bool> cmp;

  // The input is changed to be the difference between ticks
  std::vector<short> diff;
  diff.reserve(wf.size());
  // Fibonacci number are positive, so need a way to get negative number
  // We use the ZigZag approach (even numbers are positive, odd are negative)
  std::adjacent_difference(wf.begin(), wf.end(),std::back_inserter(diff), [](const short& it1, const short& it2){
                                                                            short d = it1-it2;
                                                                            if (d > 0) d =  2 * d;
                                                                            else       d = -2 * d + 1;
                                                                            return d;
                                                                          });

  // Start from 1 to avoid the first number
  for (size_t iSample=1; iSample<diff.size(); ++iSample) {

    short d = diff[iSample];
    diff_histo->Fill(d);

    try{
      add_to_sequence_terminate(table.at(d), cmp);
    } catch(std::out_of_range& e) {
      int end = d+1;
      add_to_table(end, table);
      add_to_sequence_terminate(table.at(d), cmp);
    }
  }
  
  // Now convert all this to a vector of short and it's just another day in paradise for larsoft
  size_t n_vector = cmp.size();
  while (n_vector>sizeof(short)*8) {
    // Create a bitset of the size of the short to simplify things
    std::bitset<8*sizeof(short)> this_encoded;
    // Set the bitset to match the stream
    for (size_t it=0; it<8*sizeof(short); ++it) {
      if (cmp[it]) this_encoded.set(it);
    }

    // Get rid of stuff we've dealt with 
    cmp.erase(cmp.begin(), cmp.begin()+8*sizeof(short));
    // Cast the bitset to short
    short comp_s = (short)this_encoded.to_ulong();
    
    // Store the short in the output waveform
    comp_short.push_back(comp_s);
    n_vector = cmp.size();
  }
  
  // Deal with the last part
  std::bitset<8*sizeof(short)> this_encoded;
  for (size_t it=0; it<cmp.size(); ++it) {
    if(cmp[it]) this_encoded.set(it);
  }
  short comp_s = (short)this_encoded.to_ulong();
  comp_short.push_back(comp_s);
  timing->Fill(MeasureMicroSecondTime());

  return comp_short;
}

std::vector<short> UncompressThis(const std::vector<short> comp_short) {
  // First compressed sample is the size
  size_t n_samples = comp_short[0];

  // The output vector
  std::vector<short> wf;
  // The second compressed sample is the first uncompressed sample
  wf.push_back(comp_short[1]);
  
  // The thing that we want to decode (rather than jumbled short vector)
  std::vector<bool> comp;
  
  for (size_t i=2; i<comp_short.size(); ++i) {

    std::bitset<8*sizeof(short)> this_encoded(comp_short[i]);
    for (size_t i2=0; i2<this_encoded.size(); ++i2) {
      comp.push_back(this_encoded[i2]);
    }
  }
  
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
          
        zigzag_number += current_number[it2] * FibNumbers[it2+1];
      }
      
      short decoded = 0;
      if (zigzag_number%2 == 0) decoded = zigzag_number / 2;
      else                      decoded = -(zigzag_number - 1) / 2;
      short baseline = wf.back();
      
      wf.push_back(baseline + decoded);
      
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
  timing= new TH1D("timing", ";time us", 60, 0, 30000);
  diff_histo = new TH1D("difgf", ";diff", 50, -25, 25);
  FibNumbers[0] = 1;
  FibNumbers[1] = 1;
  for (int i=2; i<25; ++i) {
    FibNumbers[i] = FibNumbers[i-1] + FibNumbers[i-2];
  }
  FibNumbers.erase(0);
  // FibNumbers[0] = 0;
  standard_fibonacci = true;
  
  // std::cout << "Printing FibNumbers size: " << FibNumbers.size() << "\n";
  // for (auto const & it: FibNumbers)
  //   std::cout << "FibNumbers[" << it.first << "] = " << it.second << "\n";
  // std::cout << "Done\n";
  
  // for (size_t i=0; i<waveforms->size(); ++i) {
  TH1D* original_decode  = new TH1D("diff", "Original-Decoded;Original-Decode(Encoded);N ticks", 41, -20, 20);
  TH1D* compression_fact = new TH1D("Compression_Factor", "Original size / encoded size;Original size/Encoded size;N waveforms (2.2 ms)", 200, 0, 5);


  std::function<void (int, std::vector<std::vector<bool>>&)> add_to_table;
  
  add_to_table = [](int end, std::vector<std::vector<bool>>& table) -> void {
                   if (table.size() > end)
                     return;

                   std::map<int, int> fibn;
                   fibn[0] = 0;
                   fibn[1] = 1;
                   fibn[2] = 1;

                   if (table.size()==0) {
                     table.push_back(std::vector<bool>());
                     table.push_back({true      });
                     table.push_back({false,true});
                   }

                   // start at the third fibonacci number
                   for (int i=2; fibn.rbegin()->second<end; ++i) {
                     fibn[i] = fibn[i-1] + fibn[i-2];
                   }
                   

                   for (int i=table.size(); i<=end; ++i) {
                     int current_number = i;
                     std::vector<bool> seq; // the final sequence of numbers
                     while (current_number>0) {
                       // find the biggest fibonacci number that is smaller than the current number
                       for (auto fib = fibn.rbegin(); fib != fibn.rend(); ++fib) {
                         if (fib->second<=current_number) {
                           // if this is the first number, create the vector
                           if (!seq.size())
                             seq=std::vector<bool>(fib->first-1, false);
                           
                           seq[fib->first-2] = true;
                           current_number-=fib->second;
                           break;
                         }
                       }
                     }
                     table.push_back(seq);
                   }
                 };
                 
  
//  for (int i=0; i<waveforms->size()/100; ++i) {
  //for (int i=0; i<1; ++i) {
   for (int i=0; i<100; ++i) {

    if (i%100==0) {
      double progress = 100. * i / waveforms->size();
      std::cout << "progress: " << progress << "%\n";
      // if (progress > 1) break;
    }
    std::vector<short> wf;
    wf.reserve((*waveforms)[i].size());
    size_t size=0;
    for (auto const& adc: (*waveforms)[i]) {
      if (size++>1000000) break;
      wf.push_back(adc);
    }
    //timing->Fill();
    // std::vector<short> compressed = CompressThis(wf, i==9);
    // std::vector<short>  uncompressed = UncompressThis(compressed, i==9);
    std::vector<short> compressed   = CompressThis(wf, add_to_table);
    std::vector<short> uncompressed = UncompressThis(compressed);

    for (size_t it=0; it<wf.size(); ++it) {
      double diff = wf[it] - uncompressed[it];
      //   if (diff != 0) {
      //     std::cout << "\n\n\ni " << i << "\n";
      //     std::cout << "it : " << it << "\n";
      //     std::cout << "wf[it] : " << wf[it] << "\n";
      //     std::cout << "uncompressed[it] : " << uncompressed[it] << "\n";
      //     std::cout << "Printing compressed size: " << compressed.size() << "\n";
      //     for (size_t it=0; it<4; ++it)
      //       std::cout << "compressed[" << it << "] = " << compressed[it] << "\n";
      //     std::cout << "Done\n";

      //     throw;
      //   }
      
      //   if      (diff >  20) diff =  19;
      //   else if (diff < -20) diff = -19;
      original_decode->Fill(diff);
     }

    compression_fact->Fill((double)wf.size()/compressed.size());
  }

  TCanvas c;
  c.Print("compression.pdf[");
  gPad->SetLogy(1);
  original_decode->SetStats(0);
  original_decode->SetMinimum(0.1);
  original_decode->SetBinContent(1, original_decode->GetBinContent(0)+original_decode->GetBinContent(1));
  original_decode->SetBinContent(original_decode->GetNbinsX(), original_decode->GetBinContent(original_decode->GetNbinsX()+1)+original_decode->GetBinContent(original_decode->GetNbinsX()));
  original_decode->Draw();
  c.Print("compression.pdf");
  gPad->SetLogy(0);
  // compression_fact->Draw();
  // c.Print("compression.pdf");
  timing->Draw();
  c.Print("compression.pdf");
  diff_histo->Draw();
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
