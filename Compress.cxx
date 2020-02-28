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
inline void add_to_sequence_terminate(std::vector<bool> array, std::vector<bool> sqce, std::vector<bool>& cmp) {
  if (sqce.size()) {
    std::copy(array.begin(), array.end(), sqce.begin());
    cmp.insert(cmp.end(), sqce.begin(), sqce.end());
  } else {
    cmp.insert(cmp.end(), array.begin(), array.end());
  }
  cmp.push_back(1);
}

std::vector<short> CompressThis(const std::vector<short> wf, bool verb=false) {
  // ProfilerStart("prof.out");
  // The format it has to be in the end
  static std::map<int,std::vector<bool>> table;
  static std::map<int,std::vector<bool>> table;
  table[1 ] = {1            };
  table[2 ] = {0,1          };
  table[3 ] = {0,0,1        };
  table[4 ] = {1,0,1        };
  table[5 ] = {0,0,0,1      };
  table[6 ] = {1,0,0,1      };
  table[7 ] = {0,1,0,1      };
  table[8 ] = {0,0,0,0,1    };
  table[9 ] = {1,0,0,0,1    };
  table[10] = {0,1,0,0,1    };
  table[11] = {0,0,1,0,1    };
  table[12] = {1,0,1,0,1    };
  table[13] = {0,0,0,0,0,1  };
  table[14] = {1,0,0,0,0,1  };
  table[15] = {0,1,0,0,0,1  };
  table[16] = {0,0,1,0,0,1  };
  table[17] = {1,0,1,0,0,1  };
  table[18] = {0,0,0,1,0,1  };
  table[19] = {1,0,0,1,0,1  };
  table[20] = {0,1,0,1,0,1  };
  // table[21] = {0,0,0,0,0,0,1};
  // table[22] = {1,0,0,0,0,0,1};
  // table[23] = {0,1,0,0,0,0,1};
  // table[24] = {0,0,1,0,0,0,1};
  // table[25] = {1,0,1,0,0,0,1};
  // table[26] = {0,0,0,1,0,0,1};
  // table[27] = {1,0,0,1,0,0,1};
  // table[28] = {0,1,0,1,0,0,1};
  // table[29] = {0,0,0,0,1,0,1};
  // table[30] = {1,0,0,0,1,0,1};
  // table[31] = {0,1,0,0,1,0,1};
  
  MeasureTimeFrom();
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
  //std::vector<bool> comp;
  std::vector<bool> cmp;

  // The input is changed to be the difference between ticks
  std::vector<short> diff;
  diff.reserve(wf.size());
  std::adjacent_difference(wf.begin(), wf.end(), std::back_inserter(diff));

  // Start from 1 to avoid the first number
  for (size_t iSample=1; iSample<diff.size(); ++iSample) {
    // std::cout << "iSample : " << iSample << "\n";

    //for (size_t iSample=0; iSample<31; ++iSample) {
    short d = diff[iSample];
    diff_histo->Fill(d);
    // Store that for posterity and checks
    const short d_orig_nozigzag = d;

    // Fibonacci number are positive, so need a way to get negative number
    // We use the ZigZag approach (even numbers are positive, odd are negative)
    if (d > 0) d =  2 * d;
    else       d = -2 * d + 1;

    // Store this zigzag number for posterity
    //const short d_orig = d;

    // This vector stores the indices of Fibonacci number used
    //std::vector<int> fib_sequence;
    //std::vector<bool> sqce;
    try{
      add_to_sequence_terminate(table.at(d), std::vector<bool>(), cmp);
    } catch(std::out_of_range& e) {
      int start = table.size();
      int end   = d;
      std::map<int, int> fibn;
      fibn[1] = 1;
      fibn[2] = 1;

      for (int iF=3; iF<=end; ++iF) {
        fibn[iF] = fibn[iF-1] + fibn[iF-2];
        if (iF>start) {
          int current_number = iF;
          std::vector<bool> seq;
          while (current_number>0) {
            for (auto fib = fibn.rbegin(); fib != fibn.rend(); ++fib) {
              if (fib->second<=current_number) {
                if (!seq.size())
                  seq=std::vector<bool>(fib->first-1, false);
                seq[fib->first-2] = true;
                current_number-=fib->second;
                break;
              }
            }
          }
          table.insert(table.end(), std::pair<int, std::vector<bool>>(iF, seq));
        }
      }
      add_to_sequence_terminate(table.at(d), std::vector<bool>(), cmp);
    }
  }
  
  timing->Fill(MeasureMicroSecondTime());
  // Now convert all this to a vector of short and it's just another day in paradise for larsoft
  size_t n_vector = cmp.size();
  while (n_vector>sizeof(short)*8) {
    // Create a bitset of the size of the short to simplify things
    std::bitset<8*sizeof(short)> this_encoded;
    // Set the bitset to match the stream
    if (verb) std::cout <<"[ ";
    for (size_t it=0; it<8*sizeof(short); ++it) {
      if (cmp[it]) this_encoded.set(it);
      if (verb) std::cout << cmp[it] << " ";
    }
    if (verb) std::cout <<"]\n";
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
  timing= new TH1D("timing", ";time ns", 60, 0, 30000);
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
  
  for (int i=0; i<waveforms->size()/100; ++i) {
  // for (int i=0; i<100; ++i) {

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
    std::vector<short> compressed   = CompressThis(wf);
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
  // gPad->SetLogy(1);
  // original_decode->SetStats(0);
  // original_decode->SetMinimum(0.1);
  original_decode->SetBinContent(1, original_decode->GetBinContent(0)+original_decode->GetBinContent(1));
  original_decode->SetBinContent(original_decode->GetNbinsX(), original_decode->GetBinContent(original_decode->GetNbinsX()+1)+original_decode->GetBinContent(original_decode->GetNbinsX()));
  original_decode->Draw();
  c.Print("compression.pdf");
  // gPad->SetLogy(0);
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
