#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <bitset>
#include <map>
#include <random>
#include <limits>
#include <chrono>
#include <stdlib.h>
#include <iomanip>

#define BIT 7
#define CHARS ((1 << BIT) - 1)

using namespace std;
using ull = unsigned long long;

extern const vector<wchar_t> symbols;
extern const vector<wchar_t> alphabet;
extern map<wchar_t, bitset<BIT> > mpCharBit;
extern unordered_map<bitset<BIT>, wchar_t > mpBitChar;
bool compare(ull a, ull b);
// Multimap lưu trữ số lần xuất hiện của các cặp bigrams
extern multimap<ull, pair<wchar_t, wchar_t>, decltype(&compare)> map_bigrams;
// Ma trận lưu trữ các trường hợp mã hóa có thể xảy ra cho 1 ký tự
extern vector<array<wchar_t, 128> > climbMatrix;

extern ull Count;
extern ull bigramNumber;
extern ull Sum;

bool binSearch(const vector<wchar_t> &symbols, wchar_t x);
void ftext_to_fbin(const string& fileInput, const string& fileOutput);
void toVieAlpb(const string& fileInput, const string& fileOutput);
void gen_vernam_key(const string& fileOutput, ull Count, int pos1, int pos2, int pos3);
void gen_vernam_key_rand(const string& fileOutput, ull Count, const int &posNum);
void vernam_crypt(map<wchar_t, bitset<BIT> >& mpi, unordered_map<bitset<BIT>, wchar_t>& mpo,
                  const string& filePlain, const string& fileKey, const string& fileOut, bool mode);
void freq_analysis(const string& fileInput, const string& fileOutput, ull &Count, ull &Sum);
void freq_analysis_key(const string& fileInput, const string& fileOutput);
void countBigrams(multimap<ull, pair<wchar_t, wchar_t>, decltype(&compare)>& map_bigrams, const string& fileIn, const string& fileOut);
void assignMap(const size_t &stride, map<wchar_t, bitset<BIT> > &tmp);
void codeAnalysis(const string& fileInput, const string& fileOutput, ull &Count, ull &Sum);
void decrypt(const string &codeFile, const string &outFile, vector<array<wchar_t, 128> >& climbMatrix);
void hillClimbing(const multimap<ull, pair<wchar_t, wchar_t>, decltype(&compare)> &map_bigrams,
                    const vector<array<wchar_t, 128> >& climbMatrix);

#endif // HEADER_H_INCLUDED
