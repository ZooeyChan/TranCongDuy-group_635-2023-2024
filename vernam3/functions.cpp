#include "header.h"

ull Count = 0;
ull Sum = 0;
ull bigramNumber = 0;
multimap<ull, pair<wchar_t, wchar_t>, decltype(&compare)> map_bigrams;
vector<array<wchar_t, 128> > climbMatrix;

bool binSearch(const vector<wchar_t> &symbols, wchar_t x) {
    ull l = 0, r = symbols.size() - 1;
    while (l <= r) {
        ull m = (l + r) / 2;
        if (symbols[m] == x) {
            return true;
            }
        else if (symbols[m] < x) {
            l = m + 1;
            }
        else {
            r = m - 1;
            }
        }
    return false;
    }

void toVieAlpb(const string& fileInput, const string& fileOutput) {
    wifstream in(fileInput);
    wofstream out(fileOutput, ios::trunc);
    if (!in || !out) {
        cerr << "Cannot open file!!!" << endl;
        return;
        }
    wchar_t ch;
    try {
        while (in.get(ch)) {
            switch (tolower(ch)) {
                case L'w':
                    out << L'v' << L'v';
                    break;
                case L'f':
                    out << L'p' << L'h';
                    break;
                case L'z':
                    out.put(L'd');
                    break;
                case L'j':
                    out << L'g' << L'i';
                    break;
                default:
                    out.put(ch);
                    break;
                }
            }
        }
    catch (const exception& e) {
        cerr << "Error during file processing: " << e.what() << endl;
        }
    in.close();
    out.close();
    }

//void gen_vernam_key_pcg(const string& fileOutput, ull Count, const int &posNum) {
//    ofstream out(fileOutput, ios::binary | ios::trunc);
//    pcg32 rng(pcg_extras::seed_seq_from<random_device>());
//    uniform_int_distribution<uint8_t> distribution(0, CHARS);
//    uint8_t ch;
//    int temp;
//    int pos[posNum];
//    // Sinh số ngẫu nhiên cho vị trí bit
//    for (int i = 0; i < posNum; i++) {
//        pos[i] = rng() % 7;
//        cout << pos[i] << " ";
//        }
//    // Sinh khóa ngẫu nhiên
//    while (Count--) {
//        ch = distribution(rng);
//        temp = posNum;
//        for (int i = 0; i < posNum; i++) {
//            ch &= ~(1 << pos[i]);
//            }
//        out.write(reinterpret_cast<const char*>(&ch), sizeof(uint8_t));
//        }
//    out.close();
//    }

void gen_vernam_key_rand(const string& fileOutput, ull Count, const int &posNum) {
    ofstream out(fileOutput, ios::binary | ios::trunc);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint8_t> distribution(0, CHARS);
    uint8_t ch;
    int temp;
    int pos[posNum];
    for(int i = 0; i < posNum; i++) {
        pos[i] = gen() % 7; // Generate random numbers from 0 to 6
        cout << pos[i] << " ";
        }
    while (Count--) {
        ch = distribution(gen);
        temp = posNum;
        // Set the bit at posNum (position) to 0
        for (int i = 0; i < posNum; i++) {
            ch &= ~(1 << pos[i]);
            }
        out.write(reinterpret_cast<const char*>(&ch), sizeof(uint8_t));
        }
    out.close();
    }

void vernam_crypt(map<wchar_t, bitset<BIT> >& mpi, unordered_map<bitset<BIT>, wchar_t>& mpo,
                  const string& filePlain, const string& fileKey, const string& fileOut, bool mode) {
    wifstream plain(filePlain);
    wifstream key(fileKey);
    wofstream out(fileOut, ios::trunc);
    if (!plain || !key || !out) {
        cerr << "Can not open file!!!" << endl;
        return;
        }
    try {
        wchar_t pl;
        wchar_t k;
        while (plain.get(pl) && key.get(k)) {
            if(pl == L'\n') {
                out.put(L'\n');
                continue;
                }
            bitset<BIT> bits = mpi[towlower(pl)];
            bits ^= static_cast<bitset<BIT> >(k);
            wchar_t result = mpo[bits];
            out.put(result);
            }
        plain.close();
        key.close();
        out.close();
        }
    catch (const exception& e) {
        cerr << "Error during file processing: " << e.what() << endl;
        plain.close();
        key.close();
        out.close();
        }
    }

void freq_analysis(const string& fileInput, const string& fileOutput, ull &Count, ull &Sum) {
    map<wchar_t, ull> mp;
    wifstream inputFile(fileInput);
    if (!inputFile.is_open()) {
        cerr << "Input error!!!" << endl;
        return;
        }
    wofstream outputFile(fileOutput, ios::trunc);
    if (!outputFile.is_open()) {
        cerr << "Output error!!!" << endl;
        inputFile.close();
        return;
        }
    wstring str;
    while (getline(inputFile, str)) {
        for (wchar_t c : str) {
            if (binSearch(symbols, tolower(c))) {
                mp[tolower(c)]++;
                Sum++;
                if(c != L' ') Count++;
                }
            }
        }
    for (wchar_t c : symbols) {
        outputFile << c << " : "<< setw(10) << mp[c] << "   " <<
                   setprecision(10) << (static_cast<double>(mp[c]) / Count)* 100 << endl;
        }
    outputFile<< "\nTotal number of characters without space: "<< Count;
    outputFile<< "\nTotal number of characters: "<< Sum;
    inputFile.close();
    outputFile.close();
    }

void freq_analysis_key(const string& fileInput, const string& fileOutput) {
    ull dem= 0;
    unordered_map<bitset<BIT>, ull> mp;
    ifstream inputFile(fileInput);
    if (!inputFile.is_open()) {
        cerr << "Input error!!!" << endl;
        return;
        }
    ofstream outputFile(fileOutput, ios::trunc);
    if (!outputFile.is_open()) {
        cerr << "Output error!!!" << endl;
        inputFile.close();
        return;
        }
    char ch;
    while (inputFile.get(ch)) {
        bitset<BIT> bits= static_cast<bitset<BIT> >(ch);
        mp[bits]++;
        dem++;
        }
    for (auto &p : mp) {
        outputFile << p.first << " : "<< p.second << " - " << (static_cast<double>(p.second) / dem) * 100 << endl;
        }
    outputFile<< "\nTotal number of key's characters: "<< dem;
    inputFile.close();
    outputFile.close();
    }

bool compare(ull a, ull b) {
    return a > b;
    }

void countBigrams(multimap<ull, pair<wchar_t, wchar_t>, decltype(&compare)>& map_bigrams, const string& fileIn, const string& fileOut) {
    int len = alphabet.size();
    ull arr[len][len] = {0};
    auto p0 = alphabet.end();
    wchar_t ch;
    wifstream in(fileIn);
    wofstream out(fileOut, ios::trunc);
    if (!in.is_open() || !out.is_open()) {
        cerr << "Can not open the file!!!" << endl;
        return;
        }
    while ((ch = in.get()) != EOF) {
        auto p1 = find(alphabet.begin(), alphabet.end(), tolower(ch));
        if (p1 != alphabet.end() && p0 != alphabet.end()) {
            arr[p0 - alphabet.begin()][p1 - alphabet.begin()]++;
            bigramNumber++;
            }
        p0 = p1;
        }
    in.close();
    for(int i = 0; i< len; i++) {
        for(int j = 0; j< len; j++) {
            map_bigrams.insert({arr[i][j], make_pair(alphabet[i], alphabet[j])});
            }
        }
    // Print header
    out << "        " << setw(9);
    for (int i = 0; i < len; i++) {
        out<< setw(9) << alphabet[i];
        }
    out << endl;
    // Print data
    for (int i = 0; i < len; i++) {
        out << setw(9) << alphabet[i];
        for (int j = 0; j < len; j++) {
            ull n = arr[i][j];
            out<< setw(9) << n;
            }
        out << endl;
        }
    out << "Total bigrams: " << bigramNumber << endl;
    multimap<ull, pair<wchar_t, wchar_t> >::iterator it = map_bigrams.begin();
    out << "Frequency of appearance" << "  ---->  " << "Bigram" << "          Appearance rate ( % )" << endl;
    while (it != map_bigrams.end()) {
        out<< "       " << setw(6) << it->first << "   ---------->     "
           << it->second.first << " - " << it->second.second << "  -------->  "
           << setprecision(10) << (static_cast<double>(it->first) / bigramNumber) * 100 << endl;
        ++it;
        }
    out.close();
    }

void assignMap(const size_t &stride, map<wchar_t, bitset<BIT> > &tmp) {
    size_t index = 0;
    while (index < symbols.size()) {
        for (size_t i = 0; i < stride && index < symbols.size(); ++i) {
            wchar_t key = symbols[index];
            tmp[key] = mpCharBit.at(key);
            ++index; // Increase index by 1 after each iteration
            }
        index += stride; // Bỏ qua stride phần tử tiếp theo
        }
    }

// Function to check if a code has been used
bool isUsed(const vector<array<wchar_t, 128> >& climbMatrix, wchar_t code, int index) {
    for (int i = 0; i < index; ++i) {
        if (climbMatrix[i][code] == code)
            return true;
        }
    return false;
    }

// Function to calculate the total number of occurrences of pairs of bigrams
ull computeScore(const vector<array<wchar_t, 128> >& climbMatrix,
                 const multimap<ull, pair<wchar_t, wchar_t>, decltype(&compare)>& map_bigrams) {
    ull score = 0;
    for (const auto& entry : map_bigrams) {
        wchar_t first = entry.second.first;
        wchar_t second = entry.second.second;
        if (climbMatrix[first][second] != first && climbMatrix[second][first] != second)
            score += entry.first;
        }
    return score;
    }

void decrypt(const string &codeFile, const string &outFile, vector<array<wchar_t, 128> >& climbMatrix) {
    map<wchar_t, bitset<BIT> > tmp;
//    size_t startKey = 0;
//    size_t endKey = 31;
//    for (size_t i = startKey; i <= endKey; ++i) {
//        wchar_t key = symbols[i];
//        tmp[key] = mpCharBit[key];
//        }
//    startKey = 63;
//    endKey = 95;
//    for (size_t i = startKey; i <= endKey; ++i) {
//        wchar_t key = symbols[i];
//        tmp[key] = mpCharBit[key];
//        }
    // stride = 2^pos, pos: position ignored in gamma.
    assignMap(32, tmp);
    wifstream in(codeFile);
    if (!in.is_open()) {
        cerr << "Input file error!!!" << endl;
        return;
        }
    wofstream out(outFile, ios::app);
    if (!out.is_open()) {
        cerr << "Output file error!!!" << endl;
        return;
        }
    wstring code;
    getline(in, code);
    vector<array<wchar_t, 128> > tempClimbMatrix(code.size());
    for(size_t i = 0; i < code.size(); i++) {
        tempClimbMatrix[i][0] = code[i];
        }
    for(size_t i = 0; i < code.size(); i++) {
        for(size_t j = 1; j < tmp.size(); j++) {
            auto it = next(tmp.begin(), j);
            tempClimbMatrix[i][j] = mpBitChar[mpCharBit[code[i]] ^ it->second];
            }
        }
    climbMatrix = tempClimbMatrix;
    for(size_t i = 0; i < tmp.size(); i++) {
        for(size_t j = 0; j < code.size(); j++) {
            out << setw(2) << tempClimbMatrix[j][i];
            }
        out << endl;
        }
    in.close();
    out.close();
    }

size_t round_off(double N, const double &n) {
    int h;
    double l, a, b, c, d, e, i, j, m, f, g;
    b = N;
    c = floor(N);
    for (i = 0; b >= 1; ++i) b = b / 10;
    d = n - i;
    b = N;
    b = b * pow(10, d);
    e = b + 0.5;
    if ((float)e == (float)ceil(b)) {
        f = (ceil(b));
        h = f - 2;
        if (h % 2 != 0) {
            e = e - 1;
            }
        }
    j = floor(e);
    m = pow(10, d);
    j = j / m;
    return j;
    }

void codeAnalysis(const string& fileInput, const string& fileOutput, ull &Count, ull &Sum) {
    map<wchar_t, ull> mp;
    wifstream inputFile(fileInput);
    if (!inputFile.is_open()) {
        cerr << "Input error!!!" << endl;
        return;
        }
    wofstream outputFile(fileOutput, ios::app);
    if (!outputFile.is_open()) {
        cerr << "Output error!!!" << endl;
        inputFile.close();
        return;
        }
    wstring str;
    while (getline(inputFile, str)) {
        for (wchar_t c : str) {
            if (binSearch(symbols, tolower(c))) {
                mp[tolower(c)]++;
                Sum++;
                if(c != L' ') Count++;
                }
            }
        }
    short dem = 1;
    size_t tmp = round_off(mp.begin()->second, 1);
    for (const auto &c : mp) {
        size_t current = round_off(double(c.second), 2);
        cout << current << endl;
        if (current != tmp) {
            tmp = current;
            dem++;
            }
        }
    outputFile << "\nThe missing bit in the Vernam key is at position: " << log2(pow(2, BIT) / dem);
    inputFile.close();
    outputFile.close();
    }

/*---------------------DEMO-----------------------------*/
// Plaintext fitness (dummy function)
double fitness(const string& plaintext) {
    // Count number of common words
    return static_cast<double>(count(plaintext.begin(), plaintext.end(), 'e'));
    }

// Decrypt the ciphertext using the key
string decipher(const string& ciphertext, const vector<string>& key) {
    string plaintext = ciphertext;
    int m = key.size();
    for (size_t i = 0; i < ciphertext.size(); ++i) {
        int pos = ciphertext[i] - 'A';
        plaintext[i] = key[i % m][pos];
        }
    return plaintext;
    }

// Generate a random key
vector<string> randomKey(int m) {
    vector<string> key(m, string(26, ' '));
    for (int i = 0; i < m; ++i) {
        string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        random_shuffle(alphabet.begin(), alphabet.end());
        key[i] = alphabet;
        }
    return key;
    }

// Function to implement the decoding algorithm
vector<string> hillClimbingAttack(const string& ciphertext, int m, int maxBigCounter) {
    srand(time(nullptr));
    // Initialize variables
    int bigCounter = 0;
    double bestFitness = fitness(ciphertext); // Fitness of the codex
    vector<string> parentKey = randomKey(m);
    vector<string> bestKey = parentKey;
    while (bigCounter < maxBigCounter) {
        for (int i = 0; i < m; ++i) {
            // Randomize the ith key alphabet
            vector<string> childKey = parentKey;
            string& keyAlphabet = childKey[i];
            random_shuffle(keyAlphabet.begin(), keyAlphabet.end());
            string plaintext = decipher(ciphertext, parentKey);
            double parentFitness = fitness(plaintext);
            int littleCounter = 0;
            while (littleCounter < 1000) {
                vector<string> tempKey = childKey;
                swap(tempKey[i][rand() % 26], tempKey[i][rand() % 26]);
                string tempPlaintext = decipher(ciphertext, tempKey);
                double childFitness = fitness(tempPlaintext);
                if (childFitness > parentFitness) {
                    parentKey = tempKey;
                    parentFitness = childFitness;
                    littleCounter = 0;
                    }
                else {
                    ++littleCounter;
                    }
                if (childFitness > bestFitness) {
                    bestFitness = childFitness;
                    bestKey = tempKey;
                    bigCounter = 0;
                    }
                }
            }
        ++bigCounter;
        }
    return bestKey;
    }
