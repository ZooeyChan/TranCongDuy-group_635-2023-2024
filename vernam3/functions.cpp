#include "header.h"
#include "monograms.h"
#include "tetragrams.h"

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

void gen_vernam_key(const string& fileOutput, ull Count, int pos1, int pos2, int pos3, int pos4, int pos5, int pos6) {
    ofstream out(fileOutput, ios::binary | ios::trunc);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint8_t> distribution(0, CHARS);
    uint8_t ch;
    while (Count--) {
        ch = distribution(gen);
        // đặt bit tại vị trí position bằng 0(position ignored in gamma)
        ch &= ~(1 << pos1);
        ch &= ~(1 << pos2);
        ch &= ~(1 << pos3);
        ch &= ~(1 << pos4);
        ch &= ~(1 << pos5);
        ch &= ~(1 << pos6);
        out.write(reinterpret_cast<const char*>(&ch), sizeof(uint8_t));
        }
    out.close();
    }

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
    cout << endl;
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
            ++index;
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
    assignMap(2, tmp);
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
    double b, c, d, e, i, j, m, f;
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

size_t lam_tron(size_t number, const int &k) {
    size_t factor = static_cast<size_t>(pow(10, k));
    size_t rounded_number = (number / factor) * factor;
    return rounded_number;
    }

void codeAnalysis2(const string& fileInput, const string& fileOutput, ull &Count, ull &Sum) {
    map<wchar_t, ull> mp;
    set<size_t> uniqueNumbers;
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
                if (c != L' ') Count++;
                }
            }
        }
    short k = 3;
    short dem1 = 1;
    size_t tmp = lam_tron(mp.begin()->second, k);
    for (const auto &c : mp) {
        size_t current = lam_tron(c.second, k);
        uniqueNumbers.insert(current);
        cout << current << endl;
        if (current != tmp) {
            tmp = current;
            dem1 ++;
            }
        }
    if(log2(pow(2, BIT) / dem1) < 0.1) {
        outputFile << "\nThe first missing bit in the Vernam key is at position: " << '0';
        }
    else outputFile << "\nThe first missing bit in the Vernam key is at position: " << ceil(log2(pow(2, BIT) / dem1));
    outputFile << "\nTotal missing bit in the Vernam key: " << int(log2(uniqueNumbers.size()));
    inputFile.close();
    outputFile.close();
    }

void findCycleLengths(const string& s) {
    int n = s.size();
    for (int len = 1; len <= n / 2; len *= 2) {
        bool isCycle = true;
        for (int i = 0; i < len; ++i) {
            if (i + len < n && s[i] != s[i + len]) {
                isCycle = false;
                break;
                }
            }
        if (!isCycle) {
            cout << len << " ";
            }
        }
    cout << endl;
    }
//    string s = "1122334411223344556677885566778811223344112233445566778855667788";
//    string s = "1122112233443344112211223344334455555555667766775555555566776677";
//void codeAnalysis3(const string& fileInput, const string& fileOutput, ull &Count, ull &Sum) {
//    map<wchar_t, ull> mp;
//    set<size_t> uniqueNumbers;
//    wifstream inputFile(fileInput);
//    if (!inputFile.is_open()) {
//        cerr << "Input error!!!" << endl;
//        return;
//        }
//    wofstream outputFile(fileOutput, ios::app);
//    if (!outputFile.is_open()) {
//        cerr << "Output error!!!" << endl;
//        inputFile.close();
//        return;
//        }
//    wstring str;
//    while (getline(inputFile, str)) {
//        for (wchar_t c : str) {
//            if (binSearch(symbols, tolower(c))) {
//                mp[tolower(c)]++;
//                Sum++;
//                if (c != L' ') Count++;
//                }
//            }
//        }
//    string s[128] = "";
//    int tmp = lam_tron(mp.begin()->second, k);
//    for (const auto &c : mp) {
//        int current = lam_tron(c.second, k);
//        s += to_string(current);
//        }
//
//    inputFile.close();
//    outputFile.close();
//    }

/*---------------------DEMO-----------------------------*/

wchar_t bitToChar(const bitset<BIT>& b) {
    auto it = mpBitChar.find(b);
    if (it != mpBitChar.end()) {
        return it->second;
        }
    return L'?';
    }

double fitness(const wchar_t* text) {
    int length, dem = 0;
    double result = 0;
    length = wcslen(text);
    if (length < 4) {
        return 0;
        }
    for (int i = 0; i < length - 3; i++) {
        int index = (mpCharBit[text[i + 0]].to_ulong() * 26 * 26 * 26) +
                    (mpCharBit[text[i + 1]].to_ulong() * 26 * 26) +
                    (mpCharBit[text[i + 2]].to_ulong() * 26) +
                    (mpCharBit[text[i + 3]].to_ulong());
        if (index >= 0 && index < 26 * 26 * 26 * 26) {
            result += tetragrams[index];
            dem++;
            }
        }
    if (dem > 0) {
        return result / dem;
        }
    else {
        return 0;
        }
    }

void decrypt(const wchar_t* c, wchar_t* p, const wchar_t* key) {
    int length = wcslen(c);
    for (int i = 0; i < length; i++) {
        bitset<BIT> c_bit = mpCharBit[c[i]];
        bitset<BIT> k_bit = mpCharBit[key[i]];
        bitset<BIT> p_bit = c_bit ^ k_bit;
        if (mpBitChar.find(p_bit) != mpBitChar.end()) {
            p[i] = mpBitChar[p_bit];
            }
        else {
            p[i] = L'?'; // Lỗi
            }
        }
    p[length] = L'\0';
    }

void random_key(wchar_t* key, int length) {
    vector<wchar_t> chars;
    for (const auto& pair : mpCharBit) {
        chars.push_back(pair.first);
        }
    for (int i = 0; i < length; i++) {
        key[i] = chars[rand() % chars.size()];
        }
    key[length] = L'\0'; // chuỗi kết thúc bằng ký tự null
    }

void copy_key(const wchar_t* source, wchar_t* target, int length) {
    for (int i = 0; i < length; i++) {
        target[i] = source[i];
        }
    target[length] = L'\0';
    }

void hill_climbing_attack(const wchar_t* ciphertext) {
    wchar_t parent_key[MAXTEXTLEN];
    wchar_t child_key[MAXTEXTLEN];
    wchar_t plaintext[MAXTEXTLEN];
    wchar_t best_plaintext[MAXTEXTLEN];
    wchar_t best_key[MAXTEXTLEN];
    double parent_fitness, child_fitness, best_fitness;
    long int dem, bigcount = 0;
    int length = wcslen(ciphertext);
    random_key(parent_key, length);
    decrypt(ciphertext, plaintext, parent_key);
    best_fitness = fitness(plaintext);
    while (bigcount < 1000000 * 40) {
        random_key(parent_key, length);
        decrypt(ciphertext, plaintext, parent_key);
        parent_fitness = fitness(plaintext);
        dem = 0;
        while (dem < 1000) {
            copy_key(parent_key, child_key, length);
            int pos = rand() % length;
            auto it = mpCharBit.begin();
            advance(it, rand() % mpCharBit.size());
            child_key[pos] = it->first;
            decrypt(ciphertext, plaintext, child_key);
            child_fitness = fitness(plaintext);
            if (child_fitness > parent_fitness) {
                copy_key(child_key, parent_key, length);
                parent_fitness = child_fitness;
                dem = 0;
                }
            else {
                dem++;
                }
            if (child_fitness > best_fitness) {
                copy_key(child_key, best_key, length);
                best_fitness = child_fitness;
                bigcount = 0;
                wcscpy(best_plaintext, plaintext);
                }
            else {
                bigcount++;
                }
            }
        }
    wcout << best_plaintext << endl;
    wcout << L"Key: " << best_key << endl;
    wcout << L"Fitness: " << best_fitness << endl;
    }
