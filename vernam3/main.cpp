#include "header.h"

int main() {
//    ios::sync_with_stdio(false);
//    cin.tie(nullptr);
    auto start = chrono::high_resolution_clock::now();

    locale::global(locale("vi_VN.utf8"));

    try {
////        sort(symbols.begin(), symbols.end());
////        toVieAlpb("plain.txt", "plain_viet.txt");
        freq_analysis("plain_viet.txt", "analysis_plain.txt", Count, Sum);
//        gen_vernam_key("key", Sum, 7, 7, 7, 7, 5, 7);
        gen_vernam_key_rand("key", Sum, 3);
        vernam_crypt(mpCharBit, mpBitChar, "plain_viet.txt", "key", "code.txt", true);
        freq_analysis("code.txt", "analysis_code.txt", Count= 0, Sum= 0);
        vernam_crypt(mpCharBit, mpBitChar, "code.txt", "key", "decode.txt", false);
        freq_analysis_key("key", "analysis_key.txt");
        codeAnalysis2("code.txt", "analysis_code.txt", Count= 0, Sum= 0);

//        multimap<ull, pair<wchar_t, wchar_t>, decltype(&compare)> map_bigrams(&compare);
//        countBigrams(map_bigrams, "plain_viet.txt", "bigrams.txt");
//        decrypt("code.txt", "outFile.txt", climbMatrix);

//        wchar_t ciphertext[MAXTEXTLEN] = L"ởp7ọ}mộặổ ệ;hẻ.ỵ{aẳếvo&càị ẩgưềhỡạtoỳìxr";
//        srand(time(0));
//        hill_climbing_attack(ciphertext);

        }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        }
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = finish - start;
    cout << "\nProgram execution time: " << duration.count() << " seconds" << endl;

    return 0;
    }
