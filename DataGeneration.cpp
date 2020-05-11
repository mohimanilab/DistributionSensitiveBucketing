/*********************************
 * Filename: DataGeneration.cpp
 * Name: Chengze Shen
 * Date: Oct 28 2019
 * Description: This is a file to generate data files for X and Y strings,
 * which have some preset of insertion, deletion and mutation rates. Each
 * one in X is designed to map to corresponding one at the same index in Y.
 *********************************/
#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <math.h>
#include <random>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>

using namespace std;
typedef chrono::high_resolution_clock Clock;

vector<char> alphabet{'A', 'C', 'G', 'T'};

ostream &operator<<(ostream &os, const vector<char> &seq) {
    for (int i = 0; i < seq.size(); i++)
        os << seq[i];
    return os;
}

void print_helper() {
    printf("Usage: ./DataGeneration -i [insertion rate] -d [deletion rate] -m "
           "[mutation rate] -n [number of sequences] -s [initial length of a "
           "sequence] {-f [num]} {-p [path]} {-vh}\n\n");
    printf("-h\t\tPrint this block of information.\n");
    printf("-v\t\tVerbose mode.\n");
    printf("-i [0<=i<1]\tInsertion rate.\n");
    printf("-d [0<=d<1]\tDeletion rate.\n");
    printf("-m [0<=e<0.5]\tMutation rate when neither insertion/deletion "
           "happens.\n");
    printf("-n [n>0]\tNumber of sequence pairs to generate.\n");
    printf("-s [s>0]\tInitial length of a sequence pair.\n");
    printf("-f [num]\tPrecision for rates that will appear on the output file. "
           "Default: 2\n");
    printf("-p [path]\tWhere generated data will be written to. Default: "
           "data/\n\n");
}

int get_index(int rand_one) {
    int index = -1;
    if (rand_one <= 25)
        index = 0;
    else if (rand_one > 25 && rand_one <= 50)
        index = 1;
    else if (rand_one > 50 && rand_one <= 75)
        index = 2;
    else
        index = 3;
    return index;
}

tuple<char, char> rand_char(default_random_engine &rng, bool mismatch) {
    char x, y;
    // mt19937 rng(time(NULL));
    uniform_int_distribution<int> gen(1, 100);
    int rand_one = gen(rng);
    int index    = -1;

    index = get_index(rand_one);
    x     = alphabet[index];

    if (!mismatch)
        y = alphabet[index];
    else {
        int y_index = index;
        while (y_index == index) {
            y_index = get_index(gen(rng));
        }
        y = alphabet[y_index];
    }

    return make_tuple(x, y);
}

// function to generate data
void gen_data(double insertion, double deletion, double mutation, int num,
              int state, vector<vector<char>> &X, vector<vector<char>> &Y) {
    int i_range = int(100 * insertion);
    int d_range = int(100 * (insertion + deletion));
    int m_range = int(100 * mutation);
    // cout << i_range << ", " << d_range << ", " << m_range << endl;

    // iterate [num] times to create
    for (int i = 0; i < num; i++) {
        int seed = 1234 * i + 5673;
        vector<char> x;
        vector<char> y;
        // mt19937 rng(time(NULL));
        default_random_engine rng(seed);
        uniform_int_distribution<int> gen(1, 100);

        for (int j = 0; j < state; j++) {
            int rand_one = gen(rng);
            // cout << rand_one << endl;
            bool mismatch = 0;
            int choice    = 0;

            // insertion
            if (rand_one <= i_range)
                choice = 0;
            // deletion
            else if (rand_one > i_range && rand_one <= d_range)
                choice = 1;
            // match/mismatch
            else {
                choice   = 2;
                mismatch = (gen(rng) <= m_range);
            }
            // cout << choice;
            auto tmp = rand_char(rng, mismatch);
            if (choice == 0)
                y.push_back(get<1>(tmp));
            else if (choice == 1)
                x.push_back(get<0>(tmp));
            else {
                x.push_back(get<0>(tmp));
                y.push_back(get<1>(tmp));
            }
        }
        // cout << endl;
        X.push_back(x);
        Y.push_back(y);
    }
}

// function to write to local
void write_to(vector<vector<char>> &data, double insertion, double deletion,
              int num_of_seq, int len_of_seq, int which, string path,
              int precision) {
    ofstream file;
    ofstream file2;
    string filename;
    string filename2;
    ostringstream oss;
    ostringstream oss2;

    auto cstart = Clock::now();
    if (which == 0) {
        cerr << "Writing X data files... ";
        oss << path << "s_" << to_string(num_of_seq) << "_"
            << to_string(len_of_seq) << "_" << fixed << setprecision(precision)
            << insertion << "_" << fixed << setprecision(precision) << deletion
            << ".fasta";
        oss2 << path << "X" << to_string(num_of_seq) << "_"
             << to_string(len_of_seq) << "_" << fixed << setprecision(precision)
             << insertion << "_" << fixed << setprecision(precision) << deletion
             << ".txt";
    } else {
        cerr << "Writing Y data files... ";
        oss << path << "q_" << to_string(num_of_seq) << "_"
            << to_string(len_of_seq) << "_" << fixed << setprecision(precision)
            << insertion << "_" << fixed << setprecision(precision) << deletion
            << ".fasta";
        oss2 << path << "Y" << to_string(num_of_seq) << "_"
             << to_string(len_of_seq) << "_" << fixed << setprecision(precision)
             << insertion << "_" << fixed << setprecision(precision) << deletion
             << ".txt";
    }
    filename  = oss.str();
    filename2 = oss2.str();
    file.open(filename);
    file2.open(filename2);

    // write as .fasta/.txt files
    for (int i = 0; i < data.size(); i++) {
        file << ">" << to_string(i) << endl << data[i] << endl;
        file2 << data[i] << endl;
    }
    file.close();
    file2.close();

    auto cend = Clock::now();
    int dur =
        chrono::duration_cast<chrono::milliseconds>(cend - cstart).count();
    cerr << dur << " ms" << endl;
}

/* main function */
int main(int argc, char **argv) {
    double insertion, deletion, mutation;
    int num_of_seq, len_of_seq;
    string path = "data/";

    char opt;
    int verbose    = 0;
    int precision  = 2;
    int error_flag = !(argc >= 11 && argc <= 17);

    while ((opt = getopt(argc, argv, "vhi:d:m:s:n:p:f:")) != -1) {
        switch (opt) {
        case 'h': // print helper text
            print_helper();
            exit(0);
        case 'v': // verbose mode
            verbose = 1;
            break;
        case 'i': // insertion rate
            insertion = atof(optarg);
            break;
        case 'd': // deletion rate
            deletion = atof(optarg);
            break;
        case 'm': // mutation rate
            mutation = atof(optarg);
            break;
        case 's': // initial length of the generated sequence
            len_of_seq = atoi(optarg);
            break;
        case 'n': // total number of sequence pairs (X, Y) to generate
            num_of_seq = atoi(optarg);
            break;
        case 'p': // path to write data to (directory)
            path = optarg;
            path = path + "/";
            break;
        case 'f': // precision (num of digits after decimal) for the rates
                  // that will appear on output files
            precision = atoi(optarg);
            break;
        default:
            error_flag = 1;
            break;
        }
    }

    if (error_flag || insertion < 0 || deletion < 0 || mutation < 0 ||
        (insertion + deletion + mutation) > 1) {
        print_helper();
        exit(error_flag);
    }

    if (verbose) {
        printf("Number of sequences: %d\tSequence length: %d\nInsertion: "
               "%.4f\tDeletion: %.4f\t\tMutation: %.4f\n\n",
               num_of_seq, len_of_seq, insertion, deletion, mutation);
    }

    // generate X,Y sequence pairs with given ins/del/mut rates
    cerr << "Generating X,Y strings... ";
    vector<vector<char>> X;
    vector<vector<char>> Y;
    auto cstart = Clock::now();
    gen_data(insertion, deletion, mutation, num_of_seq, len_of_seq, X, Y);
    auto cend = Clock::now();
    int t0 = chrono::duration_cast<chrono::milliseconds>(cend - cstart).count();
    cerr << t0 << " ms" << endl;

    // write generated sequence pairs to the designated directory
    write_to(X, insertion, deletion, num_of_seq, len_of_seq, 0, path,
             precision);
    write_to(Y, insertion, deletion, num_of_seq, len_of_seq, 1, path,
             precision);
    return 0;
}
