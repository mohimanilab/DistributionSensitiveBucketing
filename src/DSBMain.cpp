/*********************************
 * Filename: DSBMain.cpp
 * Name: Chengze Shen
 * Date: Mar 3rd, 2020
 * Description: A decision-like graph that maps data points to buckets.
 * Data points mapped to the same bucket will be considered to be candidates
 * of having overlaps, and then will be further filtered by a quick alignment
 * step. The filtered overlaps will be reported.
 *********************************/
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <math.h>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
typedef chrono::high_resolution_clock Clock;

/**************** Classes/Struct definitions ****************/
// global alignment results
struct alnresult {
    vector<char> aln_x;    // alignment of x
    vector<char> aln_y;    // alignment of y; should have same size as aln_x
    double identity = 0.0; // identity of the alignment
    int start_idx   = 0;   // starting idx in aln_x/aln_y that are not both gaps
    int aln_score   = 0; // alignment score based on match/mismatch/indel scores

    alnresult(){};
};

// simple struct for kmer tree of X/Y
struct kmernode {
    vector<char> kmer;
    map<int, vector<int>> hashed;
    bool visited = 0;

    kmernode(){};

    int size() { return hashed.size(); }
};

// struct for node in this tree (not sure what to name)
struct treenode {
    double p, q, px, py;
    vector<char> x;
    vector<char> y;
    kmernode *hashX;
    kmernode *hashY;

    treenode(double p, double px, double py, vector<char> x, vector<char> y) {
        this->p  = p;
        this->px = px;
        this->py = py;
        this->q  = px * py;
        this->x  = x;
        this->y  = y;
    }

    treenode(double p, double px, double py, vector<char> x, vector<char> y,
             kmernode *hashX, kmernode *hashY) {
        this->p     = p;
        this->px    = px;
        this->py    = py;
        this->q     = px * py;
        this->x     = x;
        this->y     = y;
        this->hashX = hashX;
        this->hashY = hashY;
    }

    treenode(treenode *copy) {
        this->p     = copy->p;
        this->px    = copy->px;
        this->py    = copy->py;
        this->q     = copy->q;
        this->x     = copy->x;
        this->y     = copy->y;
        this->hashX = copy->hashX;
        this->hashY = copy->hashY;
    }

    // flip vector x and vector y
    void flipxy() {
        vector<char> temp(this->x);
        this->x = this->y;
        this->y = temp;
    }

    // if want to compare y, flip x and y first then redo sorting
    bool operator<(const treenode &other) const {
        int i = 0;
        while (i < x.size() && i < other.x.size()) {
            if (x[i] < other.x[i])
                return 1;
            else if (x[i] > other.x[i])
                return 0;
            i++;
        }
        return (x.size() < other.x.size());
    }
};
typedef struct treenode trn;

// object comparator for vector of char (X, Y)
struct less_than {
    inline bool operator()(const vector<char> &s1, const vector<char> &s2) {
        int i = 0;
        while (i < s1.size() && i < s2.size()) {
            if (s1[i] < s2[i])
                return 1;
            else if (s1[i] > s2[i])
                return 0;
            i++;
        }
        return (s1.size() < s2.size());
    }
} less_than;

// object comparator for vector of char (X, Y)
struct greater_than {
    inline bool operator()(const vector<char> &s1, const vector<char> &s2) {
        int i = 0;
        while (i < s1.size() && i < s2.size()) {
            if (s1[i] > s2[i])
                return 1;
            else if (s1[i] < s2[i])
                return 0;
            i++;
        }
        return (s1.size() > s2.size());
    }
} greater_than;

// sort function for vector - sort by second element in pair
bool sort_by_second(const pair<int, int> &a, const pair<int, int> &b) {
    return (a.second < b.second);
}

// sort function for vector - sort by second element then by first
bool sort_by_second_then_first(const pair<int, int> &a,
                               const pair<int, int> &b) {
    if (a.second != b.second)
        return (a.second < b.second);
    else
        return (a.first < b.first);
}

/**************** End of classes/struct definitions ****************/

/**************** operator overloading for print ****************/

ostream &operator<<(ostream &os, const treenode &node) {
    for (int i = 0; i < node.x.size(); i++)
        os << node.x[i];
    os << endl;
    for (int i = 0; i < node.y.size(); i++)
        os << node.y[i];
    os << endl;
    return os;
}

ostream &operator<<(ostream &os, const vector<char> &seq) {
    for (int i = 0; i < seq.size(); i++)
        os << seq[i];
    return os;
}

ostream &operator<<(ostream &os, const vector<int> &seq) {
    for (int i = 0; i < seq.size(); i++) {
        if (i != seq.size() - 1)
            os << seq[i] << ",";
        else
            os << seq[i];
    }
    return os;
}

ostream &operator<<(ostream &os, const map<int, bool> &matches) {
    for (auto it = matches.cbegin(); it != matches.cend(); ++it)
        os << it->first << ":" << it->second << ",";
    return os;
}

ostream &operator<<(ostream &os, const pair<int, int> &pairs) {
    os << pairs.first << "," << pairs.second;
    return os;
}

ostream &operator<<(ostream &os, const set<pair<int, int>> &aset) {
    for (auto it = aset.begin(); it != aset.end(); ++it) {
        os << *it << ";";
    }
    return os;
}

ostream &operator<<(ostream &os, const vector<pair<int, int>> &avec) {
    for (auto it = avec.begin(); it != avec.end(); ++it) {
        os << *it << ";";
    }
    return os;
}

/*************** End of operator overloading ***************/

/***************** functions/variables ********************/

// The 2 sets of vectors that need to be compared
vector<vector<char>> X;
vector<vector<char>> Y;

// alphabet
vector<char> alphabet{'A', 'C', 'G', 'T'};

// Number of vectors in each set
static int N;
static int M;

// independent hash tree mapping of buckets
// map<pair<vector<char>, vector<char>>, trn *> buckets;
vector<trn *> buckets;

static double eps     = -1.;
static double p_ins   = -1.;
static double p_del   = -1.;
static double p_match = -1.;

// penalties for alignment
static const int s_indel    = -1;
static const int s_match    = 1;
static const int s_mismatch = -1;

// method to read in data from given files X and Y
// we are using real data for this one, so it is modified
// FOR {A,T,C,G} DNA alignment
void readData(string xfile, string yfile) {
    string line;
    ifstream data_x;
    data_x.open(xfile);
    ifstream data_y;
    data_y.open(yfile);
    vector<char> fx;
    vector<char> fy;

    if (data_x.is_open()) {
        while (getline(data_x, line)) {
            stringstream ss(line);
            char i;
            while (ss >> i) {
                if (i == '>') {
                    // if fx size is not 0, write fx to X
                    if (fx.size() > 0)
                        X.push_back(fx);
                    fx.clear();
                    goto skipline_x;
                }
                if (ss.peek() == ',' || ss.peek() == ' ' || ss.peek() == '\n')
                    ss.ignore();
                fx.push_back(i);
            }
        skipline_x:
            continue;
        }
    } else {
        cerr << "X file not found in " << xfile << endl;
        exit(1);
    }
    if (fx.size() > 0)
        X.push_back(fx);
    data_x.close();

    if (data_y.is_open()) {
        while (getline(data_y, line)) {
            stringstream ss(line);
            char i;
            while (ss >> i) {
                if (i == '>') {
                    if (fy.size() > 0)
                        Y.push_back(fy);
                    fy.clear();
                    goto skipline_y;
                }
                if (ss.peek() == ',' || ss.peek() == ' ' || ss.peek() == '\n')
                    ss.ignore();
                fy.push_back(i);
            }
        skipline_y:
            continue;
        }
    } else {
        cerr << "Y file not found in " << yfile << endl;
        exit(1);
    }
    if (fy.size() > 0)
        Y.push_back(fy);
    data_y.close();
}

// alignment method for comparing two given data points
// tuple<int, int, double, vector<char>*, vector<char>*>
tuple<int, int, double> alignment(vector<char> &x, vector<char> &y, int xstart,
                                  int ystart, int xrange, int yrange,
                                  bool b_backtrace) {
    // alnresult * aln_result = new alnresult();
    //  rows are x, columns are y
    int xend        = 0;
    int yend        = 0;
    double identity = 0.0; // identity score of alignment

    if (xrange != -1 && yrange != -1) {
        xend = min(xstart + xrange + 1, int(x.size() + 1));
        yend = min(ystart + yrange + 1, int(y.size() + 1));
    } else {
        xend = x.size() + 1;
        yend = y.size() + 1;
    }
    vector<vector<int>> matrix(xend - xstart, vector<int>(yend - ystart, 0));

    // initialize first column/row
    for (int i = 1; i < xend - xstart; i++) {
        matrix[i][0] = matrix[i - 1][0] + s_indel;
    }
    for (int j = 1; j < yend - ystart; j++) {
        matrix[0][j] = matrix[0][j - 1] + s_indel;
    }

    // run alignment
    for (int i = 1; i < xend - xstart; i++) {
        for (int j = 1; j < yend - ystart; j++) {
            // diagonal
            int diag = matrix[i - 1][j - 1];
            if (x[i - 1 + xstart] == y[j - 1 + ystart])
                diag += s_match;
            else
                diag += s_mismatch;

            // indels
            int from_left = matrix[i][j - 1] + s_indel;
            int from_top  = matrix[i - 1][j] + s_indel;

            // check which one is the largest
            if (diag > from_left && diag > from_top)
                matrix[i][j] = diag;
            else if (from_left > diag && from_left > from_top)
                matrix[i][j] = from_left;
            else if (from_top > diag && from_top > from_left)
                matrix[i][j] = from_top;
            else {
                if (diag == from_left)
                    matrix[i][j] = diag;
                else if (diag == from_top)
                    matrix[i][j] = diag;
                else
                    matrix[i][j] = from_top;
            }
        }
    }

    // reconstruct the alignment if backtracing
    if (b_backtrace) {
        // maximum length of alignment
        int max_aln_len = xend - xstart + yend - ystart - 2;
        int xpos        = max_aln_len;
        int ypos        = max_aln_len;
        int aln_len     = 0;
        int exact_match = 0;
        int start_idx   = 0;
        int i           = xend - xstart - 1;
        int j           = yend - ystart - 1;
        vector<char> tmp_x(max_aln_len + 1, '-');
        vector<char> tmp_y(max_aln_len + 1, '-');

        // read from the last position of the matrix to backtrace
        while (i != 0 && j != 0) {
            // diagonal - match
            if (x[xstart + i - 1] == y[ystart + j - 1]) {
                exact_match++;
                tmp_x[xpos--] = x[xstart + i - 1];
                tmp_y[ypos--] = y[ystart + j - 1];
                i--;
                j--;
            }
            // diagonal - mismatch
            else if (matrix[i][j] == matrix[i - 1][j - 1] + s_mismatch) {
                tmp_x[xpos--] = x[xstart + i - 1];
                tmp_y[ypos--] = y[ystart + j - 1];
                i--;
                j--;
            }
            // from top
            else if (matrix[i][j] == matrix[i - 1][j] + s_indel) {
                tmp_x[xpos--] = x[xstart + i - 1];
                tmp_y[ypos--] = '-';
                i--;
            }
            // from left
            else if (matrix[i][j] == matrix[i][j - 1] + s_indel) {
                tmp_x[xpos--] = '-';
                tmp_y[ypos--] = y[ystart + j - 1];
                j--;
            }
            aln_len++;
        }
        // insertions/deletions at head
        while (xpos > 0) {
            if (i > 0) {
                tmp_x[xpos--] = x[xstart + i - 1];
                aln_len++;
                i--;
            } else
                tmp_x[xpos--] = '-';
        }
        while (ypos > 0) {
            if (j > 0) {
                tmp_y[ypos--] = y[ystart + j - 1];
                aln_len++;
                j--;
            } else
                tmp_y[ypos--] = '-';
        }
        identity = ((double)(exact_match)) / ((double)aln_len);

        // get the start of the alignment (i.e., identify the first index where
        // at least one of tmp_x/tmp_y has a character)
        for (int k = max_aln_len; k >= 1; k--) {
            if (tmp_x[k] == '-' && tmp_y[k] == '-') {
                start_idx = k + 1;
                break;
            }
        }
        return make_tuple(matrix[xend - xstart - 1][yend - ystart - 1],
                          start_idx, identity);
    } else {
        return make_tuple(matrix[xend - xstart - 1][yend - ystart - 1], 0, 0.0);
    }
}

/**
 * 1.13.2023 - post-process filtering by percentage shared kmers of a pair of
 * mapped query-target
 */
double pct_shared_kmers(vector<char> &x, vector<char> &y, int xstart, int xend,
                        int ystart, int yend) {
    // use an unordered map to record the kmers of each vector
    map<vector<char>, char> kmers;
    map<vector<char>, char> shared_kmers;
    double pct_shared = 0.0;
    int l_kmer        = 6; // 6-mer

    // iterate over x on its respective range
    for (int i = xstart; i < xend - l_kmer + 1; i++) {
        vector<char> x_subvec(&x[i], &x[i + l_kmer]);
        kmers.emplace(x_subvec, 'x');
    }
    // iterate over y on its range and find the kmer in the map;
    // 1. if it exists, shared_kmers++
    // 2. else, w/e
    for (int i = ystart; i < yend - l_kmer + 1; i++) {
        vector<char> y_subvec(&y[i], &y[i + l_kmer]);
        if (kmers.find(y_subvec) != kmers.end())
            shared_kmers.emplace(y_subvec, 'y');
        kmers.emplace(y_subvec, 'y');
    }

    pct_shared = ((double)shared_kmers.size()) / ((double)kmers.size());
    return pct_shared;
}

/**
 *  1.10.2023 - post-process filtering with a global alignment to check
 *  if the mapped query-target region is valid (e.g., alignment with identity
 *  >= some threshold such as 60%)
 */
bool filter_and_write(ofstream &reports, bool b_alignment,
                      double id_threshold, int q, int t,
                      int q_start, int q_end, int t_start, int t_end) {
    // q and t are 1-based, NEED TO SUBTRACT BY 1
    double identity = 0.0;
    int aln_score   = 0;
    int start_idx   = 0;

    // compute an alignment between query q and target t, in the respectively
    // mapped region q[q_start:q_end] and t[t_start:t_end]
    if (b_alignment) {
        auto ret  = alignment(X[q - 1], Y[t - 1], q_start, t_start, q_end - q_start,
                              t_end - t_start, true);
        aln_score = get<0>(ret);
        start_idx = get<1>(ret);
        identity  = get<2>(ret);
        // identity = pct_shared_kmers(X[q - 1], Y[t - 1], q_start, q_end, t_start,
        // t_end);
    } else {
        // compute shared pct shared kmers
        identity = pct_shared_kmers(X[q - 1], Y[t - 1], q_start, q_end,
                t_start, t_end);
    }


    // if the alignment has >= id_threshold identity score, then we report it
    if (identity >= id_threshold) {
        reports << q << " " << t << " " << q_start << " " << q_end << " "
                << t_start << " " << t_end << " "
                << round(100 * identity / 0.01) * 0.01 << endl;
        return true;
        //<< " " << aln_result->aln_x
        //<< " " << aln_result->aln_y << endl;
    }
    return false;
    // delete aln_result;
}

// function to print needed arguments/documentation for this program
void print_helper() {
    printf("Usage: ./DSBMain -q [query file] -r [reference file] -i [insertion "
           "rate] -d "
           "[deletion rate] -m [mutation rate] -a [add threshold] "
           "-k [kill threshold] -o [name] {-vh}\n\n");
    printf("-h\t\tPrint this block of information.\n");
    printf("\nModel:\n");
    printf("-i [0<=i<1]\tInsertion rate.\n");
    printf("-d [0<=d<1]\tDeletion rate.\n");
    printf("-m [0<=e<0.5]\tMutation rate when neither insertion/deletion "
           "happens.\n");
    printf("-a [a>0]\tThreshold for a node to be considered a bucket.\n");
    printf("-k [k>0]\tThreshold for a node to be pruned.\n");
    printf("-b [path]\tIf specified, will use the given buckets file produced "
           "from a previously curated run\n");
    printf("-s [path]\tInclude this option to save buckets from this run to "
           "local [path]. Disabled if '-b' is specified\n");
    printf("\nOther settings:\n");
    printf("-v\t\tVerbose mode.\n");
    printf("-o [name]\tSpecify output file name.\n");
    printf("-q [path/to/q]\tPath to the query file.\n");
    printf("-r [path/to/r]\tPath to the reference file.\n\n");
}

// function to get total fp
long calc_fp_sum(map<pair<int, int>, long> &mismatches) {
    long total = 0;
    for (auto it = mismatches.cbegin(); it != mismatches.cend(); ++it) {
        total += it->second;
    }
    return total;
}

// function to write buckets to local
void write_buckets(ofstream &buckets_file) {
    for (auto it = buckets.begin(); it != buckets.end(); ++it) {
        buckets_file << (*it)->x << "," << (*it)->y << endl;
    }
}

// function to read buckets from local
// remainder -> buckets: map<pair<vector<char>, vector<char>>, trn *> buckets;
void read_buckets_from_file(string buckets_path) {
    string line;
    ifstream data_buckets;
    data_buckets.open(buckets_path);

    if (data_buckets.is_open()) {
        while (getline(data_buckets, line)) {
            stringstream ss(line);
            vector<char> xmer, ymer;
            string kmer;
            int ind = 0;

            while (getline(ss, kmer, ',')) {
                if (ind == 0) {
                    ind++;
                    copy(kmer.begin(), kmer.end(), back_inserter(xmer));
                } else {
                    copy(kmer.begin(), kmer.end(), back_inserter(ymer));
                    break;
                }
            }

            // handle empty vectors
            if (xmer.size() != 0 && ymer.size() != 0) {
                pair<vector<char>, vector<char>> to_insert =
                    make_pair(xmer, ymer);
                trn *to_add   = new trn(0, 0, 0, xmer, ymer);
                to_add->hashX = new kmernode();
                to_add->hashY = new kmernode();
                buckets.push_back(to_add);
            }
        }
    } else {
        cerr << "buckets file not found in " << buckets_path << endl;
    }
}

// function to insert X and Y to buckets
void insert_to_buckets(vector<vector<char>> &data, int t) {
    // initialize to check
    for (int ind = 0; ind < data.size(); ind++) {
        vector<char> tohash = data[ind];
        // examine every position
        for (int i = 0; i < tohash.size(); i++) {
            int front = 0, end = buckets.size(), mid = -1;
            while (front != end) {
                mid      = (front + end) / 2;
                int len  = buckets[mid]->x.size();
                int tail = min(int(tohash.size() - 1), i + len);
                vector<char> subvec(&tohash[i], &tohash[tail]);

                if (subvec == buckets[mid]->x || buckets[mid]->x.size() == 0) {
                    if (t == 0)
                        buckets[mid]->hashX->hashed[ind].push_back(i + len);
                    else
                        buckets[mid]->hashY->hashed[ind].push_back(i + len);

                    int f = mid - 1;
                    int e = mid + 1;

                    // local search for all matches
                    while (f >= front) {
                        int f_len  = buckets[f]->x.size();
                        int f_tail = min(int(tohash.size() - 1), i + f_len);
                        vector<char> temp(&tohash[i], &tohash[f_tail]);

                        if (temp == buckets[f]->x) {
                            if (t == 0)
                                buckets[f]->hashX->hashed[ind].push_back(i +
                                                                         f_len);
                            else
                                buckets[f]->hashY->hashed[ind].push_back(i +
                                                                         f_len);
                            f--;
                        } else
                            break;
                    }
                    while (e < end) {
                        int e_len  = buckets[e]->x.size();
                        int e_tail = min(int(tohash.size() - 1), i + e_len);
                        vector<char> temp(&tohash[i], &tohash[e_tail]);

                        if (temp == buckets[e]->x) {
                            if (t == 0)
                                buckets[e]->hashX->hashed[ind].push_back(i +
                                                                         e_len);
                            else
                                buckets[e]->hashY->hashed[ind].push_back(i +
                                                                         e_len);
                            e++;
                        } else
                            break;
                    }
                    break; // break from the searching for loop
                } else if (subvec > buckets[mid]->x) {
                    front = mid + 1;
                } else {
                    end = mid;
                }
            }
        }
    }
}

// method to pre-allocate a tree for X or Y-mer so that later we only
// need to refer to these nodes
void make_kmer_tree(vector<vector<char>> &target, int ind, int b,
                    map<vector<char>, kmernode *> &mapping) {
    // root, init kmer
    kmernode *root = new kmernode();
    vector<char> init_kmer;
    root->kmer = init_kmer;

    // root, init hashed
    map<int, vector<int>> init_hashed;
    for (int i = 0; i < target.size(); i++) {
        /* to insert every position into the kmer tree */
        for (int j = 0; j < target[i].size(); j++) {
            init_hashed[i].push_back(j);
        }
        /* new method end */
    }
    root->hashed = init_hashed;

    // put root to start of the queue
    mapping[init_kmer] = root;
    queue<kmernode *> q;
    q.push(root);

    // early return, only get the initial hash
    return;
}

// function to visit the kmernode and make its children
void visit(map<vector<char>, kmernode *> &prefix_tree, kmernode *cur,
           int branch) {
    cur->visited = 1;

    // 4 possible children
    for (int i = 0; i < 4; i++) {
        vector<char> next_kmer(cur->kmer);
        next_kmer.push_back(alphabet[i]);

        kmernode *next_kmernode = new kmernode();
        next_kmernode->kmer     = next_kmer;

        for (auto it = cur->hashed.cbegin(); it != cur->hashed.cend(); ++it) {
            for (int m = 0; m < it->second.size(); m++) {
                // prefix tree for X
                if (branch == 0) {
                    if (it->second[m] >= X[it->first].size())
                        continue;

                    if (X[it->first][it->second[m]] == next_kmer.back()) {
                        next_kmernode->hashed[it->first].push_back(
                            it->second[m] + 1);
                    }
                }
                // prefix tree for Y
                else {
                    if (it->second[m] >= Y[it->first].size())
                        continue;

                    if (Y[it->first][it->second[m]] == next_kmer.back()) {
                        next_kmernode->hashed[it->first].push_back(
                            it->second[m] + 1);
                    }
                }
            }
        }
        prefix_tree[next_kmer] = next_kmernode;
    }
}

// new tree construction, with no premade prefix tree
tuple<long, long, double, double, double, double>
indpt_tree(double kill_threshold, double add_threshold,
           map<vector<char>, kmernode *> &Xmer_tree,
           map<vector<char>, kmernode *> &Ymer_tree) {
    // cerr << "Debug entering independent tree construction" << endl;
    queue<trn *> q; // use queue
    map<pair<vector<char>, vector<char>>, trn *> all_nodes;

    long total      = 0;
    long num_bucket = 0;

    // initial node
    vector<char> init_x;
    vector<char> init_y;
    kmernode *init_hash_X = Xmer_tree[init_x];
    kmernode *init_hash_Y = Ymer_tree[init_y];

    //cerr << "initial hashed entries: " << \
        init_hash_X->hashed.size() << " " << init_hash_Y->hashed.size() << endl;

    trn *root =
        new trn(1.0, 1.0, 1.0, init_x, init_y, init_hash_X, init_hash_Y);
    q.push(root);
    all_nodes[make_pair(root->x, root->y)] = root;
    double total_p  = 0.0; // total p of all accepted buckets
    double total_px = 0.0; // total p_x
    double total_py = 0.0; // total p_y
    double total_q  = 0.0; // total q

    // keep inserting nodes and proceeding until no more nodes
    int counter = 0;
    while (!q.empty()) {
        counter++;
        // pop the next node
        trn *working = q.front();
        q.pop();
        total++;

        // cerr << working->p << " " << working->q << " " <<
        // working->hashX->size() << " " << working->hashY->size() << endl;
        // check current node terminating conditions
        if ((working->hashX->size() + working->hashY->size()) / working->p >
            kill_threshold)
            continue;
        /* 1.11.2023 - a failsafe for making kmers too deep */
        if (working->x.size() > 9 || working->y.size() > 9)
            continue;
        /* failsafe ends */
        // if (working->px / working->p > kill_threshold)
        if (working->hashX->size() == 0)
            continue; // if no more X prefixes
        if (working->hashY->size() == 0)
            continue; // if no more Y prefixes
        // 2) p/q > threshold_1, add to hash tree
        if ((working->p / working->q) > add_threshold) {
            pair<vector<char>, vector<char>> path =
                make_pair(working->x, working->y);
            trn *bucket = new trn(working);
            buckets.push_back(bucket);
            num_bucket++;

            total_p += bucket->p;
            total_px += bucket->px;
            total_py += bucket->py;
            total_q += bucket->q;
            // cerr << "added a node \n" << *working;
            continue;
        }

        // visit current node's prefix nodes (X and Y)
        // create their children if the prefix exists and not already visited
        // cerr << "reached here" << endl;
        if (working->hashX->size() > 0 && !(working->hashX->visited)) {
            visit(Xmer_tree, working->hashX, 0);
        }
        if (working->hashY->size() > 0 && !(working->hashY->visited)) {
            visit(Ymer_tree, working->hashY, 1);
        }

        vector<char> prev_x = working->x;
        vector<char> prev_y = working->y;
        double p, px, py;

        // construct 8 children, following by the order: ins->del->match
        // 1) insertions
        // cerr << "doing insertion" << endl;
        for (int i = 0; i < 4; i++) {
            vector<char> cur_x = prev_x;
            vector<char> cur_y = prev_y;

            // +1 y, x stays the same
            cur_y.push_back(alphabet[i]);
            p = working->p * p_ins * 0.25;

            // check if the new kmer y exists in the data
            if (Ymer_tree.find(cur_y) == Ymer_tree.end())
                continue;

            // check if there is already an existing node of same (x,y)
            // if so, update that node (no need to calculate px,py again)
            if (all_nodes.count(make_pair(cur_x, cur_y)) > 0) {
                // cerr << "updating node \n" <<
                // *all_nodes[make_pair(cur_x,cur_y)];
                // since we stored pointers, the corresponding existing
                // pointers in queue will also point to the modified node
                all_nodes[make_pair(cur_x, cur_y)]->p += p;
            } else {
                // set px, py = randomly matched with given length
                px = pow(0.25, cur_x.size());
                py = pow(0.25, cur_y.size());

                trn *insertion   = new trn(p, px, py, cur_x, cur_y);
                insertion->hashX = Xmer_tree[cur_x];
                insertion->hashY = Ymer_tree[cur_y];

                q.push(insertion);
                // next_nodes[0].push_back(insertion);
                all_nodes[make_pair(insertion->x, insertion->y)] = insertion;
            }
        }

        // cerr << "doing deletion" << endl;
        // 2) deletions
        for (int i = 0; i < 4; i++) {
            vector<char> cur_x = prev_x;
            vector<char> cur_y = prev_y;

            // +1 x, y stays the same
            cur_x.push_back(alphabet[i]);
            p = working->p * p_del * 0.25;

            // check if the new kmer x exists in the data
            if (Xmer_tree.find(cur_x) == Xmer_tree.end())
                continue;

            // check if there is already an existing node of same (x,y)
            // if so, update that node (no need to calculate px,py again)
            if (all_nodes.count(make_pair(cur_x, cur_y)) > 0) {
                // since we stored pointers, the corresponding existing
                // pointers in queue will also point to the modified node
                all_nodes[make_pair(cur_x, cur_y)]->p += p;
            } else {
                // set px, py = randomly matched with given length
                px = pow(0.25, cur_x.size());
                py = pow(0.25, cur_y.size());

                trn *deletion   = new trn(p, px, py, cur_x, cur_y);
                deletion->hashX = Xmer_tree[cur_x];
                deletion->hashY = Ymer_tree[cur_y];

                q.push(deletion);
                // next_nodes[1].push_back(deletion);
                all_nodes[make_pair(deletion->x, deletion->y)] = deletion;
            }
        }

        // cerr << "doing matches" << endl;
        // 3) matches
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                vector<char> cur_x = prev_x;
                vector<char> cur_y = prev_y;
                cur_x.push_back(alphabet[i]);
                cur_y.push_back(alphabet[j]);

                if (i == j)
                    p = working->p * p_match * (0.24 - eps);
                else
                    p = working->p * p_match * eps;

                // check if both of the new kmer exists, if not don't make
                if (Xmer_tree.find(cur_x) == Xmer_tree.end() ||
                    Ymer_tree.find(cur_y) == Ymer_tree.end())
                    continue;

                if (all_nodes.count(make_pair(cur_x, cur_y)) > 0) {
                    all_nodes[make_pair(cur_x, cur_y)]->p += p;
                } else {
                    // since this insertion will never has a case that cur_x
                    // or cur_y be empty, px/py will have all 3 incoming edges
                    // for calculation
                    pair<vector<char>, vector<char>> from_del =
                        make_pair(prev_x, cur_y);
                    pair<vector<char>, vector<char>> from_ins =
                        make_pair(cur_x, prev_y);
                    pair<vector<char>, vector<char>> from_match =
                        make_pair(prev_x, prev_y);

                    // set px, py = randomly matched with given length
                    px = pow(0.25, cur_x.size());
                    py = pow(0.25, cur_y.size());

                    trn *matched   = new trn(p, px, py, cur_x, cur_y);
                    matched->hashX = Xmer_tree[cur_x];
                    matched->hashY = Ymer_tree[cur_y];

                    q.push(matched);
                    // next_nodes[2].push_back(matched);
                    all_nodes[make_pair(matched->x, matched->y)] = matched;
                }
            }
        }
    }
    //cerr << num_bucket << "," << total << "," << total_p << "," << total_px \
    << "," << total_py << "," << total_q << endl;
    return make_tuple(num_bucket, total, total_p, total_px, total_py, total_q);
}

/***************** End of functions/variables ********************/

int main(int argc, char **argv) {
    string xfile             = "";
    string yfile             = "";
    string output_name       = "output.txt";
    string buckets_path      = "";
    string save_buckets_path = "";
    double id_threshold      = 0.7;
    double kmer_threshold    = 0.25;
    double add_threshold, kill_threshold;
    double align_multi = 0.3; // for pre-filtering alignment
    int range          = 9;   // pre-filtering alignment length
    int search_range   = 300; // maximum gap for connecting two maps
    int map_len_thres  = 200; // minimum alignment length to report
    int long_len_thres  = 1000; // alignment above this will be checked with
                                // % shared 6-kmers
    int verbose        = 0;
    bool save_buckets  = false; // whether to save buckets info for current run
    bool read_buckets  = false; // whether to use preexisting buckets

    // argument parsing
    int error_flag = 0;
    // int error_flag =
    //    !(argc == 15 || argc == 16 || argc == 17 || argc == 18 || argc == 19);
    char opt;
    while ((opt = getopt(argc, argv, "hvr:q:i:d:m:a:k:o:b:s:")) != -1) {
        switch (opt) {
        case 'q': // X datafile address
            xfile = optarg;
            break;
        case 'r': // Y datafile address
            yfile = optarg;
            break;
        case 'i': // insertion probability
            p_ins = atof(optarg);
            break;
        case 'd': // deletion probability
            p_del = atof(optarg);
            break;
        case 'm': // epsilon (mismatch prob while neither ins/del happens
            eps = atof(optarg);
            eps /= 12;
            break;
        case 'a': // threshold 1 - for a node to be considered a bucket
            add_threshold = atof(optarg);
            break;
        case 'k': // threshold 2 - for a node to be pruned
            kill_threshold = atof(optarg);
            break;
        case 'v': // verbose mode
            verbose = 1;
            break;
        case 'b': // read in buckets - from the given path
            buckets_path = optarg;
            read_buckets = true;
            break;
        case 's': // save buckets to local
            save_buckets_path = optarg;
            save_buckets      = true;
            break;
        case 'o': // specify output file name
            output_name = optarg;
            break;
        case 'h': // print all needed arguments
            print_helper();
            exit(0);
        default:
            error_flag = 1;
            break;
        }
    }

    // check if any input file is missing
    if (xfile.compare("") == 0 || yfile.compare("") == 0) {
        cerr << "Either query or target file are missing." << endl;
        print_helper();
        exit(1);
    }

    // check error flag for arguments
    if (error_flag || add_threshold <= 0 || kill_threshold <= 0) {
        cerr << "Thresholds need to be positive." << endl;
        print_helper();
        exit(1);
    }
    // check if insertion deletion and mutation has been set
    if (p_ins == -1. || p_del == -1. || eps == -1) {
        cerr << "Need to set all paramters: ins/del/mutation rates" << endl;
        print_helper();
        exit(1);
    }
    // check if buckets paths are not empty
    if (save_buckets && save_buckets_path.size() == 0) {
        cerr << "Save buckets path is empty!" << endl;
        print_helper();
        exit(1);
    }
    if (read_buckets && buckets_path.size() == 0) {
        cerr << "Read buckets path is empty!" << endl;
        print_helper();
        exit(1);
    }

    p_match = 1.0 - p_ins - p_del;

    if (verbose) {
        printf("X: %s\tY: %s\nOutput: %s\n\tUse buckets: %d\tPath: %s\n\tSave "
               "buckets: %d\tPath: %s\nInsertion: %.3f\tDeletion: "
               "%.3f\t\tMutation: "
               "%.3f\nThreshold 1: %.1f\tThreshold 2: %.1f\n\n",
               xfile.c_str(), yfile.c_str(), output_name.c_str(), read_buckets,
               buckets_path.c_str(), save_buckets, save_buckets_path.c_str(),
               p_ins, p_del, eps * 12, add_threshold, kill_threshold);
    }

    // call function to read in data
    readData(xfile, yfile);

    M = X.size();
    N = Y.size();

    // sort X and Y from A -> T
    auto cstart = Clock::now();
    // sort(X.begin(), X.end(), less_than);
    // sort(Y.begin(), Y.end(), less_than);
    auto cend = Clock::now();
    // int t_preprocess =
    // chrono::duration_cast<chrono::milliseconds>(cend-cstart).count();

    // pre-allocate two tree structures for X and Y, nodes are kmers which
    // contain the strings that have the kmers
    // run for each b, sum all up
    // unordered_map<int, bool> matches;
    // map<pair<int, int>, long> mismatches;
    map<pair<int, int>, vector<pair<int, int>>> results;
    // map<pair<int, int>, pair<int, int>> cur_mapped_low;
    ofstream reports, buckets_file;
    int t0, t1, t2, t3 = 0;
    double xmer, ymer            = 0.0;
    long num_bucket, total_nodes = 0;
    double total_p, total_px, total_py, total_q = 0.0;

    /* LOCAL ALIGNMENT enabled. Check all positions in X matching to positions
     * in Y
     */
    map<vector<char>, kmernode *> Xmer_tree;
    map<vector<char>, kmernode *> Ymer_tree;

    // (1) make kmer trees for X/Y with every position inserted
    cerr << "Making k-mer trees for query and reference... ";
    cstart = Clock::now();
    make_kmer_tree(X, 0, 0, Xmer_tree);
    make_kmer_tree(Y, 0, 0, Ymer_tree);
    cend = Clock::now();
    t0   = chrono::duration_cast<chrono::milliseconds>(cend - cstart).count();
    cerr << t0 << " ms" << endl;

    // (2) start constructing decision tree
    if (!read_buckets) {
        cerr << "Making data-dependent tree for buckets... ";
        cstart = Clock::now();
        auto tmp =
            indpt_tree(kill_threshold, add_threshold, Xmer_tree, Ymer_tree);
        sort(buckets.begin(), buckets.end(), less<trn>());
        cend = Clock::now();
        t1 = chrono::duration_cast<chrono::milliseconds>(cend - cstart).count();
        cerr << t1 << " ms" << endl;

        num_bucket  = get<0>(tmp);
        total_nodes = get<1>(tmp);
        total_p     = get<2>(tmp);
        total_px    = get<3>(tmp);
        total_py    = get<4>(tmp);
        total_q     = get<5>(tmp);
    } else {
        // or use existing buckets and insert X,Y to the buckets
        cerr << "Using existing buckets file... ";
        cstart = Clock::now();
        read_buckets_from_file(buckets_path);

        // sort buckets
        sort(buckets.begin(), buckets.end(), less<trn>());

        // insert X/Y to buckets
        cerr << "Inserting X... ";
        insert_to_buckets(X, 0);
        for (int i = 0; i < buckets.size(); i++) {
            buckets[i]->flipxy();
        }
        cerr << "Inserting Y... ";
        insert_to_buckets(Y, 1);
        for (int i = 0; i < buckets.size(); i++) {
            buckets[i]->flipxy();
        }

        // quickly sort the vector just in case
        for (auto it = buckets.begin(); it != buckets.end(); ++it) {
            for (auto xit = (*it)->hashX->hashed.begin();
                 xit != (*it)->hashX->hashed.end(); ++xit) {
                sort(xit->second.begin(), xit->second.end());
            }
            for (auto yit = (*it)->hashY->hashed.begin();
                 yit != (*it)->hashY->hashed.end(); ++yit) {
                sort(yit->second.begin(), yit->second.end());
            }
        }

        cend = Clock::now();
        t1 = chrono::duration_cast<chrono::milliseconds>(cend - cstart).count();
        cerr << t1 << " ms" << endl;
    }

    // (DEBUG) iterate through all buckets and output corresponding
    // X/Y out with correpsonding indexes
    // for (auto it = buckets.begin(); it != buckets.end(); ++it) {
    //    for (auto xit = (*it)->hashX->hashed.begin();
    //            xit != (*it)->hashX->hashed.end(); ++xit) {
    //        sort(xit->second.begin(), xit->second.end());
    //        // xout first
    //        cout << xit->first << ":";
    //        for (int i = 0; i < xit->second.size(); i++) {
    //            cout << xit->second[i] << ",";
    //        }
    //        cout << ";";
    //    }
    //    cout << endl;
    //    for (auto yit = (*it)->hashY->hashed.begin();
    //            yit != (*it)->hashY->hashed.end(); ++yit) {
    //        sort(yit->second.begin(), yit->second.end());
    //        //yout second
    //        cout << yit->first << ":";
    //        for (int i = 0; i < yit->second.size(); i++) {
    //            cout << yit->second[i] << ",";
    //        }
    //        cout << ";";
    //    }
    //    cout << endl;
    //}
    // exit(0);

    // (2.a) OPTIONAL: save bucket information to local
    if (save_buckets && !read_buckets) {
        cout << "Save buckets option detected, saving buckets to "
             << save_buckets_path << endl;
        buckets_file.open(save_buckets_path, ofstream::out);
        write_buckets(buckets_file);
        buckets_file.close();
        cout << "Exiting..." << endl;
        exit(0);
    }

    // (3) collect results from all buckets
    int align_thres = int(align_multi * range);
    cerr << "Collecting results (bucket size: " << buckets.size() << ")... ";

    // calculate memory usage from buckets: this could be the "actual"
    // memory usage if we separate tree construction from assigning
    // data into the buckets
    // long bucket_memory = 0;
    // for (auto it = buckets.cbegin(); it != buckets.cend(); ++it) {
    //    bucket_memory += sizeof(it->first);
    //    bucket_memory += sizeof(it->second->x);
    //    bucket_memory += sizeof(it->second->y);
    //    bucket_memory += 4 * sizeof(double) + 16;
    //    for (auto xit = it->second->hashX->hashed.cbegin();
    //            xit != it->second->hashX->hashed.cend(); ++xit) {
    //        bucket_memory += xit->second.size() * 4 + 4;
    //    }

    //    for (auto yit = it->second->hashY->hashed.cbegin();
    //             yit != it->second->hashY->hashed.cend(); ++yit) {
    //        bucket_memory += yit->second.size() * 4 + 4;
    //    }
    //}
    // cerr << "Bucket memory: " << bucket_memory << " bytes" << endl;
    // exit(0);

    cstart = Clock::now();
    for (auto it = buckets.begin(); it != buckets.end(); ++it) {
        // quickly sort xit and yit hashX
        for (auto xit = (*it)->hashX->hashed.cbegin();
             xit != (*it)->hashX->hashed.cend(); ++xit) {
            for (auto yit = (*it)->hashY->hashed.cbegin();
                 yit != (*it)->hashY->hashed.cend(); ++yit) {
                // if mapping within the same file, ignore self mapping
                if (xfile == yfile && xit->first == yit->first)
                    continue;
                // check alignment score (post-processing)
                int align_score       = 0;
                int aligned_i         = 0;
                int aligned_j         = 0;
                pair<int, int> cur_xy = make_pair(xit->first, yit->first);
                // pair<int, int> cur_mapped_low = make_pair(-range, -range);

                for (int i = 0; i < xit->second.size(); i++) {
                    for (int j = 0; j < yit->second.size(); j++) {
                        // if the current map point of (x,y) has been covered
                        // by previous alignment, then skip
                        // if (((cur_mapped_low.first + range) >=
                        //         xit->second[i] &&
                        //     cur_mapped_low.first < xit->second[i]) ||
                        //    ((cur_mapped_low.second + range) >=
                        //         yit->second[j] &&
                        //     cur_mapped_low.second < yit->second[j])) {
                        //    continue;
                        //}

                        // 1.30.2023 - since I changed the alignment function
                        // the returned item also include the actual alignment
                        // TESTING OUT its impact on performance
                        // pair<int, int> val = make_pair(xit->second[i],
                        //        yit->second[j]);
                        // results[cur_xy].push_back(val);
                        auto ret    = alignment(X[xit->first], Y[yit->first],
                                                xit->second[i], yit->second[j],
                                                range, range, false);
                        align_score = get<0>(ret);
                        // align_score = max(align_score, get<0>(ret));
                        //  align_score = max(
                        //      align_score,
                        //      alignment(X[xit->first], Y[yit->first],
                        //                xit->second[i], yit->second[j],
                        //                range, range, false));

                        // if alignment score is good, then end the loop
                        if (align_score >= align_thres) {
                            aligned_i = i;
                            aligned_j = j;
                            // goto endloop;
                            /* 1.12.2023 - do not end the loop but search
                             * through all pairs (can be slow but let's see) to
                             * avoid not finding all possible pairs
                             */
                            pair<int, int> val;
                            // key = make_pair(xit->first, yit->first);
                            val = make_pair(xit->second[i], yit->second[j]);
                            results[cur_xy].push_back(val);

                            //// update low end of the cur_mapped_low
                            // pair<int, int> cur_low = cur_mapped_low[cur_xy];
                            // cur_mapped_low[cur_xy] =
                            //     make_pair(max(cur_low.first, xit->second[i]),
                            //               max(cur_low.second,
                            //               yit->second[j]));
                        }
                    }
                }
                // endloop:
                //     if (align_score >= align_thres) {
                //         pair<int, int> val;
                //         // key = make_pair(xit->first, yit->first);
                //         val = make_pair(xit->second[aligned_i],
                //                         yit->second[aligned_j]);
                //         results[cur_xy].push_back(val);

                //        // update low end of the cur_mapped_low
                //        pair<int, int> cur_low = cur_mapped_low[cur_xy];
                //        cur_mapped_low[cur_xy] =
                //            make_pair(max(cur_low.first,
                //            xit->second[aligned_i]),
                //                      max(cur_low.second,
                //                      yit->second[aligned_j]));
                //    }

                //    // reports << xit->first << ":" << xit->second << ";"
                //    //    << yit->first << ":" << yit->second << endl;
            }
        }
        xmer += (*it)->x.size();
        ymer += (*it)->y.size();
    }
    cend = Clock::now();
    t2   = chrono::duration_cast<chrono::milliseconds>(cend - cstart).count();
    cerr << t2 << " ms" << endl;
    xmer /= buckets.size();
    ymer /= buckets.size();

    // (3.2) map results so they are connected
    cerr << "Mapping and writing results in buckets... ";
    cstart = Clock::now();
    reports.open(output_name, ofstream::out);
    unsigned int num_passed_filter = 0, num_checked = 0;
    unsigned int num_passed_kmers = 0, num_checked_kmers = 0;

    /////////// TEMP output to see what are the mapped pairs
    // for (auto paired = results.begin(); paired != results.end(); ++paired) {
    //     sort(paired->second.begin(), paired->second.end(),
    //     sort_by_second_then_first); int x = paired->first.first + 1; int y =
    //     paired->first.second + 1; reports << x << "," << y << ": " <<
    //     paired->second << endl;
    // }

    for (auto paired = results.begin(); paired != results.end(); ++paired) {
        map<int, char>
            unchecked_positions; // keys are ordered in ascending order
        int x           = paired->first.first + 1;
        int y           = paired->first.second + 1;
        int maps_length = paired->second.size();

        // sort the vector by second then first position in pair, ascendingly
        // sort(paired->second.begin(), paired->second.end(), sort_by_second);
        sort(paired->second.begin(), paired->second.end(),
             sort_by_second_then_first);
        unsigned int q_start = 0, q_end = 0, t_start = 0, t_end = 0;

        // vector of int to store aligned positions
        // every 4 ints: q_start, q_end, t_start, t_end
        // vector<int> alignments;

        /***************************
         * rewriting the logic for connecting mapped elements between x and y
         * new pipeline should work in worst case O(n^2), n is the maps_length
         * 1. need to look at each possible start position (i.e., q_start,
         * t_start) 1.2 starting from the start position, check all following
         *          "unchecked" position to possibly expand q_end, t_end.
         *      1.3 all positions used for expansion will be "checked" and
         *          removed (i.e., won't be used as or in other start position)
         * 2. repeat 1. until all possible start positions are checked
         */

        // populate unchecked positions
        for (int _i = 0; _i < maps_length; _i++) {
            unchecked_positions.insert({_i, '.'});
        }

        // while unchecked_positions still have items, continue searching
        while (!unchecked_positions.empty()) {
            // get the starting position and "check it out"
            int start_pos          = unchecked_positions.begin()->first;
            double pct_shared_kmer = 0.0;
            bool passed;
            unchecked_positions.erase(unchecked_positions.begin());

            // initializing q_start, q_end, t_start, t_end
            q_start = paired->second[start_pos].first;
            q_end   = q_start;
            t_start = paired->second[start_pos].second;
            t_end   = t_start;

            // loop over the unchecked positions to connect maps together
            for (auto entry = unchecked_positions.begin();
                 entry != unchecked_positions.end();
                 /* no increment */) {
                unsigned int q_pos = (paired->second[entry->first]).first;
                unsigned int t_pos = (paired->second[entry->first]).second;

                // if q_pos is close to q_end AND t_pos is close to t_end,
                // then we extend both ends and "check out" the current index
                if (((t_pos - t_end) <= search_range) &&
                    ((q_pos - q_end) <= search_range)) {
                    q_end = q_pos;
                    t_end = t_pos;
                    unchecked_positions.erase(entry++);
                }
                // if both q_pos/t_pos are out of reach, we can terminate the
                // for-loop early since it is not possible for the following
                // elements to connect
                else if (((t_pos - t_end) > search_range) &&
                         ((q_pos - q_end) > search_range)) {
                    break;
                }
                // otherwise we continue searching
                else {
                    ++entry;
                }
            }
            // cout << x << "," << y << ": " << q_start << " " << q_end << " "
            // << t_start
            //     << " " << t_end << endl;

            // now, add the connected map if it has length exceeds 100 bp
            if (((t_end - t_start) >= map_len_thres) &&
                ((q_end - q_start) >= map_len_thres)) {
                // if long mapping (e.g., > 1000 bp), do shared kmer filtering
                if ((t_end - t_start >= long_len_thres) || 
                    (q_end - q_start >= long_len_thres)) {
                    passed = filter_and_write(reports, false, kmer_threshold,
                            x, y, q_start, q_end, t_start, t_end);
                    num_checked_kmers++;
                    if (passed) num_passed_kmers++;
                } else {
                    passed = filter_and_write(reports, true, id_threshold, x, y,
                                              q_start, q_end, t_start, t_end);
                    num_checked++;
                    if (passed) num_passed_filter++;
                }
                // else
                //     cerr << x << "," << y << "," << q_start << ","
                //         << q_end << "," << t_start << "," << t_end
                //         << " --> " << shared_kmer << endl;
            }
        }
        /***************************/

        //// iterate through all elements in a paired (x, y)
        // for (auto entry = paired->second.cbegin();
        //      entry != paired->second.cend(); ++entry) {
        //     int q_pos = entry->first, t_pos = entry->second;

        //    // initialization
        //    if (q_start == -1) {
        //        q_start = q_pos, q_end = q_pos;
        //    }
        //    if (t_start == -1) {
        //        t_start = t_pos, t_end = t_pos;
        //    }

        //    // case 1 - X bp for search close one (combine together)
        //    if (t_pos - t_end <= search_range && q_pos - q_end <=
        //    search_range) {
        //        q_end = q_pos;
        //        t_end = t_pos;
        //    } else {
        //        // if next one is beyond range AND
        //        // if length has exceeds 100 bp, append such alignment
        //        // as true one
        //        if ((t_end - t_start >= map_len_thres)
        //                && (q_end - q_start >= map_len_thres)) {
        //            /* 1.10.2023 - UPDATE to have an additional post-filtering
        //             * step. That is, run a global alignment between the
        //             mapped
        //             * region to see if the alignment identity is above a set
        //             * threshold, e.g., 60%. If yes then write to output.
        //             */
        //            filter_and_write(reports, id_threshold, x, y, q_start,
        //                             q_end, t_start, t_end);
        //            /* 1.10.2023 - UPDATE ends */

        //            //reports << x << " " << y << " " << q_start << " " <<
        //            //    q_end << " " << t_start << " " << t_end << endl;
        //        }
        //        // reset alignment positions
        //        q_start = q_pos, q_end = q_pos, t_start = t_pos, t_end =
        //        t_pos;
        //    }
        //}
        //// add in the last one checked if it exists
        // if (q_start != -1 && t_start != -1) {
        //     if ((t_end - t_start >= map_len_thres)
        //             && (q_end - q_start >= map_len_thres)) {
        //         filter_and_write(reports, id_threshold, x, y, q_start,
        //                          q_end, t_start, t_end);
        //         //reports << x << " " << y << " " << q_start << " " << q_end
        //         //        << " " << t_start << " " << t_end << endl;
        //     }
        // }
    }
    cend = Clock::now();
    t3   = chrono::duration_cast<chrono::milliseconds>(cend - cstart).count();
    cerr << t3 << " ms" << endl;

    // (4) writing to local files
    // cerr << "Writing results ... ";
    // cstart = Clock::now();
    // reports.open(output_name, ofstream::out);
    // for (auto it = results.cbegin(); it != results.cend(); ++it) {
    //    // reports << it->first << ":" << *(it->second.begin()) << " -> "
    //    //        << (it->second.rbegin())->first + range / 2 << ","
    //    //        << (it->second.rbegin())->second + range / 2 << endl;
    //    reports << it->first << ":" << it->second << endl;
    //}
    // cend = Clock::now();
    // t3   = chrono::duration_cast<chrono::milliseconds>(cend -
    // cstart).count(); cerr << t3 << " ms" << endl;

    /* LOCAL ALIGNMENT END */

    int overall = t0 + t1 + t2 + t3;

    // all outputs
    reports << "align_multi=" << align_multi << ","
            << "id_threshold=" << id_threshold << ","
            << "kmer_threshold=" << kmer_threshold << ","
            << "search_range=" << search_range << ","
            << "map_len_thres=" << map_len_thres << ","
            << "long_len_thres=" << long_len_thres << ","
            << "num_checked_alignment=" << num_checked << ","
            << "num_passed_alignment=" << num_passed_filter << ","
            << "num_checked_kmers=" << num_checked_kmers << ","
            << "num_passed_kmers=" << num_passed_kmers << ","
            << add_threshold << "," << kill_threshold << "," << overall << ","
            << t0 << "," << t1 << "," << t2 << "," << t3 << "," << total_p
            << "," << total_px << "," << total_py << "," << total_q << ","
            << num_bucket << "," << total_nodes << "," << xmer << "," << ymer
            << endl;
    reports.close();
    return 0;
}
