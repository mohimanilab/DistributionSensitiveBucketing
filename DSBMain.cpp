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

// simple struct for kmer tree of X/Y
struct kmernode {
    bool visited = 0;
    vector<char> kmer;
    map<int, vector<int>> hashed;

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
map<pair<vector<char>, vector<char>>, trn *> buckets;

static double eps;
static double p_ins;
static double p_del;
static double p_match;

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
        cout << "X file not found in " << xfile << endl;
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
        cout << "Y file not found in " << yfile << endl;
        exit(1);
    }
    if (fy.size() > 0)
        Y.push_back(fy);
    data_y.close();
}

// alignment method for comparing two given data points
int alignment(vector<char> &x, vector<char> &y, int xstart, int ystart,
              int range) {
    // rows are x, columns are y
    int xend = 0;
    int yend = 0;
    if (range != -1) {
        xend = min(xstart + range + 1, int(x.size() + 1));
        yend = min(ystart + range + 1, int(y.size() + 1));
    } else {
        xend = x.size() + 1;
        yend = y.size() + 1;
    }
    vector<vector<int>> matrix(xend - xstart, vector<int>(yend - ystart, 0));

    // initialize first column/row
    for (int i = 1; i < xend - xstart; i++) {
        matrix[i][0] = matrix[i - 1][0] - 1;
    }
    for (int j = 1; j < yend - ystart; j++) {
        matrix[0][j] = matrix[0][j - 1] - 1;
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
    return matrix[xend - xstart - 1][yend - ystart - 1];
}

// function to print needed arguments/documentation for this program
void print_helper() {
    printf("Usage: ./DSBMain -q [query file] -r [reference file] -i [insertion "
           "rate] -d "
           "[deletion rate] -m [mutation rate] -a [add threshold] "
           "-k [kill threshold] -o [name] {-vh}\n\n");
    printf("-h\t\tPrint this block of information.\n");
    printf("-v\t\tVerbose mode.\n");
    printf("-q [path/to/q]\tPath to the query file.\n");
    printf("-r [path/to/r]\tPath to the reference file.\n");
    printf("-i [0<=i<1]\tInsertion rate.\n");
    printf("-d [0<=d<1]\tDeletion rate.\n");
    printf("-m [0<=e<0.5]\tMutation rate when neither insertion/deletion "
           "happens.\n");
    printf("-a [a>0]\tThreshold for a node to be considered a bucket.\n");
    printf("-k [k>0]\tThreshold for a node to be pruned.\n");
    printf("-o [name]\tSpecify output file name.\n\n");
}

// function to get total fp
long calc_fp_sum(map<pair<int, int>, long> &mismatches) {
    long total = 0;
    for (auto it = mismatches.cbegin(); it != mismatches.cend(); ++it) {
        total += it->second;
    }
    return total;
}

// method to pre-allocate a tree for X or Y-mer so that later we only
// need to refer to these nodes
void make_kmer_tree(vector<vector<char>> &target, int ind, int b,
                    map<vector<char>, kmernode *> &mapping) {
    // cerr << "starting making kmer tree..." << endl;

    // root, init kmer
    kmernode *root = new kmernode();
    vector<char> init_kmer;
    root->kmer = init_kmer;

    // root, init hashed
    map<int, vector<int>> init_hashed;
    for (int i = 0; i < target.size(); i++) {
        /* to insert every position into the kmer tree */
        for (int j = 0; j < target[i].size() - 10; j++) {
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

    // iterate till no more nodes can be added (the new nodes will have no
    // kmers inside)
    int total = 0;
    while (!q.empty()) {
        kmernode *cur = q.front();
        q.pop();
        total++;
        // cerr << "on node " << cur->kmer << endl;

        vector<char> cur_kmer = cur->kmer;

        // else we continue on expanding the tree
        for (int i = 0; i < 4; i++) {
            // copy over the current alphabet and expand it
            vector<char> next_kmer(cur_kmer);
            next_kmer.push_back(alphabet[i]);

            // create new node with expanded kmer
            kmernode *next_kmernode = new kmernode();
            next_kmernode->kmer     = next_kmer;

            // iterate through current hashed sequences
            // push in the ones that matched current kmer
            for (auto it = cur->hashed.cbegin(); it != cur->hashed.cend();
                 ++it) {
                for (int m = 0; m < it->second.size(); m++) {
                    // failsafe, if index out of bound
                    if (it->second[m] >= X[it->first].size())
                        continue;

                    // check if the current position fits the last position
                    // of next_kmer (assuming we have everything correct before)
                    if (X[it->first][it->second[m]] == next_kmer.back()) {
                        next_kmernode->hashed[it->first].push_back(
                            it->second[m] + 1);
                    }
                }
            }

            // finish up by pushing the next kmernode to the map&queue
            // IF the node is not empty
            if (next_kmernode->hashed.size() > 0 &&
                next_kmernode->kmer.size() < 30) {
                // cerr << next_kmernode->hashed.cbegin()->first << endl;
                mapping[next_kmer] = next_kmernode;
                q.push(next_kmernode);
            } else
                delete next_kmernode;
        }
    }
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
    while (!q.empty()) {
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
        if (working->hashX->size() == 0)
            continue; // if no more X prefixes
        if (working->hashY->size() == 0)
            continue; // if no more Y prefixes
        // 2) p/q > threshold_1, add to hash tree
        if ((working->p / working->q) > add_threshold) {
            pair<vector<char>, vector<char>> path =
                make_pair(working->x, working->y);
            trn *bucket   = new trn(working);
            buckets[path] = bucket;
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
    string xfile       = "";
    string yfile       = "";
    string output_name = "output.txt";
    double add_threshold, kill_threshold;
    double align_multi = 0.35;
    int range          = 12;
    int verbose        = 0;

    // argument parsing
    int error_flag =
        !(argc == 15 || argc == 16 || argc == 17 || argc == 18 || argc == 19);
    char opt;
    while ((opt = getopt(argc, argv, "hvr:q:i:d:m:a:k:o:")) != -1) {
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

    // check error flag for arguments
    if (error_flag || add_threshold <= 0 || kill_threshold <= 0) {
        print_helper();
        exit(error_flag);
    }

    p_match = 1.0 - p_ins - p_del;

    if (verbose) {
        printf("X: %s\tY: %s\nInsertion: %.3f\tDeletion: %.3f\t\tMutation: "
               "%.3f\nThreshold 1: %.1f\tThreshold 2: %.1f\n\n",
               xfile.c_str(), yfile.c_str(), p_ins, p_del, eps * 12,
               add_threshold, kill_threshold);
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
    map<pair<int, int>, set<pair<int, int>>> results;
    map<pair<int, int>, pair<int, int>> cur_mapped_low;
    ofstream reports;
    int t0           = 0;
    int t1           = 0;
    int t2           = 0;
    int t3           = 0;
    double xmer      = 0.0;
    double ymer      = 0.0;
    long num_bucket  = 0;
    long total_nodes = 0;
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
    cerr << "Making data-dependent tree for buckets... ";
    cstart   = Clock::now();
    auto tmp = indpt_tree(kill_threshold, add_threshold, Xmer_tree, Ymer_tree);
    cend     = Clock::now();
    t1 = chrono::duration_cast<chrono::milliseconds>(cend - cstart).count();
    cerr << t1 << " ms" << endl;

    num_bucket  = get<0>(tmp);
    total_nodes = get<1>(tmp);
    total_p     = get<2>(tmp);
    total_px    = get<3>(tmp);
    total_py    = get<4>(tmp);
    total_q     = get<5>(tmp);

    // (3) collect results from all buckets
    int align_thres = int(align_multi * range);
    cerr << "Collecting results (bucket size: " << buckets.size() << ")... ";
    cstart = Clock::now();
    for (auto it = buckets.cbegin(); it != buckets.cend(); ++it) {
        for (auto xit = it->second->hashX->hashed.cbegin();
             xit != it->second->hashX->hashed.cend(); ++xit) {
            for (auto yit = it->second->hashY->hashed.cbegin();
                 yit != it->second->hashY->hashed.cend(); ++yit) {
                // if mapping within the same file, ignore self mapping
                if (xfile == yfile && xit->first == yit->first)
                    continue;
                // check alignment score (post-processing)
                int align_score = 0;
                int aligned_i, aligned_j;
                for (int i = 0; i < xit->second.size(); i++) {
                    for (int j = 0; j < yit->second.size(); j++) {
                        // if the current map point of (x,y) has been covered
                        // by previous alignment, then skip
                        pair<int, int> cur_xy =
                            make_pair(xit->first, yit->first);
                        // if the pair is checked already (exists in the
                        // result), then skip
                        // if (results.find(cur_xy) != results.end()) {
                        //    goto endloop;
                        //}
                        if (((cur_mapped_low[cur_xy].first + range) >=
                                 xit->second[i] &&
                             cur_mapped_low[cur_xy].first < xit->second[i]) ||
                            ((cur_mapped_low[cur_xy].second + range) >=
                                 yit->second[j] &&
                             cur_mapped_low[cur_xy].second < yit->second[j])) {
                            continue;
                        }

                        align_score = max(
                            align_score,
                            alignment(X[xit->first], Y[yit->first],
                                      xit->second[i], yit->second[j], range));

                        // if alignment score is good, then end the loop
                        if (align_score >= align_thres) {
                            aligned_i = i;
                            aligned_j = j;
                            goto endloop;
                        }
                    }
                }
            endloop:
                if (align_score >= align_thres) {
                    pair<int, int> key, val;
                    key = make_pair(xit->first, yit->first);
                    val = make_pair(xit->second[aligned_i],
                                    yit->second[aligned_j]);
                    results[key].insert(val);

                    // update low end of the cur_mapped_low
                    pair<int, int> cur_low = cur_mapped_low[key];
                    cur_mapped_low[key] =
                        make_pair(max(cur_low.first, xit->second[aligned_i]),
                                  max(cur_low.second, yit->second[aligned_j]));
                }

                // reports << xit->first << ":" << xit->second << ";"
                //    << yit->first << ":" << yit->second << endl;
            }
        }
        xmer += (it->first.first).size();
        ymer += (it->first.second).size();
    }
    cend = Clock::now();
    t2   = chrono::duration_cast<chrono::milliseconds>(cend - cstart).count();
    cerr << t2 << " ms" << endl;
    xmer /= buckets.size();
    ymer /= buckets.size();

    // (4) writing to local files
    cerr << "Writing results ... ";
    cstart = Clock::now();
    reports.open(output_name, ofstream::out);
    for (auto it = results.cbegin(); it != results.cend(); ++it) {
        // reports << it->first << ":" << *(it->second.begin()) << " -> "
        //        << (it->second.rbegin())->first + range / 2 << ","
        //        << (it->second.rbegin())->second + range / 2 << endl;
        reports << it->first << ":" << it->second << endl;
    }
    cend = Clock::now();
    t3   = chrono::duration_cast<chrono::milliseconds>(cend - cstart).count();
    cerr << t3 << " ms" << endl;

    /* LOCAL ALIGNMENT END */

    int overall = t0 + t1 + t2 + t3;

    // all outputs
    reports << add_threshold << "," << kill_threshold << "," << overall << ","
            << t0 << "," << t1 << "," << t2 << "," << t3 << "," << total_p
            << "," << total_px << "," << total_py << "," << total_q << ","
            << num_bucket << "," << total_nodes << "," << xmer << "," << ymer
            << endl;
    reports.close();
    return 0;
}
