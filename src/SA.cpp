#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <random>
#include <chrono>
#include <cmath>
#include <algorithm>

using namespace std;

// global grid size
int R = 0, C = 0;

struct pin_info
{
    string name;
    double x;
    double y;
};

struct logic_info
{
    string name;
    int x;                    // column (0..C-1)
    int y;                    // row    (0..R-1)
    vector<int> connect_nets; // indices of nets that touch this block
};

struct net_info
{
    string name;
    int degree;
    vector<string> connections; // names of blocks or pins
};
struct SAParams
{
    double alpha;          // cost = alpha * CC + (1-alpha) * HPWL
    double T_init;         // initial temperature
    double T_end;          // end temperature
    double cooling;        // cooling rate
    int outer_iter;        // max number of temperature levels
    int moves_per_T_min;   // minimum moves per temperature (baseline)
    double time_limit_sec; // time limit in seconds
    unsigned int seed;     // RNG seed
};

// -------------------- parsing --------------------

void parse_input(
    ifstream &parse,
    vector<logic_info> &logic_blocks_table,
    vector<pin_info> &pin_table,
    vector<net_info> &net_table,
    unordered_map<string, int> &name2idx,
    unordered_map<string, int> &pin2idx,
    unordered_map<string, int> &net2idx)
{
    string line;
    if (!getline(parse, line))
    {
        cerr << "Empty input file or failed to read header.\n";
        exit(1);
    }

    istringstream info(line);
    int logic_block_R, logic_block_C, num_logic_blocks, num_io_pins, num_nets;
    info >> logic_block_R >> logic_block_C >> num_logic_blocks >> num_io_pins >> num_nets;

    R = logic_block_R;
    C = logic_block_C;

    if (R <= 0 || C <= 0)
    {
        cerr << "Invalid grid size R=" << R << " C=" << C << endl;
        exit(1);
    }

    // logic blocks
    logic_blocks_table.reserve(num_logic_blocks);
    for (int i = 0; i < num_logic_blocks; ++i)
    {
        string name_b;
        getline(parse, line);
        if (line.empty())
        {
            --i;
            continue;
        }
        istringstream logic_inf(line);
        logic_inf >> name_b;
        logic_info inst;
        inst.name = name_b;
        inst.x = 0;
        inst.y = 0;
        logic_blocks_table.push_back(inst);
        name2idx[name_b] = i;
    }

    // IO pins
    pin_table.reserve(num_io_pins);
    for (int i = 0; i < num_io_pins; ++i)
    {
        string name_p;
        double x_p, y_p;
        getline(parse, line);
        if (line.empty())
        {
            --i;
            continue;
        }
        istringstream pin_inf(line);
        pin_inf >> name_p >> x_p >> y_p;
        pin_info inst;
        inst.name = name_p;
        inst.x = x_p;
        inst.y = y_p;
        pin_table.push_back(inst);
        pin2idx[name_p] = i;
    }

    // nets
    net_table.reserve(num_nets);
    for (int i = 0; i < num_nets; ++i)
    {
        string name_n;
        int degree;
        getline(parse, line);
        if (line.empty())
        {
            --i;
            continue;
        }
        istringstream net_inf(line);
        net_inf >> name_n >> degree;
        net_info inst;
        inst.name = name_n;
        inst.degree = degree;
        inst.connections.reserve(degree);
        for (int j = 0; j < degree; ++j)
        {
            string item;
            net_inf >> item;
            inst.connections.push_back(item);
            auto it = name2idx.find(item);
            if (it != name2idx.end())
            {
                // this net connects to movable block it->second
                logic_blocks_table[it->second].connect_nets.push_back(i);
            }
        }
        net_table.push_back(inst);
        net2idx[name_n] = i;
    }
}

// -------------------- initial placement with occupancy --------------------

// occ[y][x] = index of block in logic_blocks_table, or -1 if empty
void initial_place(vector<logic_info> &logic_blocks_table,
                   vector<vector<int>> &occ)
{
    int idx = 0;
    for (int y = 0; y < R; ++y)
    {
        for (int x = 0; x < C; ++x)
        {
            if (idx < (int)logic_blocks_table.size())
            {
                occ[y][x] = idx;
                logic_blocks_table[idx].x = x;
                logic_blocks_table[idx].y = y;
                ++idx;
            }
            else
            {
                occ[y][x] = -1; // empty
            }
        }
    }
}

// -------------------- cost computation: HPWL + CC --------------------

// compute total HPWL and congestion coefficient CC (coverage-based)
void compute_hpwl_CC(
    const vector<logic_info> &logic_blocks_table,
    const vector<pin_info> &pin_table,
    const vector<net_info> &net_table,
    const unordered_map<string, int> &name2idx,
    const unordered_map<string, int> &pin2idx,
    double &total_hpwl_out,
    double &CC_out)
{
    // coverage grid U[y][x]
    vector<vector<int>> U(R, vector<int>(C, 0));

    double total_hpwl = 0.0;

    for (const auto &net : net_table)
    {
        bool first = true;
        double xmin = 0.0, xmax = 0.0;
        double ymin = 0.0, ymax = 0.0;

        for (const auto &item : net.connections)
        {
            auto it_blk = name2idx.find(item);
            auto it_pin = pin2idx.find(item);

            if (it_blk != name2idx.end())
            {
                int idx = it_blk->second;
                double xb = static_cast<double>(logic_blocks_table[idx].x);
                double yb = static_cast<double>(logic_blocks_table[idx].y);
                double this_xmin = xb;
                double this_xmax = xb + 1.0;
                double this_ymin = yb;
                double this_ymax = yb + 1.0;

                if (first)
                {
                    xmin = this_xmin;
                    xmax = this_xmax;
                    ymin = this_ymin;
                    ymax = this_ymax;
                    first = false;
                }
                else
                {
                    xmin = std::min(xmin, this_xmin);
                    xmax = std::max(xmax, this_xmax);
                    ymin = std::min(ymin, this_ymin);
                    ymax = std::max(ymax, this_ymax);
                }
            }
            else if (it_pin != pin2idx.end())
            {
                int idx = it_pin->second;
                double xp = pin_table[idx].x;
                double yp = pin_table[idx].y;
                double this_xmin = xp;
                double this_xmax = xp;
                double this_ymin = yp;
                double this_ymax = yp;

                if (first)
                {
                    xmin = this_xmin;
                    xmax = this_xmax;
                    ymin = this_ymin;
                    ymax = this_ymax;
                    first = false;
                }
                else
                {
                    xmin = std::min(xmin, this_xmin);
                    xmax = std::max(xmax, this_xmax);
                    ymin = std::min(ymin, this_ymin);
                    ymax = std::max(ymax, this_ymax);
                }
            }
        }

        if (first)
        {
            // this net has no valid terminal; ignore
            continue;
        }

        double hpwl = (xmax - xmin) + (ymax - ymin);
        total_hpwl += hpwl;

        // update coverage U: integer sites within bounding box
        int x_start = std::max(0, (int)std::floor(xmin));
        int x_end = std::min(C, (int)std::ceil(xmax)); // [x_start, x_end)
        int y_start = std::max(0, (int)std::floor(ymin));
        int y_end = std::min(R, (int)std::ceil(ymax)); // [y_start, y_end)

        for (int y = y_start; y < y_end; ++y)
        {
            for (int x = x_start; x < x_end; ++x)
            {
                U[y][x] += 1;
            }
        }
    }

    // compute CC = (N * sum(U^2)) / (sum(U)^2)
    long long N = (long long)R * (long long)C;
    long long sumU = 0;
    long long sumU2 = 0;
    for (int y = 0; y < R; ++y)
    {
        for (int x = 0; x < C; ++x)
        {
            int v = U[y][x];
            sumU += v;
            sumU2 += 1LL * v * v;
        }
    }

    double CC;
    if (sumU == 0)
    {
        CC = 1.0; // no net coverage -> treat as perfectly uniform
    }
    else
    {
        CC = (double)N * (double)sumU2 / ((double)sumU * (double)sumU);
    }

    total_hpwl_out = total_hpwl;
    CC_out = CC;
}

// -------------------- simulated annealing --------------------

void simulated_annealing(
    vector<logic_info> &logic_blocks_table,
    const vector<pin_info> &pin_table,
    const vector<net_info> &net_table,
    const unordered_map<string, int> &name2idx,
    const unordered_map<string, int> &pin2idx,
    vector<vector<int>> &occ,
    const SAParams &params)
{
    // RNG & timer
    std::mt19937 rng(params.seed);
    std::uniform_real_distribution<double> urand(0.0, 1.0);

    auto start_time = std::chrono::high_resolution_clock::now();

    // --- initial HPWL & CC (for normalization only) ---
    double hpwl0, CC0;
    compute_hpwl_CC(logic_blocks_table, pin_table, net_table,
                    name2idx, pin2idx, hpwl0, CC0);

    if (hpwl0 <= 0.0) hpwl0 = 1.0;
    if (CC0   <= 0.0) CC0   = 1.0;

    auto cost_from = [&](double hpwl, double CC)
    {
        double f_hpwl = hpwl / hpwl0;
        double f_CC   = CC   / CC0;
        return params.alpha * f_CC + (1.0 - params.alpha) * f_hpwl;
    };

    // --- current solution cost ---
    double cur_hpwl, cur_CC;
    compute_hpwl_CC(logic_blocks_table, pin_table, net_table,
                    name2idx, pin2idx, cur_hpwl, cur_CC);
    double cur_cost = cost_from(cur_hpwl, cur_CC);

    // --- best solution backup ---
    vector<logic_info> best_solution = logic_blocks_table;
    double best_hpwl = cur_hpwl;
    double best_CC   = cur_CC;
    double best_cost = cur_cost;

    // SA parameters
    double T = params.T_init;
    const double T_end          = params.T_end;
    const double cooling        = params.cooling;
    const int    outer_iter     = params.outer_iter;
    const int    moves_per_T    = std::max(params.moves_per_T_min,
                                           (int)logic_blocks_table.size());
    const double time_limit_sec = params.time_limit_sec;

    std::uniform_int_distribution<int> x_dist(0, C - 1);
    std::uniform_int_distribution<int> y_dist(0, R - 1);

    // --- SA main loop (textbook style) ---
    for (int it = 0; it < outer_iter && T > T_end; ++it)
    {
        for (int m = 0; m < moves_per_T; ++m)
        {
            auto now = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(now - start_time).count();
            if (elapsed > time_limit_sec)
            {
                it = outer_iter; // break outer
                break;
            }

            // pick two random grid positions
            int x1 = x_dist(rng);
            int y1 = y_dist(rng);
            int x2 = x_dist(rng);
            int y2 = y_dist(rng);

            // skip trivial same cell
            if (x1 == x2 && y1 == y2)
                continue;

            int b1 = occ[y1][x1]; // -1 or block index
            int b2 = occ[y2][x2];

            // skip empty <-> empty
            if (b1 == -1 && b2 == -1)
                continue;

            // backup
            int old_x1 = x1, old_y1 = y1;
            int old_x2 = x2, old_y2 = y2;

            // apply swap in occ and block coordinates
            occ[y1][x1] = b2;
            occ[y2][x2] = b1;
            if (b1 != -1)
            {
                logic_blocks_table[b1].x = x2;
                logic_blocks_table[b1].y = y2;
            }
            if (b2 != -1)
            {
                logic_blocks_table[b2].x = x1;
                logic_blocks_table[b2].y = y1;
            }

            // evaluate new cost
            double new_hpwl, new_CC;
            compute_hpwl_CC(logic_blocks_table, pin_table, net_table,
                            name2idx, pin2idx, new_hpwl, new_CC);
            double new_cost = cost_from(new_hpwl, new_CC);
            double delta    = new_cost - cur_cost;

            bool accept = false;

            // textbook SA rule:
            // downhill: always accept; uphill: e^{-Δ/T}
            if (delta <= 0.0)
            {
                accept = true;
            }
            else
            {
                double prob = std::exp(-delta / T);
                if (urand(rng) < prob)
                    accept = true;
            }

            if (accept)
            {
                cur_cost = new_cost;
                cur_hpwl = new_hpwl;
                cur_CC   = new_CC;

                if (cur_cost < best_cost)
                {
                    best_cost   = cur_cost;
                    best_hpwl   = cur_hpwl;
                    best_CC     = cur_CC;
                    best_solution = logic_blocks_table;
                }
            }
            else
            {
                // rollback swap
                occ[y1][x1] = b1;
                occ[y2][x2] = b2;
                if (b1 != -1)
                {
                    logic_blocks_table[b1].x = old_x1;
                    logic_blocks_table[b1].y = old_y1;
                }
                if (b2 != -1)
                {
                    logic_blocks_table[b2].x = old_x2;
                    logic_blocks_table[b2].y = old_y2;
                }
            }
        }

        T *= cooling;

        // 如果要 TEMP log，留著這行：
        cerr << "[TEMP] T=" << T
             << " best_cost=" << best_cost
             << " best_HPWL=" << best_hpwl
             << " best_CC=" << best_CC << "\n";
    }

    // restore best solution
    logic_blocks_table = best_solution;

    // 最後把 best 的 HPWL / CC 印出來
    compute_hpwl_CC(logic_blocks_table, pin_table, net_table,
                    name2idx, pin2idx, best_hpwl, best_CC);
    cerr << "[SA] best_cost = " << best_cost
         << "  HPWL = " << best_hpwl
         << "  CC = " << best_CC << endl;
}


// -------------------- output --------------------

void write_output(const string &out_file,
                  const vector<logic_info> &logic_blocks_table)
{
    ofstream fout(out_file);
    if (!fout)
    {
        cerr << "Failed to open output file: " << out_file << endl;
        return;
    }
    for (const auto &lb : logic_blocks_table)
    {
        fout << lb.name << " " << lb.x << " " << lb.y << "\n";
    }
    fout.close();
}

// -------------------- main --------------------

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>\n";
        return 1;
    }

    ifstream parse_in(argv[1]);
    if (!parse_in.is_open())
    {
        cerr << "Failed to open input file: " << argv[1] << endl;
        return 1;
    }

    vector<logic_info> logic_blocks_table;
    vector<pin_info> pin_table;
    vector<net_info> net_table;
    unordered_map<string, int> name2idx;
    unordered_map<string, int> pin2idx;
    unordered_map<string, int> net2idx;

    parse_input(parse_in,
                logic_blocks_table,
                pin_table,
                net_table,
                name2idx,
                pin2idx,
                net2idx);

    // occupancy grid
    vector<vector<int>> occ(R, vector<int>(C, -1));
    initial_place(logic_blocks_table, occ);
    double init_hpwl, init_cc;
    compute_hpwl_CC(logic_blocks_table, pin_table, net_table,
                    name2idx, pin2idx, init_hpwl, init_cc);
    cerr << "[INIT] HPWL=" << init_hpwl << " CC=" << init_cc << "\n";

    SAParams params;
    params.alpha = 0.05;  // CC 0.4, HPWL 0.6
    params.T_init = 2.0; // 高一點，方便跳出爛區
    params.T_end = 1e-4;
    params.cooling = 0.95;
    params.outer_iter = 120;       // 多一點溫度層
    params.moves_per_T_min = 1200; // 每層多一點 move
    params.time_limit_sec = 180.0; // 留時間給 I/O
    params.seed = 1234567;

    // run SA to improve placement
    simulated_annealing(logic_blocks_table,
                        pin_table,
                        net_table,
                        name2idx,
                        pin2idx,
                        occ, params);

    // write final placement
    write_output(string(argv[2]), logic_blocks_table);

    return 0;
}
