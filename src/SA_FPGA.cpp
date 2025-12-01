#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <math.h>
#include <random>
#include <chrono>
#include <unordered_set>
#include <algorithm>
#include <climits>

std::mt19937 rng(123456);

using namespace std;
using ns = std::chrono::nanoseconds;
using Clock = std::chrono::high_resolution_clock;

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
    int x;
    int y;
    vector<int> connect_nets;
};

struct rec
{
    double xmax;
    double xmin;
    double ymax;
    double ymin;
};

struct NetTerm
{
    bool is_pin; // false => logic block, true => IO pin
    int idx;     // index in logic_blocks_table 或 pin_table
};

struct net_info
{
    string name;
    int degree;
    vector<string> connections;
    vector<NetTerm> terms;
};

struct net_hpwl_cong
{
    double hpwl;
    rec box;
};

struct change_info
{
    string nameA, nameB;
    int Ax, Ay;
    int Bx, By;
};

struct OptimalRegion
{
    int lx, ly;
    int ux, uy;
};

void initial_place(vector<vector<string>> &array2D, int R, int C,
                   vector<logic_info> &logic_blocks_table)
{
    int N = (int)logic_blocks_table.size();
    int idx = 0;
    for (int i = 0; i < R && idx < N; ++i)
    {
        for (int j = 0; j < C && idx < N; ++j)
        {
            array2D[i][j] = logic_blocks_table[idx].name;
            logic_blocks_table[idx].x = i;
            logic_blocks_table[idx].y = j;
            ++idx;
        }
    }
}

rec calculate_bounding_box(const net_info &net,
                           const vector<logic_info> &logic_blocks_table,
                           const vector<pin_info> &pin_table)
{
    bool first = true;
    double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;

    for (const auto &term : net.terms)
    {
        if (term.idx < 0) continue;

        double this_xmin, this_xmax, this_ymin, this_ymax;

        if (!term.is_pin)
        {
            const auto &blk = logic_blocks_table[term.idx];
            double xb = (double)blk.x;
            double yb = (double)blk.y;
            this_xmin = xb;
            this_xmax = xb + 1.0;
            this_ymin = yb;
            this_ymax = yb + 1.0;
        }
        else
        {
            const auto &p = pin_table[term.idx];
            double xp = p.x;
            double yp = p.y;
            this_xmin = xp;
            this_xmax = xp;
            this_ymin = yp;
            this_ymax = yp;
        }

        if (first)
        {
            xmin = this_xmin; xmax = this_xmax;
            ymin = this_ymin; ymax = this_ymax;
            first = false;
        }
        else
        {
            xmin = min(xmin, this_xmin);
            xmax = max(xmax, this_xmax);
            ymin = min(ymin, this_ymin);
            ymax = max(ymax, this_ymax);
        }
    }

    rec Rect;
    if (first)
    {
        Rect.xmin = Rect.xmax = Rect.ymin = Rect.ymax = 0.0;
        return Rect;
    }
    Rect.xmin = xmin; Rect.xmax = xmax;
    Rect.ymin = ymin; Rect.ymax = ymax;
    return Rect;
}

void compute_hpwl_cong(const vector<logic_info> &logic_blocks_table,
                       const vector<pin_info> &pin_table,
                       const vector<net_info> &net_table,
                       vector<vector<int>> &U,
                       vector<net_hpwl_cong> &net_cost,
                       double &total_hpwl,
                       double &CC)
{
    for (int y = 0; y < R; ++y)
        fill(U[y].begin(), U[y].end(), 0);

    total_hpwl = 0.0;

    for (size_t i = 0; i < net_table.size(); ++i)
    {
        const net_info &net = net_table[i];
        rec box = calculate_bounding_box(net, logic_blocks_table, pin_table);
        double hpwl = (box.xmax - box.xmin) + (box.ymax - box.ymin);
        total_hpwl += hpwl;

        net_cost[i].hpwl = hpwl;
        net_cost[i].box  = box;

        int x_start = max(0, (int)floor(box.xmin));
        int x_end   = min(C, (int)ceil(box.xmax));
        int y_start = max(0, (int)floor(box.ymin));
        int y_end   = min(R, (int)ceil(box.ymax));

        for (int y = y_start; y < y_end; ++y)
        {
            for (int x = x_start; x < x_end; ++x)
            {
                U[y][x] += 1;
            }
        }
    }

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

    if (sumU == 0) CC = 1.0;
    else CC = (double)N * (double)sumU2 / ((double)sumU * (double)sumU);
}

double compute_cost(double total_hpwl, double CC, double lambda)
{
    return total_hpwl * CC;
}

void undo_move(vector<logic_info> &logic_blocks_table,
               unordered_map<string, int> &name2idx,
               vector<vector<string>> &array2D,
               const change_info &chg)
{
    int idxA = name2idx[chg.nameA];
    auto &A = logic_blocks_table[idxA];

    if (chg.nameB != "none")
    {
        int idxB = name2idx[chg.nameB];
        auto &B = logic_blocks_table[idxB];
        A.x = chg.Ax; A.y = chg.Ay;
        B.x = chg.Bx; B.y = chg.By;
        array2D[chg.Ax][chg.Ay] = chg.nameA;
        array2D[chg.Bx][chg.By] = chg.nameB;
    }
    else
    {
        A.x = chg.Ax; A.y = chg.Ay;
        array2D[chg.Ax][chg.Ay] = chg.nameA;
        array2D[chg.Bx][chg.By] = "none";
    }
}

OptimalRegion compute_optimal_region_for_block(
    const logic_info &blk,
    const vector<logic_info> &logic_blocks_table,
    const vector<net_info> &net_table,
    const vector<pin_info> &pin_table,
    const vector<net_hpwl_cong> &net_cost,
    double scale)
{
    (void)logic_blocks_table;
    (void)pin_table;

    double sum_cx = 0.0;
    double sum_cy = 0.0;
    int cnt = 0;

    for (int nid : blk.connect_nets)
    {
        const rec &box = net_cost[nid].box;
        if (box.xmax <= box.xmin && box.ymax <= box.ymin) continue;

        double cx = 0.5 * (box.xmin + box.xmax);
        double cy = 0.5 * (box.ymin + box.ymax);
        sum_cx += cx;
        sum_cy += cy;
        cnt++;
    }

    OptimalRegion region;
    if (cnt == 0)
    {
        region.lx = 0; region.ly = 0;
        region.ux = max(0, R - 1); region.uy = max(0, C - 1);
        return region;
    }

    double mx = sum_cx / cnt;
    double my = sum_cy / cnt;

    int center_x = (int)std::round(mx);
    int center_y = (int)std::round(my);
    center_x = std::max(0, std::min(R - 1, center_x));
    center_y = std::max(0, std::min(C - 1, center_y));

    int base_w = std::max(1, R / 10);
    int base_h = std::max(1, C / 10);

    int half_w = std::max(2, (int)(base_w * scale));
    int half_h = std::max(2, (int)(base_h * scale));

    region.lx = std::max(0, center_x - half_w);
    region.ux = std::min(R - 1, center_x + half_w);
    region.ly = std::max(0, center_y - half_h);
    region.uy = std::min(C - 1, center_y + half_h);

    return region;
}

// 取樣幾個 swap 估 T0（比較接近工業版作法）
double EstimateInitialTemperature(const vector<logic_info> &logic_blocks_table,
                                  const vector<pin_info> &pin_table,
                                  const vector<net_info> &net_table,
                                  double total_hpwl,
                                  double CC,
                                  double lambda)
{
    double base_cost = compute_cost(total_hpwl, CC, lambda);
    int num_cells = (int)logic_blocks_table.size();
    if (num_cells <= 1) return 1.0;

    // 複製一份 blocks 用來做暫時測試
    vector<logic_info> tmp_blocks = logic_blocks_table;
    vector<vector<int>> tmpU(R, vector<int>(C, 0));
    vector<net_hpwl_cong> tmp_net_cost(net_table.size());

    std::uniform_int_distribution<int> dist_block(0, num_cells - 1);
    int num_samples = min(60, num_cells * 2);

    double sum_pos_delta = 0.0;
    int    cnt_pos       = 0;

    for (int s = 0; s < num_samples; ++s)
    {
        int a = dist_block(rng);
        int b = dist_block(rng);
        if (a == b) continue;

        swap(tmp_blocks[a].x, tmp_blocks[b].x);
        swap(tmp_blocks[a].y, tmp_blocks[b].y);

        double t_hpwl = 0.0, t_CC = 1.0;
        compute_hpwl_cong(tmp_blocks, pin_table, net_table, tmpU, tmp_net_cost, t_hpwl, t_CC);
        double new_cost = compute_cost(t_hpwl, t_CC, lambda);
        double delta = new_cost - base_cost;
        if (delta > 0.0)
        {
            sum_pos_delta += delta;
            cnt_pos++;
        }

        // swap back
        swap(tmp_blocks[a].x, tmp_blocks[b].x);
        swap(tmp_blocks[a].y, tmp_blocks[b].y);
    }

    if (cnt_pos == 0) return max(1.0, base_cost * 0.01);

    double avg_pos_delta = sum_pos_delta / (double)cnt_pos;
    double p0 = 0.8;
    double T0 = -avg_pos_delta / std::log(p0);
    if (T0 <= 0.0) T0 = max(1.0, base_cost * 0.01);
    return T0;
}

// 在 U 上套用一個 bbox 的 delta（+1 或 -1），並更新 sumU / sumU2
static inline void apply_rect_delta_on_U(const rec &box,
                                         int delta,
                                         vector<vector<int>> &U,
                                         long long &sumU,
                                         long long &sumU2)
{
    int x_start = max(0, (int)floor(box.xmin));
    int x_end   = min(C, (int)ceil(box.xmax));
    int y_start = max(0, (int)floor(box.ymin));
    int y_end   = min(R, (int)ceil(box.ymax));

    if (x_start >= x_end || y_start >= y_end) return;

    for (int y = y_start; y < y_end; ++y)
    {
        for (int x = x_start; x < x_end; ++x)
        {
            int old_v = U[y][x];
            int new_v = old_v + delta;
            U[y][x] = new_v;
            sumU  += (new_v - old_v);
            sumU2 += 1LL * new_v * new_v - 1LL * old_v * old_v;
        }
    }
}

void SA(vector<logic_info> &logic_blocks_table,
        vector<vector<string>> &array2D,
        vector<pin_info> &pin_table,
        vector<net_info> &net_table,
        unordered_map<string, int> &name2idx,
        unordered_map<string, int> &pin2idx,
        vector<net_hpwl_cong> &net_cost,
        vector<vector<int>> &U,
        double &total_hpwl,
        double &CC)
{
    (void)pin2idx;

    auto start_SA = Clock::now();
    const int TIME_LIMIT_MS = 225000;

    int num_cells = (int)logic_blocks_table.size();
    if (num_cells == 0) return;

    double lambda = 3.0;

    // 先把 sumU / sumU2 算出來
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
    long long Nsite = (long long)R * (long long)C;

    double cur_cost = compute_cost(total_hpwl, CC, lambda);
    double best_cost = cur_cost;
    auto best_logic_blocks = logic_blocks_table;
    auto best_array2D = array2D;
    double best_hpwl = total_hpwl;
    double best_CC = CC;

    cerr << "[INIT] HPWL=" << total_hpwl
         << " CC=" << CC
         << " cost=" << cur_cost << endl;

    double T = EstimateInitialTemperature(logic_blocks_table, pin_table, net_table,
                                          total_hpwl, CC, lambda);
    double initial_T = T;
    double T_end  = 1e-5;


    int base_iter;
    if (num_cells < 3000)
    {
        base_iter = 200000;
    }
    else if (num_cells < 8000)
    {
        base_iter = std::min(num_cells * 30, 100000);
    }
    else
    {
        base_iter = 50000;
    }
    int iter_per_T = base_iter;

    std::uniform_real_distribution<double> dist01(0.0, 1.0);
    std::uniform_int_distribution<int> dist_block(0, num_cells - 1);
    std::uniform_int_distribution<int> distR(0, R - 1);
    std::uniform_int_distribution<int> distC(0, C - 1);

    int netN = (int)net_table.size();
    std::vector<int> net_mark(netN, 0);
    int net_tag = 0;

    double current_region_prob = 0.3;
    const double max_region_prob = 0.95;
    const double min_region_prob = 0.05;
    const double prob_step_up   = 0.02;
    const double prob_step_down = 0.02;

    int round_cnt = 0;
    int low_acc_rounds = 0;

    cerr << "[SA] Estimated T0 = " << T << " Iter/T = " << iter_per_T << endl;

    while (T > T_end)
    {
        auto now = Clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_SA);
        double time_ratio = (double)duration.count() / TIME_LIMIT_MS;
        if (time_ratio > 1.0)
        {
            cerr << "[SA] Time limit reached.\n";
            break;
        }

        double scale_factor = T / initial_T;
        double window_scale = std::max(0.01, std::min(1.0, sqrt(max(scale_factor, 1e-8))));

        int accepted_moves = 0;

        for (int it = 0; it < iter_per_T; ++it)
        {
            net_tag++;
            if (net_tag == INT_MAX)
            {
                std::fill(net_mark.begin(), net_mark.end(), 0);
                net_tag = 1;
            }

            int idxA = dist_block(rng);
            logic_info &A = logic_blocks_table[idxA];
            int Ax = A.x;
            int Ay = A.y;

            int Bx, By;
            bool use_region = (dist01(rng) < current_region_prob);
            OptimalRegion region;

            if (use_region)
            {
                region = compute_optimal_region_for_block(A, logic_blocks_table, net_table, pin_table, net_cost, window_scale);
                if (region.lx == region.ux && region.ly == region.uy &&
                    region.lx == Ax && region.ly == Ay)
                {
                    use_region = false;
                }
            }

            int search_rx = std::max(1, (int)std::round(window_scale * R ));
            int search_ry = std::max(1, (int)std::round(window_scale * C ));

            do {
                if (use_region) {
                    std::uniform_int_distribution<int> dist_rx(region.lx, region.ux);
                    std::uniform_int_distribution<int> dist_ry(region.ly, region.uy);
                    Bx = dist_rx(rng);
                    By = dist_ry(rng);
                } else {
                    int left   = std::max(0, Ax - search_rx);
                    int right  = std::min(R - 1, Ax + search_rx);
                    int bottom = std::max(0, Ay - search_ry);
                    int top    = std::min(C - 1, Ay + search_ry);
                    std::uniform_int_distribution<int> dist_rx(left, right);
                    std::uniform_int_distribution<int> dist_ry(bottom, top);
                    Bx = dist_rx(rng);
                    By = dist_ry(rng);
                }
            } while (Bx == Ax && By == Ay);

            string swap_target = array2D[Bx][By];
            string nameA = A.name;
            string temp_name;

            // 做 swap（可能跟 none 換）
            if (swap_target != "none")
            {
                int idxB = name2idx[swap_target];
                auto &B = logic_blocks_table[idxB];
                temp_name = B.name;
                A.x = Bx; A.y = By;
                B.x = Ax; B.y = Ay;
                array2D[Ax][Ay] = swap_target;
                array2D[Bx][By] = nameA;
            }
            else
            {
                temp_name = "none";
                A.x = Bx; A.y = By;
                array2D[Bx][By] = nameA;
                array2D[Ax][Ay] = "none";
            }

            change_info chg;
            chg.nameA = nameA; chg.Ax = Ax; chg.Ay = Ay;
            chg.nameB = temp_name; chg.Bx = Bx; chg.By = By;

            // 收集 affected nets
            std::vector<int> affected_nets;
            affected_nets.reserve(32);

            auto collect_nets = [&](const std::string &blk_name) {
                if (blk_name == "none") return;
                auto it_blk = name2idx.find(blk_name);
                if (it_blk == name2idx.end()) return;
                const auto &blk = logic_blocks_table[it_blk->second];
                for (int nid : blk.connect_nets)
                {
                    if (net_mark[nid] != net_tag)
                    {
                        net_mark[nid] = net_tag;
                        affected_nets.push_back(nid);
                    }
                }
            };
            collect_nets(chg.nameA);
            collect_nets(chg.nameB);

            if (affected_nets.empty())
            {
                // 這個 move 沒有任何 net 受影響 → 沒差，直接略過
                undo_move(logic_blocks_table, name2idx, array2D, chg);
                continue;
            }

            struct NetUpdate
            {
                int net_idx;
                double old_hpwl;
                rec old_box;
                double new_hpwl;
                rec new_box;
            };
            std::vector<NetUpdate> updates;
            updates.reserve(affected_nets.size());

            double trial_hpwl = total_hpwl;
            for (int nid : affected_nets)
            {
                const net_info &net = net_table[nid];
                double old_hp   = net_cost[nid].hpwl;
                rec    old_rect = net_cost[nid].box;
                rec    new_rect = calculate_bounding_box(net, logic_blocks_table, pin_table);
                double new_hp   = (new_rect.xmax - new_rect.xmin) + (new_rect.ymax - new_rect.ymin);
                trial_hpwl += (new_hp - old_hp);
                updates.push_back({nid, old_hp, old_rect, new_hp, new_rect});
            }

            // 對 U / sumU / sumU2 套用 old -> -1, new -> +1
            for (const auto &u : updates)
            {
                if (u.old_hpwl > 0.0) apply_rect_delta_on_U(u.old_box, -1, U, sumU, sumU2);
                if (u.new_hpwl > 0.0) apply_rect_delta_on_U(u.new_box, +1, U, sumU, sumU2);
            }

            double new_CC;
            if (sumU == 0 || Nsite == 0) new_CC = 1.0;
            else new_CC = (double)Nsite * (double)sumU2 / ((double)sumU * (double)sumU);

            double new_cost = compute_cost(trial_hpwl, new_CC, lambda);
            double delta_cost = new_cost - cur_cost;

            bool accept = false;
            if (delta_cost <= 0.0) accept = true;
            else if (dist01(rng) < exp(-delta_cost / T)) accept = true;

            if (accept)
            {
                cur_cost   = new_cost;
                total_hpwl = trial_hpwl;
                CC         = new_CC;

                for (const auto &u : updates)
                {
                    net_cost[u.net_idx].hpwl = u.new_hpwl;
                    net_cost[u.net_idx].box  = u.new_box;
                }

                if (cur_cost < best_cost)
                {
                    best_cost = cur_cost;
                    best_hpwl = total_hpwl;
                    best_CC   = CC;
                    best_logic_blocks = logic_blocks_table;
                    best_array2D      = array2D;
                }
                accepted_moves++;
            }
            else
            {
                // U / sumU / sumU2 revert：new -1, old +1
                for (const auto &u : updates)
                {
                    if (u.new_hpwl > 0.0) apply_rect_delta_on_U(u.new_box, -1, U, sumU, sumU2);
                    if (u.old_hpwl > 0.0) apply_rect_delta_on_U(u.old_box, +1, U, sumU, sumU2);
                }
                undo_move(logic_blocks_table, name2idx, array2D, chg);
            }
        }

        // 動態調整 region prob
        double acceptance_rate = (double)accepted_moves / (double)iter_per_T;

        if (acceptance_rate > 0.85) {
            current_region_prob -= prob_step_down;
        } else if (acceptance_rate > 0.6) {
            current_region_prob += prob_step_up;
        } else if (acceptance_rate > 0.15) {
            current_region_prob += prob_step_up;
        } else {
            current_region_prob -= prob_step_down;
        }
        if (current_region_prob < min_region_prob) current_region_prob = min_region_prob;
        if (current_region_prob > max_region_prob) current_region_prob = max_region_prob;

        // 溫度更新
        double alpha;
        if (time_ratio > 0.9)       alpha = 0.5;
        else if (time_ratio > 0.6)  alpha = 0.82;
        else {
            if      (acceptance_rate > 0.9)  alpha = 0.7;
            else if (acceptance_rate > 0.7)  alpha = 0.85;
            else if (acceptance_rate < 0.15) alpha = 0.95;
            else                             alpha = 0.9;
        }
        T *= alpha;

        // moves / T 動態
        if (acceptance_rate > 0.9)
        {
            int nm = (int)(iter_per_T * 0.9);
            iter_per_T = max(nm, 4000);
        }
        else if (acceptance_rate < 0.2)
        {
            int nm = (int)(iter_per_T * 1.1);
            iter_per_T = min(nm, 60000);
        }

        round_cnt++;
        if (acceptance_rate < 0.02) low_acc_rounds++;
        else low_acc_rounds = 0;
        if (low_acc_rounds >= 6)
        {
            cerr << "[SA] Early stop due to very low acceptance.\n";
            break;
        }

        cerr << "Round=" << round_cnt
             << " T=" << T
             << " Acc=" << acceptance_rate
             << " HPWL=" << best_hpwl
             << " CC=" << best_CC
             << " Lam=" << lambda
             << " Time=" << time_ratio << endl;
    }

    logic_blocks_table = best_logic_blocks;
    array2D = best_array2D;
    total_hpwl = best_hpwl;
    CC = best_CC;

    // 重新算一次，確保狀態一致
    compute_hpwl_cong(logic_blocks_table, pin_table, net_table, U, net_cost, total_hpwl, CC);
}

void parse_input(ifstream &parse,
                 vector<logic_info> &logic_blocks_table,
                 vector<pin_info> &pin_table,
                 vector<net_info> &net_table,
                 unordered_map<string, int> &name2idx,
                 unordered_map<string, int> &pin2idx,
                 unordered_map<string, int> &net2idx)
{
    string line;
    getline(parse, line);
    istringstream info(line);
    int logic_block_R, logic_block_C, num_logic_blocks, num_io_pins, num_nets;
    info >> logic_block_R >> logic_block_C >> num_logic_blocks >> num_io_pins >> num_nets;
    R = logic_block_R;
    C = logic_block_C;
    for (int i = 0; i < num_logic_blocks; i++)
    {
        string name;
        getline(parse, line);
        istringstream logic_inf(line);
        logic_inf >> name;
        logic_info inst;
        inst.name = name; inst.x = 0; inst.y = 0;
        logic_blocks_table.push_back(inst);
        name2idx.insert({name, i});
    }
    for (int i = 0; i < num_io_pins; i++)
    {
        string name;
        double x, y;
        getline(parse, line);
        istringstream pin_inf(line);
        pin_inf >> name >> x >> y;
        pin_info inst;
        inst.name = name; inst.x = x; inst.y = y;
        pin_table.push_back(inst);
        pin2idx.insert({name, i});
    }
    for (int i = 0; i < num_nets; i++)
    {
        string name;
        int degree;
        getline(parse, line);
        istringstream net_inf(line);
        net_inf >> name >> degree;
        net_info inst;
        inst.name = name; inst.degree = degree;

        for (int j = 0; j < degree; j++)
        {
            string item;
            net_inf >> item;
            inst.connections.push_back(item);
            NetTerm term;
            auto it_blk = name2idx.find(item);
            if (it_blk != name2idx.end())
            {
                term.is_pin = false;
                term.idx = it_blk->second;
                logic_blocks_table[it_blk->second].connect_nets.push_back(i);
            }
            else
            {
                auto it_pin = pin2idx.find(item);
                if (it_pin != pin2idx.end())
                {
                    term.is_pin = true;
                    term.idx = it_pin->second;
                }
                else
                {
                    term.is_pin = true;
                    term.idx = -1;
                }
            }
            inst.terms.push_back(term);
        }
        net_table.push_back(inst);
        net2idx.insert({name, i});
    }
}

void write_output(const string &out_file, const vector<logic_info> &logic_blocks_table)
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

int main(int argc, char *argv[])
{
    auto start = Clock::now();

    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " <input_file> <out_file>" << endl;
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

    double total_hpwl, CC;

    parse_input(parse_in, logic_blocks_table, pin_table, net_table,
                name2idx, pin2idx, net2idx);

    vector<vector<string>> array2D(R, vector<string>(C, "none"));
    vector<vector<int>> U(R, vector<int>(C, 0));
    vector<net_hpwl_cong> net_cost(net_table.size());

    initial_place(array2D, R, C, logic_blocks_table);

    compute_hpwl_cong(logic_blocks_table, pin_table, net_table,
                      U, net_cost, total_hpwl, CC);

    SA(logic_blocks_table, array2D,
       pin_table, net_table,
       name2idx, pin2idx,
       net_cost, U,
       total_hpwl, CC);

    cout << "[FINAL] HPWL=" << total_hpwl
         << " CC=" << CC << endl;

    write_output(string(argv[2]), logic_blocks_table);

    auto end = Clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Execution time: " << duration.count() << " milliseconds" << endl;

    return 0;
}
